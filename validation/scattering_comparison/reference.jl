using LinearAlgebra: I, Diagonal, cond, cross, dot, norm
using StaticArrays: SVector

# This reference uses surface quadrature and explicit segment visibility. It
# deliberately does not call ArchimedLight's projections, ray stacks, or graph
# construction. Each object has a uniform scattered radiance on its two faces.

function _reference_rectangle(object, subdivisions::Int)
    vertices = object.vertices
    length(vertices) == 4 || throw(ArgumentError("Reference objects must have four ordered rectangle vertices"))
    origin = SVector{3,Float64}(vertices[1])
    u = SVector{3,Float64}(vertices[2]) - origin
    v = SVector{3,Float64}(vertices[4]) - origin
    normal0 = cross(u, v)
    geometric_area = norm(normal0)
    geometric_area > 0.0 || throw(ArgumentError("Reference rectangles must have positive area"))
    area = Float64(object.area)
    isfinite(area) && area > 0.0 || throw(ArgumentError("Reference object area must be finite and positive"))
    isapprox(area, geometric_area; rtol=1e-6, atol=1e-12) ||
        throw(ArgumentError("Object area does not match its rectangle geometry"))
    isapprox(SVector{3,Float64}(vertices[3]), origin + u + v; rtol=1e-6, atol=1e-10) ||
        throw(ArgumentError("Reference vertices must form a parallelogram in perimeter order"))
    uu, uv, vv = dot(u, u), dot(u, v), dot(v, v)
    abs(uv) <= 1e-6 * sqrt(uu * vv) ||
        throw(ArgumentError("Reference objects must be rectangles"))
    inverse_gram_determinant = inv(uu * vv - uv * uv)
    points = Vector{SVector{3,Float64}}(undef, subdivisions^2)
    k = 1
    for j in 1:subdivisions, i in 1:subdivisions
        points[k] = origin + ((i - 0.5) / subdivisions) * u + ((j - 0.5) / subdivisions) * v
        k += 1
    end
    return (
        origin=origin,
        u=u,
        v=v,
        normal=normal0 / geometric_area,
        uu=uu,
        uv=uv,
        vv=vv,
        inverse_gram_determinant=inverse_gram_determinant,
        area=area,
        points=points,
        point_area=area / length(points),
    )
end

@inline function _reference_segment_intersects_rectangle(origin, displacement, rectangle)
    denominator = dot(rectangle.normal, displacement)
    # A segment parallel to a blocker cannot cross its plane. Source and target
    # rectangles are excluded by the caller; coplanar edge contacts have zero
    # area in the surface integral.
    abs(denominator) > 1e-12 * norm(displacement) || return false
    t = dot(rectangle.normal, rectangle.origin - origin) / denominator
    1e-10 < t < 1.0 - 1e-10 || return false
    offset = origin + t * displacement - rectangle.origin
    du, dv = dot(offset, rectangle.u), dot(offset, rectangle.v)
    a = (du * rectangle.vv - dv * rectangle.uv) * rectangle.inverse_gram_determinant
    b = (dv * rectangle.uu - du * rectangle.uv) * rectangle.inverse_gram_determinant
    return -1e-10 <= a <= 1.0 + 1e-10 && -1e-10 <= b <= 1.0 + 1e-10
end

@inline function _reference_segment_visible(origin, displacement, rectangles, source::Int, target::Int)
    for blocker in eachindex(rectangles)
        (blocker == source || blocker == target) && continue
        _reference_segment_intersects_rectangle(origin, displacement, rectangles[blocker]) && return false
    end
    return true
end

"""
    reference_transfer(objects; subdivisions=24)

Compute an independent Lambertian transfer matrix `F[to, from]` for nonintersecting
rectangular surfaces in open space. Every object provides `vertices` (four ordered
`SVector{3,Float64}` vertices), `area`, `scatter`, and `label`. The latter two fields
are not used in the geometric calculation. All listed objects block visibility.

Composite midpoint quadrature uses `subdivisions^2` points per rectangle. Each
source's outgoing power is split equally between its two Lambertian faces, so
the double-area integral has denominator `2π * area_from * distance^4` when its
numerator uses normal dot products with the displacement vector. It includes
escape by leaving each column's missing fraction unassigned.

Each unordered pair is integrated only once, enforcing area reciprocity. Compare
`subdivisions=24` and `48` to quantify reference quadrature error; do not infer
convergence from energy closure alone. Partial occlusion can converge more slowly
because the visibility indicator is discontinuous.
"""
function reference_transfer(objects; subdivisions::Integer=24)
    subdivisions > 0 || throw(ArgumentError("subdivisions must be positive"))
    rectangles = [_reference_rectangle(object, Int(subdivisions)) for object in objects]
    n = length(rectangles)
    F = zeros(Float64, n, n)
    for source in 1:n, target in (source + 1):n
        a, b = rectangles[source], rectangles[target]
        kernel_sum = 0.0
        for source_point in a.points, target_point in b.points
            displacement = target_point - source_point
            distance2 = dot(displacement, displacement)
            distance2 > 0.0 || continue
            numerator = abs(dot(a.normal, displacement)) * abs(dot(b.normal, displacement))
            numerator > 0.0 || continue
            _reference_segment_visible(source_point, displacement, rectangles, source, target) || continue
            kernel_sum += numerator / (distance2 * distance2)
        end
        exchange_area = kernel_sum * a.point_area * b.point_area / (2pi)
        F[target, source] = exchange_area / a.area
        F[source, target] = exchange_area / b.area
    end
    return F
end

"""
    reference_parallel_square_fraction(; side=1.0, separation=1.0)

Analytical two-face transfer fraction for equal, aligned parallel squares.
At unit side and separation it is `0.199824895698 / 2` (approximately).
This provides an independent check on the numerical surface quadrature.

The one-face expression is Eq. 71 in Salazar et al. (2016), NASA CHAR verification:
https://ntrs.nasa.gov/citations/20160006076
"""
function reference_parallel_square_fraction(; side::Real=1.0, separation::Real=1.0)
    side > 0 && separation > 0 || throw(ArgumentError("side and separation must be positive"))
    x = Float64(side / separation)
    s = sqrt(1 + x^2)
    return (0.5 * log((1 + x^2)^2 / (1 + 2x^2)) +
            2x * s * atan(x / s) - 2x * atan(x)) / (pi * x^2)
end

"""
    reference_parallel_square_directional_fraction(sectors; side=1.0, separation=1.0)

Exact projected-overlap reference for the candidate's angular discretization of
two equal, aligned horizontal squares. Each sector supplies `direction` and a
nonnegative solid-angle `weight`; all directions must point upward. The projected
overlap is `(side - abs(separation*s_x/s_z))₊` times its y counterpart.

The result uses the same discrete source normalization and equal two-face split
as the candidate, but has no pixels or triangle rasterization. Compare it with
`reference_parallel_square_fraction()` to measure angular error, and with the
candidate to measure the additional raster error. Multiply by source power and
its scattering coefficient to obtain the black receiver's power.
"""
function reference_parallel_square_directional_fraction(sectors; side::Real=1.0, separation::Real=1.0)
    side > 0 && separation > 0 || throw(ArgumentError("side and separation must be positive"))
    width, distance = Float64(side), Float64(separation)
    weighted_overlap = 0.0
    weighted_source_area = 0.0
    for sector in sectors
        sx, sy, sz = Float64.(sector.direction)
        weight = Float64(sector.weight)
        all(isfinite, (sx, sy, sz, weight)) && sz > 0 && weight >= 0 ||
            throw(ArgumentError("Sectors must have finite upward directions and nonnegative weights"))
        overlap = max(0.0, width - abs(distance * sx / sz)) *
                  max(0.0, width - abs(distance * sy / sz))
        measure = weight * sz
        weighted_overlap += measure * overlap
        weighted_source_area += measure * width^2
    end
    weighted_source_area > 0 || throw(ArgumentError("Sector weights must have positive total measure"))
    return weighted_overlap / (2 * weighted_source_area)
end

"""
    solve_reference(F, b, c)

Solve repeated exchanges for geometric transfer matrix `F[to, from]`, initial
intercepted powers `b`, and scattering coefficients `c`. The assumption of uniform
redistribution within each object is retained. Thus `incident = b + F*Diagonal(c)*incident`.

`scattered` is the added received power (`incident - b`), while `order1` is the
received power from the first exchange alone. `escaped` is a scene total; neither
cumulative incident power nor virtual observations belong in an energy closure.
No column normalization or clipping is applied to hide quadrature error.
"""
function solve_reference(F::AbstractMatrix{<:Real}, b::AbstractVector{<:Real}, c::AbstractVector{<:Real})
    n = length(b)
    size(F) == (n, n) && length(c) == n || throw(DimensionMismatch("F, b, and c must describe the same objects"))
    all(isfinite, F) && all(x -> x >= 0, F) || throw(ArgumentError("F must be finite and nonnegative"))
    all(x -> isfinite(x) && x >= 0, b) || throw(ArgumentError("Initial powers must be finite and nonnegative"))
    all(x -> isfinite(x) && 0 <= x <= 1, c) || throw(ArgumentError("Scattering coefficients must be in [0, 1]"))
    initial = Float64.(b)
    coefficients = Float64.(c)
    transfer = Matrix{Float64}(F) * Diagonal(coefficients)
    system = Matrix{Float64}(I, n, n) - transfer
    incident = system \ initial
    absorbed = (1 .- coefficients) .* incident
    scattered = incident - initial
    order1 = transfer * initial
    escaped_per_node = (1 .- vec(sum(F; dims=1))) .* coefficients .* incident
    escaped = sum(escaped_per_node)
    return (
        initial=initial,
        incident=incident,
        absorbed=absorbed,
        scattered=scattered,
        order1=order1,
        escaped=escaped,
        escaped_per_node=escaped_per_node,
        balance_error=sum(absorbed) + escaped - sum(initial),
        residual_norm=norm(system * incident - initial),
        condition_number=cond(system),
        solver=:direct_linear_system,
    )
end
