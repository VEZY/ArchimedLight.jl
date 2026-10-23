# Shared geometry for the bundled inputs and the interactive documentation.
# Coordinates are metres; the original scene uses unit voxels on [0,8]×[0,8]×[0,9].
function _voxel_docs_elevation(x, y)
    broad_slope = 0.055x + 0.025y
    central_mound = 0.52exp(-((x - 4.8)^2 + 0.7(y - 4.2)^2) / 5.5)
    shallow_hollow = 0.28exp(-((x - 2.0)^2 + (y - 6.0)^2) / 2.0)
    ripples = 0.13sinpi(x / 4) * cospi(y / 4)
    # Keep the full terrain below the first trunk voxel, so the base is visible.
    return 0.60 * (
        0.42 + broad_slope + central_mound - shallow_hollow + ripples
    )
end

function _voxel_docs_tree_trunk(x, y, z)
    return (x - 4.35)^2 + (y - 4.65)^2 <= 0.42^2 &&
           z > _voxel_docs_elevation(x, y) + 0.10 && z <= 4.60
end

function _voxel_docs_tree_pad(nx, ny, nz)
    pad = zeros(nx, ny, nz)
    tree_x, tree_y = 4.35, 4.65
    for k in 1:nz, j in 1:ny, i in 1:nx
        x, y, z = i - 0.5, j - 0.5, k - 0.5
        trunk = _voxel_docs_tree_trunk(x, y, z)

        # An ellipsoid with a deterministic, slightly ruffled boundary gives
        # the crown an organic outline without introducing randomness.
        crown_metric = ((x - tree_x) / 2.55)^2 +
                       ((y - tree_y) / 2.25)^2 +
                       ((z - 6.30) / 2.40)^2
        crown_irregularity = 0.10sin(1.7x + 0.7y + 0.4z) +
                             0.07cos(1.2y - 1.1z)
        boundary_hole = crown_metric > 0.62 && mod(5i + 3j + 7k, 17) == 0
        crown = crown_metric + crown_irregularity <= 1.0 && !boundary_hole

        crown_density = crown ?
            1.02 + 0.68clamp(1.0 - crown_metric, 0.0, 1.0) : 0.0
        trunk_density = trunk ? 0.34 : 0.0
        # Give the trunk precedence where it enters the crown.
        pad[i, j, k] = trunk ? trunk_density : crown_density
    end
    return pad
end
