# Build-time overlay only: WasmGC vectors do not expose native capacity. Ignoring
# a size hint preserves values, size, ordering, and all transport calculations.
# WasmTarget 0.5.3 otherwise inlines Base's unsupported memmove for these buffers.
using Base.Experimental: @overlay

const VoxelWasmHintVector = Union{
    Vector{VoxelWasmKernel.AL.VoxelRaySegment{Float64}},
    Vector{NTuple{2,Float64}},
}

@overlay WasmTarget.WASM_METHOD_TABLE Base.sizehint!(v::VoxelWasmHintVector, n::Integer) = v

# The generic array-field/new path currently chooses the vector registry even
# for Array{T,N}. Reuse its existing N-dimensional registry so :size is NTuple{N}.
function WasmTarget.register_vector_type!(module_::WasmTarget.WasmModule,
                                         registry::WasmTarget.TypeRegistry,
                                         type::Type{Array{T,N}}) where {T,N}
    N > 1 && return WasmTarget.register_matrix_type!(module_, registry, type)
    return invoke(WasmTarget.register_vector_type!,
                  Tuple{WasmTarget.WasmModule,WasmTarget.TypeRegistry,Type},
                  module_, registry, type)
end
