@doc raw"""
    get_crystal_coords!(crystal :: AbstractMatrix{T}, cartesian :: AbstractMatrix{T}, cell :: AbstractMatrix{T}; buffer = default_buffer())

Convert cartesian coordinates into crystal ones. 
Optionally, a Bumper buffer can be provided to avoid memory allocations.

"""
function get_crystal_coords!(crystal :: AbstractMatrix{T}, cartesian :: AbstractMatrix{T}, cell :: AbstractMatrix{T}; buffer = default_buffer()) where T
    n_atoms = size(cartesian, 2)
    @no_escape buffer begin
        metric_tensor = @alloc(T, 3, 3)
        inv_metric_tensor = @alloc(T, 3, 3)
        tmp_vectors = @alloc(T, 3, n_atoms)

        mul!(metric_tensor, cell', cell)
        inv_metric_tensor .= inv(metric_tensor)

        mul!(tmp_vectors, cell', cartesian)
        mul!(crystal, inv_metric_tensor, tmp_vectors)
    end
end

function get_cartesian_coords!(cartesian :: AbstractMatrix{T}, crystal :: AbstractMatrix{T}, cell :: AbstractMatrix{T}) where T
    mul!(cartesian, cell, crystal)
end
