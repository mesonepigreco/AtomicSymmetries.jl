# This files init the symmetries from the spglib library
# Initialize a symmetry group starting from spglib

@doc raw"""
    get_symmetry_group_from_spglib(positions :: AbstractMatrix, cell :: AbstractMatrix, types :: Array{Int}; symprec :: Float64 = 1e-6)

Build a symmetry group from the spglib library.

# Arguments

- `positions::AbstractMatrix`: The atomic positions in the cell (crystallographic coordinates), with shape `(3, N)`.
- `cell::AbstractMatrix`: The cell matrix with shape `(3, 3)`.
- `types::Array{Int}`: The atomic types.
- `symprec::Float64`: The symmetry precision.
"""
function get_symmetry_group_from_spglib(positions::AbstractMatrix{<: Real}, cell::AbstractMatrix{<:Real}, types::Vector{<:Int}; symprec::Float64 = 1e-6, type::Type = Float64) :: Symmetries
    nat = size(positions, 2)
    ndim = size(positions, 1)

    @assert ndim == 3 "Only 3D systems are supported for SPGLIB symmetry groups."

    # Build the SPGLIB cell object
    cell = Cell(cell, positions, types)
    
    # Get the symmetry operations
    R, T = Spglib.get_symmetry(cell, symprec)

    # Create the new symmetry group
    sym_group = get_identity_symmetry_group(type; dims=ndim, n_atoms=nat)

    # Add the symmetry operations excluding the first identity
    for i in 2:length(R)
        irt = zeros(Int, nat)
        get_irt!(irt, positions, R[i], T[i])

        # Add the symmetry operation
        add_symmetry!(sym_group, R[i], update=false; irt=irt)
    end

    update_symmetry_functions!(sym_group)

    return sym_group
end


