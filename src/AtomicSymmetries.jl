module AtomicSymmetries

using Spglib
using LinearAlgebra
using SparseArrays
using Bumper

# Bumper.allow_ptr_array_to_escape() = true

include("symmetries_core.jl")
include("spglib_init.jl")
include("generators.jl")
include("asr.jl")
include("filter_symmetries.jl")
include("crystal.jl")
include("sparsify.jl")

export get_symmetry_group_from_spglib,
       get_nsymmetries,
       symmetrize_fc!,
       symmetrize_vector!,
       symmetrize_positions!,
       ASRConstraint!,
       Symmetries,
       get_empty_symmetry_group,
       complete_symmetry_group!,
       filter_invariant_symmetries!,
       get_crystal_coords!,
       get_cartesian_coords!,
       to_primitive_cell_cryst!,
       to_primitive_cell_cart!,
       apply_sparse_symmetry


end # module AtomicSymmetries
