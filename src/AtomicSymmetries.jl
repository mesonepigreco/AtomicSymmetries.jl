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
include("fourier_transform.jl")
include("symmetrize_qspace.jl")

export get_symmetry_group_from_spglib,
       get_nsymmetries,
       symmetrize_fc!,
       symmetrize_vector!,
       symmetrize_positions!,
       ASRConstraint!,
       translation_mask!,
       Symmetries,
       get_empty_symmetry_group,
       complete_symmetry_group!,
       filter_invariant_symmetries!,
       get_crystal_coords!,
       get_cartesian_coords!,
       to_primitive_cell_cryst!,
       to_primitive_cell_cart!,
       apply_sparse_symmetry,
       apply_sparse_symmetry!,
       vector_r2q!, vector_q2r!,
       matrix_r2q!, matrix_q2r!,
       SymmetriesQSpace, 
       apply_symmetry_vectorq!, apply_symmetry_matrixq!,
       get_translations,
       apply_translations!,
       symmetrize_vector_q!, symmetrize_matrix_q!,
       symmetrize_vector_cartesian_q!,
       symmetrize_matrix_cartesian_q!




end # module AtomicSymmetries
