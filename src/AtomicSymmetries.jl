module AtomicSymmetries

using Spglib
using LinearAlgebra
using Bumper

# Bumper.allow_ptr_array_to_escape() = true

include("symmetries_core.jl")
include("spglib_init.jl")
include("generators.jl")
include("asr.jl")
include("filter_symmetries.jl")
include("crystal.jl")

export get_symmetry_group_from_spglib,
       get_nsymmetries,
       ASRVectorConstraint!,
       ASRMatrixConstraint!,
       Symmetries,
       get_empty_symmetry_group,
       complete_symmetry_group!,
       filter_invariant_symmetries!,
       get_crystal_coords!,
       get_cartesian_coords!,
       symmetrize_positions!

end # module AtomicSymmetries
