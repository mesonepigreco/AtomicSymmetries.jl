module AtomicSymmetries

using Spglib
using LinearAlgebra


include("symmetries_core.jl")
include("spglib_init.jl")

export get_symmetry_group_from_spglib,
       get_nsymmetries

end # module AtomicSymmetries
