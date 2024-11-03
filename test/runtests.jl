using Test
using AtomicSymmetries

# Include generic files
include("define_cell.jl")

include("test_bcc.jl")
@testset "BCC" test_bcc()

include("test_r3m.jl")
@testset "R3M" test_r3m()

include("test_spglib.jl")
@testset "Spglib simple" test_load_spglib()
@testset "Spglib supercell" test_symmetries_supercell()

include("test_filter.jl")
@testset "Filter" test_filter()
