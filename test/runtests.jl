using Test
using AtomicSymmetries

input("test_bcc.jl")
@testset "BCC" test_bcc()

input("test_r3m.jl")
@testset "R3M" test_r3m()
