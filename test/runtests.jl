using Test
using AtomicSymmetries
using PyCall

# Include generic files
include("define_cell.jl")

include("test_bcc.jl")
@testset "BCC" test_bcc()

include("test_r3m.jl")
@testset "R3M" test_r3m()

@testset "Spglib tests" begin
    include("test_spglib.jl")
    test_load_spglib()
    test_symmetries_supercell()
    test_gold_unit_cell_spglib()
    test_gold_supercell_spglib()

    # Check that the python spglib module works
    try
        pyimport("spglib")
        test_python_spglib()
    catch
        @warn "The python spglib module is not installed. Skipping the test."
    end
end

include("test_filter.jl")
@testset "Filter" begin 
    test_filter()
    test_filter_full_vector()
end

@testset "Cart to Crystal" begin
    include("test_crystal_cart.jl")
    test_crystal_to_cart()
end

@testset "exact fc symmetrization" begin
    include("test_fc_symmetrization.jl")
    test_fc_sym_primitive_cell()
    test_fc_sym_pbte_uc() # This test symmetries without IRT
    test_fc_sym_pbte_444() # This test also pure translations
end

@testset "translations and irt and ASR" begin
    include("test_pbte_symmetrization.jl")
    test_pbte_supercell_symmetrization()
    test_asr_impose()
end


