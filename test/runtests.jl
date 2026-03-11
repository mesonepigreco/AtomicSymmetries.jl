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
    test_filter_fourier()
end

include("test_coset_decomposition.jl")
@testset "Coset decomposition" begin
    test_coset_decomposition()
    test_coset_decomposition_qspace()
end

@testset "Cart to Crystal" begin
    include("test_crystal_cart.jl")
    test_crystal_to_cart()
end

@testset "ASR and Vector symmetrization (non-orthogonal cell)" begin
    include("test_symmetrize_vector_cart.jl")
    test_symmetrize_vector()
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


@testset "generators of a non-orthogonal cell" begin
    include("test_pbte_generators.jl")
    test_find_generators_pbte()
end

@testset "sparsify" begin
    include("test_sparsify.jl")
    test_symmetry_sparsification()
end


@testset "Fourier symmetries" begin
    include("test_symmetrize_qspace.jl")
    test_fourier_matrix()
    test_symmetrize_q_space()

    # A test on a complex noncubic cell
    include("test_symmetrize_cartesian_qspace.jl")
    test_symmetrize_cartesian_qspace()

    # A test with fractional translations
    include("test_fractional_symmetries.jl")
    test_fractional_symmetries_qspace()
end

include("test_general_interface.jl")
@testset "General interface" begin
    test_rotate_vector_real()
    test_rotate_matrix_real()
    test_rotate_dynamical_matrix_real()
    test_rotate_vector_qspace()
    test_rotate_dynamical_matrix_qspace()
    test_rotate_centroid_real()
    test_rotate_centroid_identity()
end

include("test_efficient_generators.jl")
include("test_fcc_supercell_generators.jl")
@testset "Efficient generators" begin
    test_index_helpers()
    test_apply_symmetry_tensor()
    test_rank2_consistency()
    test_rank3_generators()
    test_contraction()
    test_compact_reconstruction()
    test_fast_vs_standard_generators()
    test_fcc_444_rotated_generators()
end
