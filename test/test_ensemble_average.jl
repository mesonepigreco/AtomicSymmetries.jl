using AtomicSymmetries
using Test
using LinearAlgebra
using Random

"""
Test rank-2 ensemble average: compute Φ_ab = -⟨v_a f_b⟩ from random ensemble data,
compare coefficients from ensemble vs from the explicit full tensor.
"""
function test_ensemble_rank2()
    # PbTe unit cell (2 atoms)
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'
    atomic_numbers = [1, 2]

    symmetry_group = get_symmetry_group_from_spglib(positions, cell, atomic_numbers)
    dim = 3
    nat = 2
    n_modes = dim * nat
    n_configs = 200

    generators = get_tensor_generators_fast(symmetry_group, cell; rank=2)
    @test length(generators) > 0

    # Generate random ensemble data
    Random.seed!(42)
    v = randn(Float64, dim, nat, n_configs)
    f = randn(Float64, dim, nat, n_configs)

    # Build the full rank-2 tensor: Φ[a,b] = -(1/Nc) Σ_I v_flat[a,I] * f_flat[b,I]
    tensor = zeros(Float64, n_modes, n_modes)
    for I in 1:n_configs
        v_flat = reshape(view(v, :, :, I), n_modes)
        f_flat = reshape(view(f, :, :, I), n_modes)
        for a in 1:n_modes, b in 1:n_modes
            tensor[a, b] -= v_flat[a] * f_flat[b]
        end
    end
    tensor ./= n_configs

    # Reference: coefficients from full tensor
    coeffs_ref = zeros(Float64, length(generators))
    get_coefficients_from_tensor!(coeffs_ref, tensor, generators, cell)

    # Check that the generators actually satisfy the permutation symmetry.
    coeffs_ref_sym = zeros(Float64, length(generators))
    t2 = tensor + tensor'
    t2 ./= 2.0
    get_coefficients_from_tensor!(coeffs_ref_sym, t2, generators, cell)

    # Explicitly check the permutation invariance of the generators: coefficients should be the same for the symmetrized tensor
    @test coeffs_ref ≈ coeffs_ref_sym atol=1e-10

    # Ensemble: coefficients directly from v, f
    coeffs_ens = zeros(Float64, length(generators))
    get_coefficients_from_ensemble!(coeffs_ens, v, f, generators)

    @test coeffs_ens ≈ coeffs_ref atol=1e-10
end


"""
Test rank-3 ensemble average: compute Φ_abc = -⟨v_a v_b f_c⟩.
Uses zinc-blende (F-43m, non-centrosymmetric) to have non-trivial rank-3 generators.
"""
function test_ensemble_rank3()
    positions, cell, types = get_zincblende()

    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)
    dim = 3
    nat = 2
    n_modes = dim * nat
    n_configs = 500

    generators = get_tensor_generators_fast(symmetry_group, cell; rank=3)
    @test length(generators) > 0

    Random.seed!(123)
    v = randn(Float64, dim, nat, n_configs)
    f = randn(Float64, dim, nat, n_configs)

    # Build full rank-3 tensor: Φ[a,b,c] = -(1/Nc) Σ_I v_flat[a]*v_flat[b]*f_flat[c]
    tensor = zeros(Float64, n_modes, n_modes, n_modes)
    for I in 1:n_configs
        v_flat = reshape(view(v, :, :, I), n_modes)
        f_flat = reshape(view(f, :, :, I), n_modes)
        for a in 1:n_modes, b in 1:n_modes, c in 1:n_modes
            tensor[a, b, c] -= v_flat[a] * v_flat[b] * f_flat[c]
        end
    end
    tensor ./= n_configs

    coeffs_ref = zeros(Float64, length(generators))
    get_coefficients_from_tensor!(coeffs_ref, tensor, generators, cell)

    coeffs_ens = zeros(Float64, length(generators))
    get_coefficients_from_ensemble!(coeffs_ens, v, f, generators)

    @test coeffs_ens ≈ coeffs_ref atol=1e-10
end


"""
Test rank-4 ensemble average: compute Φ_abcd = -⟨v_a v_b v_c f_d⟩.
"""
function test_ensemble_rank4()
    positions, cell, types = get_zincblende()

    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)
    dim = 3
    nat = 2
    n_modes = dim * nat
    n_configs = 500

    generators = get_tensor_generators_fast(symmetry_group, cell; rank=4)
    @test length(generators) > 0

    Random.seed!(456)
    v = randn(Float64, dim, nat, n_configs)
    f = randn(Float64, dim, nat, n_configs)

    # Build full rank-4 tensor: Φ[a,b,c,d] = -(1/Nc) Σ_I v[a]*v[b]*v[c]*f[d]
    tensor = zeros(Float64, n_modes, n_modes, n_modes, n_modes)
    for I in 1:n_configs
        v_flat = reshape(view(v, :, :, I), n_modes)
        f_flat = reshape(view(f, :, :, I), n_modes)
        for a in 1:n_modes, b in 1:n_modes, c in 1:n_modes, d in 1:n_modes
            tensor[a, b, c, d] -= v_flat[a] * v_flat[b] * v_flat[c] * f_flat[d]
        end
    end
    tensor ./= n_configs

    coeffs_ref = zeros(Float64, length(generators))
    get_coefficients_from_tensor!(coeffs_ref, tensor, generators, cell)

    coeffs_ens = zeros(Float64, length(generators))
    get_coefficients_from_ensemble!(coeffs_ens, v, f, generators)

    @test coeffs_ens ≈ coeffs_ref atol=1e-10
end


"""
Test that the streaming API (accumulate_ensemble_config!) gives the same result
as the batch API (get_coefficients_from_ensemble!).
"""
function test_streaming_api()
    positions, cell, types = get_zincblende()

    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)
    dim = 3
    nat = 2
    n_configs = 500

    generators = get_tensor_generators_fast(symmetry_group, cell; rank=3)
    @test length(generators) > 0

    Random.seed!(789)
    v = randn(Float64, dim, nat, n_configs)
    f = randn(Float64, dim, nat, n_configs)

    # Batch
    coeffs_batch = zeros(Float64, length(generators))
    get_coefficients_from_ensemble!(coeffs_batch, v, f, generators)

    # Streaming
    coeffs_stream = zeros(Float64, length(generators))
    for I in 1:n_configs
        v_config = @view v[:, :, I]
        f_config = @view f[:, :, I]
        accumulate_ensemble_config!(coeffs_stream, v_config, f_config, generators)
    end
    coeffs_stream .*= -1.0 / n_configs

    @test coeffs_stream ≈ coeffs_batch atol=1e-10
end


if abspath(PROGRAM_FILE) == @__FILE__
    include("define_cell.jl")
    test_ensemble_rank2()
    println("Ensemble rank-2 passed")
    test_ensemble_rank3()
    println("Ensemble rank-3 passed")
    test_ensemble_rank4()
    println("Ensemble rank-4 passed")
    test_streaming_api()
    println("Streaming API passed")
    println("All ensemble average tests passed!")
end
