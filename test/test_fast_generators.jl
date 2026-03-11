using AtomicSymmetries
using LinearAlgebra
using Random
using Test

include("define_cell.jl")

@testset "Fast Rank-4 Generators" begin
    # 1. Verification against slow version for Rank 2 (Elastic Constants like)
    # Use BCC Iron (simple)
    pos_bcc = [0.0 0.5; 0.0 0.5; 0.0 0.5]
    cell_bcc = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] .* 2.86
    types_bcc = [1, 1]
    sym_bcc = get_symmetry_group_from_spglib(pos_bcc, cell_bcc, types_bcc)

    println("--- Rank 2 Verification ---")
    gens_slow = get_tensor_generators(sym_bcc, cell_bcc; rank=2)
    gens_fast = get_tensor_generators_fast(sym_bcc, cell_bcc; rank=2)

    println("Rank 2: slow=$(length(gens_slow)), fast=$(length(gens_fast))")
    @test length(gens_slow) == length(gens_fast)

    # 2. Performance and Robustness Rank 4 (High Rank)
    # Use any supercell, e.g. 2x2x2 Gold
    println("\n--- Rank 4 Performance (Rotated Cell) ---")
    pos_au, cell_au, types_au = get_gold_unit_cell()

    # Rotate cell
    Random.seed!(123)
    R_rand = qr(randn(3, 3)).Q
    rot_cell = R_rand * cell_au

    sym_au = get_symmetry_group_from_spglib(pos_au, rot_cell, types_au)

    gens_slow4 = get_tensor_generators(sym_au, rot_cell; rank=4)
    time_fast = @elapsed gens_fast4 = get_tensor_generators_fast(sym_au, rot_cell; rank=4)
    println("Rank 4 (Gold FCC): slow=$(length(gens_slow4)), fast=$(length(gens_fast4)) in $(round(time_fast, digits=2))s")

    @test length(gens_slow4) == length(gens_fast4)

    # 3. Orthogonality test
    if !isempty(gens_fast4)
        # Using the internal dot product helper for generators
        overlap = AtomicSymmetries._generator_dot(gens_fast4[1], gens_fast4[2])
        println("Orthogonality check (1,2): $overlap")
        @test abs(overlap) < 1e-10

        # Norm check
        n1 = AtomicSymmetries._generator_dot(gens_fast4[1], gens_fast4[1])
        println("Norm check (1): $n1")
        @test abs(n1 - 1.0) < 1e-10
    end
end
