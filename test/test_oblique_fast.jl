using AtomicSymmetries
using LinearAlgebra
using Random
using Test

include("define_cell.jl")

@testset "Oblique Stress Robustness (Fast vs Slow)" begin
    # === Highly oblique triclinic cell with 2 atoms ===
    a, b, c = 3.0, 4.5, 5.2
    α, β, γ = 70.0, 65.0, 80.0
    αr, βr, γr = deg2rad.([α, β, γ])

    cell = zeros(3, 3)
    cell[:, 1] = [a, 0, 0]
    cell[:, 2] = [b * cos(γr), b * sin(γr), 0]
    cx = c * cos(βr)
    cy = c * (cos(αr) - cos(βr) * cos(γr)) / sin(γr)
    cz = sqrt(c^2 - cx^2 - cy^2)
    cell[:, 3] = [cx, cy, cz]

    positions = [0.0 0.3; 0.0 0.4; 0.0 0.2]
    types = [1, 2]

    sym = get_symmetry_group_from_spglib(positions, cell, types)
    println("Triclinic: N_sym = ", length(sym.symmetries))

    gens_slow = get_tensor_generators(sym, cell; rank=2)
    gens_fast = get_tensor_generators_fast(sym, cell; rank=2)

    println("Triclinic Rank 2: slow=$(length(gens_slow)), fast=$(length(gens_fast))")
    @test length(gens_slow) == length(gens_fast)
    @test length(gens_fast) == 21  # 6(1,1) + 6(2,2) + 9(1,2)

    # === Rotated Triclinic ===
    Random.seed!(1)
    θ, φ, ψ = 0.5, 0.4, 0.3
    Rz1 = [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]
    Ry = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
    Rz2 = [cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1]
    R_rand = Rz1 * Ry * Rz2
    rotated_triclinic = R_rand * cell

    sym_rot_tri = get_symmetry_group_from_spglib(positions, rotated_triclinic, types)
    gens_fast_rot_tri = get_tensor_generators_fast(sym_rot_tri, rotated_triclinic; rank=2)
    println("Rotated Triclinic Rank 2: fast=$(length(gens_fast_rot_tri))")
    @test length(gens_fast_rot_tri) == length(gens_fast)

    # === Randomly rotated FCC Gold ===
    pos_au, cell_au, types_au = get_gold_unit_cell()
    Random.seed!(42)
    θ = 1.23
    φ = 0.87
    ψ = 2.45
    Rz1 = [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]
    Ry = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
    Rz2 = [cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1]
    R_rand = Rz1 * Ry * Rz2
    rotated_cell = R_rand * cell_au

    sym_rot = get_symmetry_group_from_spglib(pos_au, rotated_cell, types_au)
    println("\nRotated FCC: N_sym = ", length(sym_rot.symmetries))

    gens_slow_rot = get_tensor_generators(sym_rot, rotated_cell; rank=2)
    gens_fast_rot = get_tensor_generators_fast(sym_rot, rotated_cell; rank=2)

    println("Rotated FCC Rank 2: slow=$(length(gens_slow_rot)), fast=$(length(gens_fast_rot))")
    @test length(gens_slow_rot) == length(gens_fast_rot)

    # Orthogonality check for fast generators in rotated cell
    for i in 1:length(gens_fast_rot), j in i+1:length(gens_fast_rot)
        dot_val = AtomicSymmetries._generator_dot(gens_fast_rot[i], gens_fast_rot[j])
        @test abs(dot_val) < 1e-10
    end

    # === Rotated Gold Supercell 3x3x3 (Rank 2 and 3) ===
    pos_au, cell_au, types_au = get_gold_unit_cell()
    pos_sup, cell_sup, types_sup = get_supercell(pos_au, cell_au, types_au, [3, 3, 3])

    sym_sup = get_symmetry_group_from_spglib(pos_sup, cell_sup, types_sup)
    gens2_fast = get_tensor_generators_fast(sym_sup, cell_sup; rank=2)
    gens3_fast = get_tensor_generators_fast(sym_sup, cell_sup; rank=3)

    println("\nGold Supercell 3x3x3: Rank 2 = $(length(gens2_fast)), Rank 3 = $(length(gens3_fast))")

    # Rotate by a small random angle
    Random.seed!(777)
    θ, φ, ψ = 0.1, 0.15, 0.05
    Rz1 = [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]
    Ry = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
    Rz2 = [cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1]
    R_rand = Rz1 * Ry * Rz2

    rotated_cell_sup = R_rand * cell_sup

    sym_rot_sup = get_symmetry_group_from_spglib(pos_sup, rotated_cell_sup, types_sup)
    gens2_rot_fast = get_tensor_generators_fast(sym_rot_sup, rotated_cell_sup; rank=2)
    gens3_rot_fast = get_tensor_generators_fast(sym_rot_sup, rotated_cell_sup; rank=3)

    println("Rotated Gold Supercell Rank 2: $(length(gens2_rot_fast)) (orig: $(length(gens2_fast)))")
    println("Rotated Gold Supercell Rank 3: $(length(gens3_rot_fast)) (orig: $(length(gens3_fast)))")

    @test length(gens2_rot_fast) == length(gens2_fast)
    @test length(gens3_rot_fast) == length(gens3_fast)
end
