using AtomicSymmetries
using LinearAlgebra
using Random
using Test

include("define_cell.jl")

@testset "Little Group Supercell" begin
    # Setup a Perovskite supercell (difficult because of high symmetry and multiple atom types)
    supercell_size = [2, 2, 2]
    positions_uc, cell_uc, types_uc = get_pm3m_perovskite()
    positions_sc, cell_sc, types_sc = get_supercell(positions_uc, cell_uc, types_uc, supercell_size)

    # Rotate the supercell randomly to make it "difficult" (robustness check)
    Random.seed!(42)
    θ = rand() * 2π
    φ = rand() * π
    ψ = rand() * 2π
    Rz1 = [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]
    Ry = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
    Rz2 = [cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1]
    Rot = Rz1 * Ry * Rz2

    rotated_cell = Rot * cell_sc
    # Coordinates in supercell are already crystal, so rotation only affects the cell matrix
    # but we should ensure symmetries are computed for the rotated cell.

    sym_group = get_symmetry_group_from_spglib(positions_sc, rotated_cell, types_sc)
    n_sym = length(sym_group)
    println("Number of symmetries in 2x2x2 Perovskite supercell: $n_sym")

    # Test case 1: Single atom (site symmetry)
    # The identity should always be in the little group
    at1 = [1]
    lg1 = get_little_group(at1, sym_group)
    @test 1 in lg1

    # Test case 2: Stabilizer of a pair
    # Take two atoms and check if the little group is a subgroup
    at2 = [1, 2]
    lg2 = get_little_group(at2, sym_group)
    @test 1 in lg2

    # Test case 3: Rank-4 tuple (as requested for high-rank tensors)
    at4 = [1, 2, 3, 4]
    lg4 = get_little_group(at4, sym_group)
    @test 1 in lg4

    # Verify group property: if s1 and s2 are in LG, then s1*s2 should be in LG
    # (Note: sym_group.symmetries stores rotation matrices. Composing them
    # and finding the index in the group is a good check)
    if length(lg4) > 1
        s1_idx = lg4[1]
        s2_idx = lg4[end]
        R1 = sym_group.symmetries[s1_idx]
        R2 = sym_group.symmetries[s2_idx]
        R_prod = R1 * R2

        # Find which symmetry index corresponds to R_prod
        prod_idx = findfirst(S -> isapprox(S, R_prod, atol=1e-8), sym_group.symmetries)
        @test prod_idx !== nothing
        @test prod_idx in lg4
    end

    # Test permutational vs non-permutational
    at_perm = [1, 2]
    # If the symmetry swaps 1 and 2, it's in the permutational LG but not in the non-permutational one
    lg_p = get_little_group(at_perm, sym_group, is_permutational=true)
    lg_np = get_little_group(at_perm, sym_group, is_permutational=false)

    @test issubset(lg_np, lg_p)

    println("Little Group for tuple $at4 has size $(length(lg4))")
end
