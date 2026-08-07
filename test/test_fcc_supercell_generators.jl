using AtomicSymmetries
using Test
using LinearAlgebra
using Random

"""
Test that rank-2 generators on a randomly rotated FCC 4×4×4 supercell produce
the same symmetrized matrix as direct symmetrization via symmetrize_fc!.

A random rotation is applied to the FCC primitive cell so that no lattice vector
is aligned with the Cartesian axes, exercising the coordinate-conversion paths
in the generator machinery.
"""
function test_fcc_444_rotated_generators()
    Random.seed!(42)

    # ── FCC primitive cell (columns = lattice vectors) ──
    a = 4.05  # Å (aluminium-like)
    fcc_prim = (a / 2) * [0.0 1.0 1.0;
                           1.0 0.0 1.0;
                           1.0 1.0 0.0]'

    # Random proper rotation via QR decomposition
    Q, _ = qr(randn(3, 3))
    R = Matrix(Q)
    if det(R) < 0
        R[:, 1] .*= -1
    end
    prim_cell = R * fcc_prim

    # ── 4×4×4 supercell ──
    sc = [4, 4, 4]
    nat_sc = prod(sc)  # 64
    sc_cell = hcat([sc[i] * prim_cell[:, i] for i in 1:3]...)

    # Fractional positions in the supercell
    positions = zeros(3, nat_sc)
    idx = 0
    for ix in 0:sc[1]-1, iy in 0:sc[2]-1, iz in 0:sc[3]-1
        idx += 1
        positions[:, idx] = [ix / sc[1], iy / sc[2], iz / sc[3]]
    end
    types = ones(Int, nat_sc)

    # ── Symmetry group ──
    symmetry_group = get_symmetry_group_from_spglib(positions, sc_cell, types)
    n_sym = get_nsymmetries(symmetry_group)
    @test n_sym > 1

    dim = 3
    n_modes = dim * nat_sc  # 192

    # ── Rank-2 generators ──
    generators = get_tensor_generators(symmetry_group, sc_cell; rank=2)
    @test length(generators) > 0

    # ── Random symmetric force-constant matrix ──
    fc_random = randn(n_modes, n_modes)
    fc_random = (fc_random + fc_random') / 2

    # Method 1: project onto generator basis and reconstruct
    coeffs = zeros(length(generators))
    get_coefficients_from_tensor!(coeffs, fc_random, generators, sc_cell)
    fc_from_generators = zeros(n_modes, n_modes)
    reconstruct_tensor!(fc_from_generators, generators, coeffs, sc_cell)

    # Method 2: direct symmetrization
    fc_direct_sym = copy(fc_random)
    symmetrize_fc!(fc_direct_sym, sc_cell, symmetry_group)

    # Both methods must agree
    @test fc_from_generators ≈ fc_direct_sym atol = 1e-8

    # Verify that the result is actually symmetric (fc symmetry)
    @test fc_from_generators ≈ fc_from_generators' atol = 1e-10

    # Verify idempotency: symmetrizing again should not change the result
    fc_idem = copy(fc_from_generators)
    symmetrize_fc!(fc_idem, sc_cell, symmetry_group)
    @test fc_from_generators ≈ fc_idem atol = 1e-10
end

if abspath(PROGRAM_FILE) == @__FILE__
    include("define_cell.jl")
    test_fcc_444_rotated_generators()
    println("FCC 4×4×4 rotated generators test passed!")
end
