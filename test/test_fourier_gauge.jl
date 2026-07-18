using AtomicSymmetries
using LinearAlgebra
using Random
using Test

# Test suite for the atomic-position phase gauge of the Fourier transform:
#
#   v_a(q) ~ e^{-2πi q·(R + τ_a)},   M_ab(q) ~ e^{+2πi q·(R_a + τ_a - R_b - τ_b)}
#
# The tests use a CsCl-like structure (simple cubic, atom 1 at the origin,
# atom 2 at (0.5, 0.5, 0.5) with a different type) on a 2x2x2 supercell.
# This structure has the full Oh point group (including inversion), and since
# τ₂ = (½,½,½) is nontrivial, the q-grid foldings S^{-T}q = q' + G generate
# nontrivial phase factors e^{2πi G·τ} = -1: the sign conventions of all the
# phases are therefore actually exercised.
function get_cscl_supercell()
    a = 3.1
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    positions = [0.0 0.5
                 0.0 0.5
                 0.0 0.5]
    types = [1, 2]

    supercell = [2, 2, 2]
    nat = size(positions, 2)
    ndims = 3
    n_sc = prod(supercell)
    nat_sc = nat * n_sc

    R_lat = zeros(Float64, ndims, nat_sc)
    q_vec = zeros(Float64, ndims, n_sc)
    super_cell = zeros(Float64, ndims, ndims)
    super_positions = zeros(Float64, ndims, nat_sc)
    super_types = ones(Int, nat_sc)
    super_itau = zeros(Int, nat_sc)

    for i in 1:ndims
        @views super_cell[:, i] .= cell[:, i] * supercell[i]
    end

    counter = 1
    for i in 1:supercell[1]
        for j in 1:supercell[2]
            for k in 1:supercell[3]
                latvec = [i-1, j-1, k-1] ./ supercell
                @views q_vec[:, counter] .= latvec

                for h in 1:nat
                    R_lat[:, nat * (counter - 1) + h] = latvec .* supercell
                    super_positions[:, nat * (counter - 1) + h] = positions[:, h] ./ supercell + latvec
                    super_itau[nat * (counter - 1) + h] = h
                end
                counter += 1
            end
        end
    end

    for iat in 1:nat_sc
        super_types[iat] = types[super_itau[iat]]
    end

    return (cell = cell, positions = positions, types = types,
            supercell = supercell, super_cell = super_cell,
            super_positions = super_positions, super_types = super_types,
            super_itau = super_itau, R_lat = R_lat, q_vec = q_vec)
end


# Build a random translation-invariant symmetric force-constant matrix
function get_random_invariant_fc(super_positions, super_itau, supercell, super_cell, super_types)
    Random.seed!(5678)
    ndims, nat_sc = size(super_positions)

    sc_group = get_symmetry_group_from_spglib(super_positions, super_cell, super_types)
    translations = get_translations(sc_group)

    fc = randn(Float64, ndims * nat_sc, ndims * nat_sc)
    fc .+= fc'
    apply_translations!(fc, translations)
    return fc, translations
end


# Check that the new gauge is related to the lattice gauge (phases with R only)
# by the rephasing M̃_ab(q) = e^{2πi q·(τ_a - τ_b)} M_ab(q)
function test_gauge_phase_relation(; verbose=false)
    s = get_cscl_supercell()
    ndims = 3
    nat = size(s.positions, 2)
    nat_sc = size(s.super_positions, 2)
    nq = size(s.q_vec, 2)

    fc, translations = get_random_invariant_fc(s.super_positions, s.super_itau, s.supercell,
                                               s.super_cell, s.super_types)

    phi_q = zeros(ComplexF64, ndims*nat, ndims*nat, nq)
    matrix_r2q!(phi_q, fc, s.q_vec, s.super_itau, s.R_lat, s.positions)

    # Reference: lattice-gauge Fourier transform computed by hand
    phi_q_lattice = zeros(ComplexF64, ndims*nat, ndims*nat, nq)
    for iq in 1:nq
        for k in 1:nat
            for h in 1:nat_sc
                h_uc = s.super_itau[h]
                ΔR = s.R_lat[:, k] .- s.R_lat[:, h]
                phase = exp(-2im * π * (s.q_vec[:, iq]' * ΔR))
                @views phi_q_lattice[ndims*(h_uc-1)+1:ndims*h_uc, ndims*(k-1)+1:ndims*k, iq] .+=
                    phase .* fc[ndims*(h-1)+1:ndims*h, ndims*(k-1)+1:ndims*k]
            end
        end
    end

    for iq in 1:nq
        for b in 1:nat
            for a in 1:nat
                gauge_phase = exp(2im * π * (s.q_vec[:, iq]' * (s.positions[:, a] .- s.positions[:, b])))
                @views expected = gauge_phase .* phi_q_lattice[ndims*(a-1)+1:ndims*a, ndims*(b-1)+1:ndims*b, iq]
                @views @test isapprox(phi_q[ndims*(a-1)+1:ndims*a, ndims*(b-1)+1:ndims*b, iq],
                                      expected; atol=1e-10)
            end
        end
    end
end


# Round trip r -> q -> r -> q for the matrix and vector transforms
function test_gauge_roundtrip(; verbose=false)
    s = get_cscl_supercell()
    ndims = 3
    nat = size(s.positions, 2)
    nat_sc = size(s.super_positions, 2)
    nq = size(s.q_vec, 2)

    fc, translations = get_random_invariant_fc(s.super_positions, s.super_itau, s.supercell,
                                               s.super_cell, s.super_types)

    phi_q = zeros(ComplexF64, ndims*nat, ndims*nat, nq)
    phi_back = zeros(Float64, ndims*nat_sc, ndims*nat_sc)
    phi_q_bis = similar(phi_q)

    matrix_r2q!(phi_q, fc, s.q_vec, s.super_itau, s.R_lat, s.positions)
    matrix_q2r!(phi_back, phi_q, s.q_vec, s.super_itau, s.R_lat, s.positions; translations=translations)
    @test isapprox(phi_back, fc; atol=1e-10)

    matrix_r2q!(phi_q_bis, phi_back, s.q_vec, s.super_itau, s.R_lat, s.positions)
    @test isapprox(phi_q_bis, phi_q; atol=1e-10)

    # Vector roundtrip
    Random.seed!(91011)
    u = randn(Float64, ndims * nat_sc)
    u_q = zeros(ComplexF64, ndims*nat, nq)
    u_back = zeros(Float64, ndims * nat_sc)
    vector_r2q!(u_q, u, s.q_vec, s.super_itau, s.R_lat, s.positions)
    vector_q2r!(u_back, u_q, s.q_vec, s.super_itau, s.R_lat, s.positions)
    @test isapprox(u_back, u; atol=1e-10)

    # absolute_positions: transforming (equilibrium + displacement) with
    # absolute_positions=true must match transforming the displacement alone.
    # Here everything is in (primitive) crystal coordinates: the equilibrium
    # positions are R_lat + tau.
    u_abs = copy(u)
    for k in 1:nat_sc
        @views u_abs[ndims*(k-1)+1:ndims*k] .+= s.R_lat[:, k] .+ s.positions[:, s.super_itau[k]]
    end
    u_q_abs = zeros(ComplexF64, ndims*nat, nq)
    vector_r2q!(u_q_abs, u_abs, s.q_vec, s.super_itau, s.R_lat, s.positions; absolute_positions=true)
    @test isapprox(u_q_abs, u_q; atol=1e-10)

    # ... and the backward transform restores the absolute positions
    u_back_abs = zeros(Float64, ndims * nat_sc)
    vector_q2r!(u_back_abs, u_q_abs, s.q_vec, s.super_itau, s.R_lat, s.positions; absolute_positions=true)
    @test isapprox(u_back_abs, u_abs; atol=1e-10)
end


# Apply each symmetry both in real space and in q space and compare.
# Since τ₂ = (0.5, 0.5, 0.5) and the Oh group folds the q grid, the
# G-folding phases are nontrivial here.
function test_gauge_symmetry_application(; verbose=false)
    s = get_cscl_supercell()
    ndims = 3
    nat = size(s.positions, 2)
    nat_sc = size(s.super_positions, 2)
    nq = size(s.q_vec, 2)

    fc, translations = get_random_invariant_fc(s.super_positions, s.super_itau, s.supercell,
                                               s.super_cell, s.super_types)

    uc_group = get_symmetry_group_from_spglib(s.positions, s.cell, s.types)
    n_sym = length(uc_group)
    if verbose
        println("Number of symmetries: $n_sym")
    end
    @test n_sym == 48

    phi_q = zeros(ComplexF64, ndims*nat, ndims*nat, nq)
    matrix_r2q!(phi_q, fc, s.q_vec, s.super_itau, s.R_lat, s.positions)

    irt_q = zeros(Int, nq)
    irt_sc = zeros(Int, nat_sc)
    trans_vect = zeros(Float64, ndims, nat_sc)

    phi_q_sym = similar(phi_q)
    fc_sym = zeros(Float64, ndims*nat_sc, ndims*nat_sc)
    phi_q_ref = similar(phi_q)

    n_nontrivial_foldings = 0

    for i_sym in 1:n_sym
        sym_mat = uc_group.symmetries[i_sym]
        AtomicSymmetries.get_irt_q!(irt_q, s.q_vec, sym_mat)

        # Count the q points folded back with G != 0 (to make sure the phases matter)
        sym_rec = inv(sym_mat)'
        for iq in 1:nq
            G = sym_rec * s.q_vec[:, iq] .- s.q_vec[:, irt_q[iq]]
            if maximum(abs.(G)) > 1e-6
                n_nontrivial_foldings += 1
            end
        end

        # q-space application
        phi_q_sym .= 0
        AtomicSymmetries.apply_symmetry_matrixq!(phi_q_sym, phi_q, sym_mat,
                                                 uc_group.irt[i_sym], irt_q,
                                                 s.positions, s.q_vec)

        # Real-space application
        AtomicSymmetries.get_irt!(irt_sc, trans_vect, s.super_positions, sym_mat,
                                  uc_group.translations[i_sym] ./ s.supercell)
        fc_sym .= 0
        AtomicSymmetries.apply_sym_fc!(fc_sym, fc, sym_mat, ndims, irt_sc)
        apply_translations!(fc_sym, translations)

        phi_q_ref .= 0
        matrix_r2q!(phi_q_ref, fc_sym, s.q_vec, s.super_itau, s.R_lat, s.positions)

        @test isapprox(phi_q_sym, phi_q_ref; atol=1e-8)
    end

    # The Oh group on the 2x2x2 grid must fold many q points outside the grid
    @test n_nontrivial_foldings > 0
    if verbose
        println("Number of foldings with G != 0: $n_nontrivial_foldings")
    end

    # Full symmetrization: q space vs real space
    q_symmetries = SymmetriesQSpace(uc_group, s.q_vec, s.positions)
    phi_q_sym .= 0
    symmetrize_matrix_q!(phi_q_sym, phi_q, q_symmetries)

    sc_group = get_symmetry_group_from_spglib(s.super_positions, s.super_cell, s.super_types)
    fc_sym .= fc
    sc_group.symmetrize_fc!(fc_sym)
    phi_q_ref .= 0
    matrix_r2q!(phi_q_ref, fc_sym, s.q_vec, s.super_itau, s.R_lat, s.positions)

    @test isapprox(phi_q_sym, phi_q_ref; atol=1e-8)

    # Negative control: ignoring the atomic positions (i.e. using the
    # lattice-gauge symmetrization on a new-gauge matrix) must give a
    # different (wrong) result: the folding phases really matter here.
    fake_positions = zeros(Float64, ndims, nat)
    fake_q_symmetries = SymmetriesQSpace(uc_group, s.q_vec, fake_positions)
    phi_q_wrong = similar(phi_q)
    phi_q_wrong .= 0
    symmetrize_matrix_q!(phi_q_wrong, phi_q, fake_q_symmetries)
    @test maximum(abs.(phi_q_wrong .- phi_q_ref)) > 1e-6
end


# Apply each symmetry to a displacement field both in real space and in
# q space (with the folding + fractional translation phases) and compare.
function test_gauge_vector_symmetry(; verbose=false)
    s = get_cscl_supercell()
    ndims = 3
    nat = size(s.positions, 2)
    nat_sc = size(s.super_positions, 2)
    nq = size(s.q_vec, 2)

    uc_group = get_symmetry_group_from_spglib(s.positions, s.cell, s.types)
    n_sym = length(uc_group)

    Random.seed!(1213)
    u = randn(Float64, ndims * nat_sc)
    u_q = zeros(ComplexF64, ndims*nat, nq)
    vector_r2q!(u_q, u, s.q_vec, s.super_itau, s.R_lat, s.positions)

    irt_q = zeros(Int, nq)
    irt_sc = zeros(Int, nat_sc)
    trans_vect = zeros(Float64, ndims, nat_sc)
    u_q_sym = similar(u_q)
    u_sym = zeros(Float64, ndims * nat_sc)
    u_q_ref = similar(u_q)

    for i_sym in 1:n_sym
        sym_mat = uc_group.symmetries[i_sym]
        AtomicSymmetries.get_irt_q!(irt_q, s.q_vec, sym_mat)

        # q-space application (with all the gauge phases)
        u_q_sym .= 0
        AtomicSymmetries.apply_symmetry_vectorq!(u_q_sym, u_q, sym_mat,
                                                 uc_group.irt[i_sym], irt_q;
                                                 positions=s.positions,
                                                 q_points=s.q_vec,
                                                 translation=uc_group.translations[i_sym])

        # Real-space application (displacements: no translations applied to the values)
        AtomicSymmetries.get_irt!(irt_sc, trans_vect, s.super_positions, sym_mat,
                                  uc_group.translations[i_sym] ./ s.supercell)
        u_sym .= 0
        AtomicSymmetries.apply_sym_centroid!(u_sym, u, sym_mat, ndims, irt_sc)

        u_q_ref .= 0
        vector_r2q!(u_q_ref, u_sym, s.q_vec, s.super_itau, s.R_lat, s.positions)

        @test isapprox(u_q_sym, u_q_ref; atol=1e-8)
    end
end


# The Fourier transform of a real, symmetric, translation-invariant
# force-constant matrix must already satisfy hermitianity + time reversal:
# impose_hermitianity_q! must leave it unchanged (and be idempotent).
function test_gauge_hermitianity(; verbose=false)
    s = get_cscl_supercell()
    ndims = 3
    nat = size(s.positions, 2)
    nq = size(s.q_vec, 2)

    fc, translations = get_random_invariant_fc(s.super_positions, s.super_itau, s.supercell,
                                               s.super_cell, s.super_types)

    uc_group = get_symmetry_group_from_spglib(s.positions, s.cell, s.types)
    q_symmetries = SymmetriesQSpace(uc_group, s.q_vec, s.positions)

    phi_q = zeros(ComplexF64, ndims*nat, ndims*nat, nq)
    matrix_r2q!(phi_q, fc, s.q_vec, s.super_itau, s.R_lat, s.positions)

    phi_q_fixed = copy(phi_q)
    impose_hermitianity_q!(phi_q_fixed, q_symmetries)
    @test isapprox(phi_q_fixed, phi_q; atol=1e-10)

    # Idempotency on a generic (random) matrix
    Random.seed!(1415)
    random_q = randn(ComplexF64, ndims*nat, ndims*nat, nq)
    once = copy(random_q)
    impose_hermitianity_q!(once, q_symmetries)
    twice = copy(once)
    impose_hermitianity_q!(twice, q_symmetries)
    @test isapprox(twice, once; atol=1e-10)

    # Negative control: without the positions, the folding phases are missed
    # and the FT of a valid real-space matrix is (wrongly) modified
    phi_q_wrong = copy(phi_q)
    impose_hermitianity_q!(phi_q_wrong, q_symmetries.minus_q_index)
    @test maximum(abs.(phi_q_wrong .- phi_q)) > 1e-6
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_gauge_phase_relation(; verbose=true)
    test_gauge_roundtrip(; verbose=true)
    test_gauge_symmetry_application(; verbose=true)
    test_gauge_vector_symmetry(; verbose=true)
    test_gauge_hermitianity(; verbose=true)
end
