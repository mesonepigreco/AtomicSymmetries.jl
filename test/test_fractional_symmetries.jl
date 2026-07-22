using AtomicSymmetries
using Test
using DelimitedFiles
using Unitful, UnitfulAtomic
using LinearAlgebra
using PhysicalConstants


function test_fractional_symmetries_qspace(; verbose=false)
    # Load the data files on LaAlO3
    Φ_sc = readdlm(joinpath(@__DIR__, "data", "LaAlO3_phisc.txt"))
    unit_cell = readdlm(joinpath(@__DIR__, "data", "LaAlO3_unitcell.txt"))
    coords_sc_ = readdlm(joinpath(@__DIR__, "data", "LaAlO3_coords_sc.txt"))
    itau_ = readdlm(joinpath(@__DIR__, "data", "LaAlO3_itau.txt"))
    qpoints = readdlm(joinpath(@__DIR__, "data", "LaAlO3_qpoints.txt"))
    Rlat = readdlm(joinpath(@__DIR__, "data", "LaAlO3_Rlat.txt"))
    types_sc_ = readdlm(joinpath(@__DIR__, "data", "LaAlO3_types_sc.txt"))
    types_sc = [Int.(types_sc_[:])...]
    itau = [Int.(itau_[:])...]

    supercell_dim = 4

    coords_sc = reshape(coords_sc_, 3, :)
    primitive_cell = copy(unit_cell)
    unit_cell .*= supercell_dim
    cryst_coords = similar(coords_sc)
    reciprocal_vectors = similar(unit_cell)
    get_reciprocal_lattice!(reciprocal_vectors, unit_cell)
    q_points = similar(qpoints)
    rec_vect = similar(reciprocal_vectors)
    get_reciprocal_lattice!(rec_vect, unit_cell)


    cryst_cart_conv!(cryst_coords, coords_sc, unit_cell, reciprocal_vectors, false; q_space=false)
    cryst_cart_conv!(q_points, qpoints, primitive_cell, rec_vect, false; q_space=true) 

    println(cryst_coords')
    @show q_points

    
    # Get the symmetries in the supercell
    symmetries_supercell = get_symmetry_group_from_spglib(cryst_coords, unit_cell, types_sc)

    # Get the dynamical matrix in q space
    n_q = size(qpoints, 2)
    n_dims = size(qpoints, 1)
    nat_sc = size(coords_sc, 2)
    nat = nat_sc ÷ n_q
    
    # Get the unit cell coordinates
    coords_uc = coords_sc[:, 1:nat]
    cryst_coords = similar(coords_uc)
    types_uc = types_sc[1:nat]

    get_reciprocal_lattice!(reciprocal_vectors, primitive_cell)
    cryst_coords = similar(coords_uc)
    cryst_cart_conv!(cryst_coords, coords_uc, primitive_cell, reciprocal_vectors, false; q_space=false)

    

    symmetries_uc = get_symmetry_group_from_spglib(cryst_coords, primitive_cell, types_uc)

    Φ_q = zeros(ComplexF64, n_dims * nat, n_dims * nat, n_q)
    matrix_r2q!(Φ_q, Φ_sc, qpoints, itau, Rlat)
    symmetries_qspace = SymmetriesQSpace(symmetries_uc, q_points)


    # Perform the symmetrization in real space and in q_space, then compare
    symmetrize_fc!(Φ_sc, unit_cell, symmetries_supercell)
    symmetrize_matrix_cartesian_q!(Φ_q, primitive_cell, symmetries_qspace)

    if verbose
        println("symmetries in supercell: $(length(symmetries_supercell))")
        println("primitive cell = ", primitive_cell)
        println("IRTQ:")
        println("q_points = $qpoints")
        n_sym = length(symmetries_qspace.symmetries)
        for i in 1:n_sym
            println("sym = $i ; irt_q = $(symmetries_qspace.irt_q[i])")
        end
    end

    # Convert the supercell force constants to q space for comparison
    Φ_q_real_space_sym = zeros(ComplexF64, n_dims * nat, n_dims * nat, n_q)
    matrix_r2q!(Φ_q_real_space_sym, Φ_sc, qpoints, itau, Rlat)

    hplack = PhysicalConstants.CODATA2018.h

    for i in 1:n_q
        if verbose

            println("Comparing q-point $i / $n_q")
            println("Real space symmetrized FCs =")
            ω_sc = eigvals(Φ_q_real_space_sym[:, :, i])
            println(sqrt.(abs.(ω_sc)) * ustrip(uconvert(u"c/cm", 1.0u"hartree/h")))

            println("Q-space symmetrized FCs =")
            ω_q = eigvals(Φ_q[:, :, i])
            println(sqrt.(abs.(ω_q)) * ustrip(uconvert(u"c/cm", 1.0u"hartree/h")))

        end
        @test isapprox(Φ_q[:, :, i], Φ_q_real_space_sym[:, :, i]; atol=1e-10, rtol=1e-6)

        if verbose
            delta_minus_q = maximum(abs.(Φ_q[:, :, i] - conj.(Φ_q[:, :, symmetries_qspace.minus_q_index[i]])))
            println("iq = $i; -q mapped to iq = ", symmetries_qspace.minus_q_index[i])
            println("q = $(symmetries_qspace.q_points[:, i])")
            println("minus q = $(symmetries_qspace.q_points[:, symmetries_qspace.minus_q_index[i]])")
            println("Δ = $delta_minus_q")
            println()
        end

        @test isapprox(Φ_q[:, :, i], conj.(Φ_q[:, :, symmetries_qspace.minus_q_index[i]]); rtol=1e-10, atol=1e-12)
    end

end


if abspath(PROGRAM_FILE) == @__FILE__
    test_fractional_symmetries_qspace(; verbose=true)
end


