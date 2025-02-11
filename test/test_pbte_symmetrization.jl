using Test
using LinearAlgebra
using AtomicSymmetries

"""
Test the symmetrization of the PbTe Fm3m group in a 2x2x2 supercell.
This is a task that gives a bug in 0.3.0 version. 
"""
function test_pbte_supercell_symmetrization(; verbose=false)
    atomic_positions_cartesian = [-0.0030913163479886186 0.004423255302789579 0.0025311681149686113; 6.108708886136721 6.112556792389377 6.107500252702037; -6.110258448315757 6.116668060145084 0.0050761535898710824; -0.0027353674835996844 12.21735062841415 6.107479044669904; 0.0051811198475670344 6.104438957319992 6.116427043749321; 6.104006323695241 12.22480694997538 12.216334846421883; -6.108366774622958 12.219181329722256 6.110899178291157; -0.006475079250780129 18.329211612098593 12.220325313413129; -6.108482242351837 0.003518330124867612 6.108705900826501; 7.69635531461288e-5 6.113872750260061 12.217597003523203; -12.209386549212592 6.105574927399339 6.118048977203069; -6.11372570959497 12.220628756341823 12.216440802255619; -6.109033948640472 6.1176075743427605 12.218363995022791; -0.0004241918044713746 12.215662067872604 18.326548153018294; -12.220757822321128 12.220273192306285 12.218793382222819; -6.102919425297092 18.328313408523517 18.3261662099498]'
    cell = [-12.215000876671429 0.0 12.215000876671429; 0.0 12.215000876671429 12.215000876671429; -12.215000876671429 12.215000876671429 0.0]'
    types = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]

    # Convert the coordinates in the primitive cell
    to_primitive_cell_cart!(atomic_positions_cartesian, cell)

    # Convert in crystal coordinates
    crystal_coords = zeros(Float64, size(atomic_positions_cartesian)...)
    get_crystal_coords!(crystal_coords, atomic_positions_cartesian, cell)

    # Extract the symmetries 
    symmetry_group = get_symmetry_group_from_spglib(crystal_coords, cell, types; symprec = 0.1)
    n_symmetries = length(symmetry_group)

    if verbose
        println("Number of symmetries: ", n_symmetries)
        println("Before symmetrization:", atomic_positions_cartesian')
        println()
        println("Before symmetrization crystal:", crystal_coords')
    end

    @test n_symmetries == 384

    # Perform the symmetrization
    symmetrize_positions!(atomic_positions_cartesian, cell, symmetry_group)

    get_crystal_coords!(crystal_coords, atomic_positions_cartesian, cell)

    if verbose
        println("Symmetrized positions: ", crystal_coords')
    end

    # Perform the test. All atoms are in good positions in a 1/4 fractional coordinates
    nat = size(atomic_positions_cartesian, 2)
    for i in 1:nat
        for j in 1:3
            δ = abs(crystal_coords[j, i]*4 - round(crystal_coords[j, i]*4))
            @test δ < 1e-12
        end
    end
end


function test_asr_impose(; verbose=false)
    # Define the PbTe cell 
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'

    # Get the symmetry group and ASR
    sym_group = get_symmetry_group_from_spglib(positions, cell, [1, 2])
    asr! = ASRConstraint!(3)

    if verbose 
        println("Found $(length(sym_group)) symmetries")
    end

    # Initialize a random Hermitian dynamical matrix
    fc = randn(Float64, 6, 6)
    fc += fc'
    asr!(fc)
    symmetrize_fc!(fc, cell, sym_group)

    ω = real.(eigvals(fc))

    if verbose
        println()
        for i in 1:length(sym_group)
            println("Symmetry ", i)
            println("IRT: ", sym_group.irt[i])
        end

        println("Symmetrized force constants:")
        println(fc)
        println("Eigenvalues: ", ω)
    end

    # Check if at least 3 eigenvalues are zero
    zerocount = 0
    for ω_val in ω
        if abs(ω_val) < 1e-12
            zerocount += 1
        end
    end
    @test zerocount == 3

    # Check if the 3 non-zero eigenvalues are the same
    value = 0.0
    for ω_val in ω
        if abs(ω_val) < 1e-12
            continue
        end
        if value == 0
            value = ω_val
        else
            @test ω_val ≈ value atol = 1e-6
        end
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_pbte_supercell_symmetrization(verbose=true)
    test_asr_impose(verbose=true)
end
