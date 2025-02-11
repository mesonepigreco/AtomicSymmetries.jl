using AtomicSymmetries
using Test
using DelimitedFiles

function test_fc_sym_primitive_cell(; verbose=false)
    # Load the force constants before symmetrization and after symmetrization
    fc_nosym = readdlm(joinpath(@__DIR__, "fc_cart_small_nosym.txt"))
    fc_sym = readdlm(joinpath(@__DIR__, "fc_cart_small_sym.txt"))


    cell = [1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0]
    positions = [0.0 0.5
                 0.0 0.5
                 0.0 0.5]
    itypes = [1, 1]

    # Get the symmetry group
    sym_group = get_symmetry_group_from_spglib(positions, cell, itypes)
    asr! = ASRConstraint!(3)

    # Apply the symmetries on the force constant matrix
    asr!(fc_nosym)
    symmetrize_fc!(fc_nosym, cell, sym_group)

    if verbose
        println("Symmetrized force constants:")
        println(fc_nosym)
    end
    
    # Check if the force constant matrix is the same
    # as the correct one
    
    for i in 1:size(fc_sym, 1)
        for j in 1:size(fc_sym, 2)
            @test fc_nosym[i, j] ≈ fc_sym[i, j] atol = 1e-6
        end
    end
end

function test_fc_sym_pbte_uc(; verbose=false)
    # Define the PbTe cell 
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'
    atoms_types = [1, 2]
    nat = length(atoms_types)

    # Load the force constants before symmetrization and after symmetrization
    fc_nosym = readdlm(joinpath(@__DIR__, "fc_cart_small_nosym.txt"))
    fc_sym = readdlm(joinpath(@__DIR__, "fc_cart_pbte_sym.txt"))

    # Convert the matrix to crystal coordinates
    fc_nosym_crystal = zeros(Float64, size(fc_nosym)...)
    test_cartesian = zeros(Float64, size(fc_nosym)...)
    AtomicSymmetries.cart_cryst_matrix_conversion!(fc_nosym_crystal, fc_nosym, cell; cart_to_cryst = true)
    AtomicSymmetries.cart_cryst_matrix_conversion!(test_cartesian, fc_nosym_crystal, cell; cart_to_cryst = false)

    # Test that the conversion is correct
    for i in 1:size(fc_nosym, 1)
        for j in 1:size(fc_nosym, 2)
            @test fc_nosym[i, j] ≈ test_cartesian[i, j] atol = 1e-6
        end
    end


    # Perform the symmetrization
    symmetry_group = get_symmetry_group_from_spglib(positions, cell, atoms_types)
    asr! = ASRConstraint!(3)
    asr!(fc_nosym)
    symmetrize_fc!(fc_nosym, cell, symmetry_group)

    if verbose
        println()
        println("TEST PBTE FC SYMMETRIZATION")
        println("===========================")
        println("Number of symmetries: ", length(symmetry_group))
        println("FC in crystal coordinates:")
        for i in 1:nat
            for j in 1:nat
                println("atom i = $i, atom j = $j")
                println(fc_nosym_crystal[3*(i-1)+1:3*i, 3*(j-1)+1:3*j])
            end
        end
        println("Symmetrized force constants:")
        println(fc_nosym)
    end

    # Check if the force constant matrix is the same
    # as the correct one
    for i in 1:size(fc_sym, 1)
        for j in 1:size(fc_sym, 2)
            @test fc_nosym[i, j] ≈ fc_sym[i, j] atol = 1e-6
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_fc_sym_primitive_cell(verbose=true)
    test_fc_sym_pbte_uc(verbose=true)
end
