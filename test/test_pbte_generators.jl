using AtomicSymmetries
using Test


function test_find_generators_pbte()
    # Get the supercell structure
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'
    nat = 2

    # Generate the supercell (2x2x2)
    new_pos = zeros(Float64, 3, 16)
    for i in 1:2
        for j in 1:2
            for k in 1:2
                for iat in 1:nat
                    new_pos[:, (i-1)*8 + (j-1)*4 + (k-1)*2 + iat] = positions[:, iat] + [i-1, j-1, k-1]
                end
            end
        end
    end
    new_pos /= 2.0
    cell *= 2.0
    atomic_numbers = [(i - 1)%nat + 1 for i in 1:nat*8]

    # Get the symmetry group
    symmetry_group = get_symmetry_group_from_spglib(new_pos, cell, atomic_numbers)

    fc_generators = AtomicSymmetries.get_matrix_generators(symmetry_group, cell)

    # Try the symmetrization using both generators and symmetry group
    fc = randn(Float64, 16*3, 16*3)
    fc += fc'
    fc_nosym = copy(fc)
    fc_sym = similar(fc)

    symmetrize_fc!(fc, cell, symmetry_group)

    # Now project in the generator basis
    coefficients = zeros(Float64, length(fc_generators))
    AtomicSymmetries.get_coefficients_from_fc!(coefficients, fc_nosym, fc_generators, symmetry_group, cell)
    AtomicSymmetries.get_fc_from_generators!(fc_sym, fc_generators, coefficients, symmetry_group, cell)

    # Compare the two results
    for i in 1:16*3
        for j in 1:16*3
            @test fc_sym[i, j] ≈ fc[i, j] atol = 1e-8 
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_find_generators_pbte()
end
