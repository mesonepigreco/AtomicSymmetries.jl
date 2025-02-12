using AtomicSymmetries
using Test

function test_r3m(; verbose=false)
    a = 2.87
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    positions = [0.0 0.4
                 0.0 0.4
                 0.0 0.4]
    types = [1, 1]


    # Get the symmetry group
    r3m_group = get_symmetry_group_from_spglib(positions, cell, types)
    if verbose
        println("Number of symmetries: ", get_nsymmetries(r3m_group))
    end



    # Now check if we constrain the symmetry we get a null vector
    vector_input = rand(6)
    @show vector_input
    r3m_group.symmetrize_centroid!(vector_input)
    @show vector_input

    for i in 1:2
        @test maximum(abs.(vector_input[3*(i-1)+1:3*i] .- (sum(vector_input[3*(i-1)+1:3*i])/3))) < 1e-12
    end


    # Extract the generators of this group
    generators = AtomicSymmetries.get_vector_generators(r3m_group, cell)
                                                       


    if verbose
        println("This group has $(length(generators)) generators for the vectors")
        for i in 1:length(generators)

            AtomicSymmetries.get_vector_generator!(vector_input, generators[i], r3m_group)
            println("Generator $i: ", reshape(vector_input, 3, 2))
        end
    end


    # Now check the generators of the 2-rank tensors
    tensor_input = rand(6, 6)
    generators = AtomicSymmetries.get_matrix_generators(r3m_group, cell)
    if verbose
        println("This group has $(length(generators)) generators for the matrices")
        for i in 1:length(generators)
            AtomicSymmetries.get_matrix_generator!(tensor_input, generators[i], r3m_group, cell)
            println("Generator $i: ", tensor_input)
        end
    end

    # Try to symmetrize a dynamical matrix and check the eigenvectors
    tensor_input = rand(6, 6)
    tensor_input .+= tensor_input'
    second_input = copy(tensor_input)
    r3m_group.symmetrize_fc!(tensor_input)

    # Now perform the symmetrization projecting into the generators
    coeffs = zeros(length(generators))
    AtomicSymmetries.get_coefficients_from_fc!(coeffs, second_input, generators, r3m_group, cell)
    AtomicSymmetries.get_fc_from_generators!(second_input, generators, coeffs, r3m_group, cell)

    if verbose
        @show tensor_input
        @show second_input
    end

    # Compare the two results
    @test maximum(abs.(tensor_input .- second_input)) < 1e-10
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_r3m(verbose=true)
end
