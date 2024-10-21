using AtomicSymmetries
using Test

function test_r3m()
    a = 2.87
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    positions = [0.0 0.4
                 0.0 0.4
                 0.0 0.4]
    types = [1, 2]


    # Get the symmetry group
    r3m_group = get_symmetry_group_from_spglib(positions, cell, types)
    println("Number of symmetries: ", get_nsymmetries(r3m_group))



    # Now check if we constrain the symmetry we get a null vector
    vector_input = rand(3)
    @show vector_input
    r3m_group.symmetrize_centroid!(vector_input)
    @show vector_input

    @test maximum(abs.(vector_input .- (sum(vector_input)/3))) < 1e-12
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_r3m()
end
