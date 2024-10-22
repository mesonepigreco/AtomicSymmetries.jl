using AtomicSymmetries
using Test

function test_bcc()
    a = 2.87
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    positions = [0.0 0.5
                 0.0 0.5
                 0.0 0.5]
    types = [1, 1]


    # Get the symmetry group
    bcc_group = get_symmetry_group_from_spglib(positions, cell, types)
    println("Number of symmetries: ", get_nsymmetries(bcc_group))


    # Now check if we constrain the symmetry we get a null vector
    vector_input = rand(6)
    @show vector_input
    bcc_group.symmetrize_centroid!(vector_input)
    @show vector_input

    @test maximum(abs.(vector_input)) < 1e-12
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_bcc()
end
