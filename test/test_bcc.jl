using AtomicSymmetries
using Test

function test_bcc(;verbose=false)
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
    if verbose
        println("Number of symmetries: ", get_nsymmetries(bcc_group))
    end


    # Now check if we constrain the symmetry we get a null vector
    vector_input = rand(6)
    @show vector_input
    bcc_group.symmetrize_centroid!(vector_input) # This is good for crystal coordinates
    #symmetrize_vector!(vector_input, bcc_group, cell) # This is good for cartesian coordinates
    @show vector_input

    @test maximum(abs.(vector_input)) < 1e-12

    # add the noise to the position and check if 
    # applying the translations the vector remains the same
    cartesian_positions = zeros(Float64, size(positions)...)
    get_cartesian_coords!(cartesian_positions, positions, cell)
    
    new_positions = cartesian_positions + 0.1 * rand(size(positions)...)
    symmetrize_positions!(new_positions, cell, bcc_group)
    if verbose
        println("Positions: ", cartesian_positions)
        println("New positions: ", new_positions)
    end

    for i in 1:size(positions, 2)
        @test maximum(abs.(cartesian_positions[:, i] - new_positions[:, i])) < 1e-12
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_bcc(verbose=true)
end
