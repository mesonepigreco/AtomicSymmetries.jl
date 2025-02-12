using AtomicSymmetries
using LinearAlgebra
using Test

function test_symmetrize_vector(; verbose=false)
    # Define the PbTe cell 
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.4 0.4 0.4]'

    # Get the symmetry group and ASR
    sym_group = get_symmetry_group_from_spglib(positions, cell, [1, 2])

    if verbose
        println("Found $(length(sym_group)) symmetries")
    end

    vector = randn(Float64, 6)
    asr! = ASRConstraint!(3)

    # Apply the acoustic sum rule
    asr!(vector)

    # Symmetrize in cartesian space
    symmetrize_vector!(vector, cell, sym_group)


    # Now the new vector
    cartesian_pos = zeros(Float64, size(positions)...)
    get_cartesian_coords!(cartesian_pos, positions, cell)
    cartesian_pos[:, 1] = - cartesian_pos[:, 2] / 2
    cartesian_pos[:, 2] = - cartesian_pos[:, 1] 

    # Check if the angle between the two vectors is 0 (or π)
    cosθ = reshape(cartesian_pos, :)' * vector / (norm(reshape(cartesian_pos, :)) * norm(vector))

    if verbose
        println("Vector after symmetrization: ", vector)
        println("Proposed vector:", cartesian_pos')
        println("Angle between the two vectors: ", cosθ)
    end
    @test abs(abs(cosθ) - 1) < 1e-8
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_symmetrize_vector(verbose=true)
end
