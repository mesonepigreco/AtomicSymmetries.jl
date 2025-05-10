using Test
using AtomicSymmetries
using SparseArrays

function test_symmetry_sparsification(; verbose=false)
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

    # Get the sparsified representation of the symmetry group
    sparse_symmetry_group = sparse(symmetry_group)

    # Generate a random vector, and check that the application of symmetries is the same with
    # the sparse module and with the exact module
    initial_vector = randn(Float64, length(atomic_positions_cartesian))
    final_vector = zeros(Float64, length(atomic_positions_cartesian))
    sparse_final_vector = zeros(Float64, length(atomic_positions_cartesian))
    for i in 1:length(symmetry_group)
        final_vector .= 0.0

        # Apply the symmetry to the vector (standard)
        AtomicSymmetries.apply_sym_centroid!(final_vector, initial_vector, symmetry_group.symmetries[i], 3, symmetry_group.irt[i])
        sparse_final_vector .= 0.0
        # Apply the symmetry to the vector (sparse)
        sparse_final_vector = (initial_vector' * sparse_symmetry_group.symmetries[i])'

        

        # Check that the two vectors are the same
        if verbose
            println()
            println("Symmetry $i")
            @show symmetry_group.symmetries[i]
            @show sparse_symmetry_group.symmetries[i]
            println("Initial vector: ", initial_vector)
            println("Final vector: ", final_vector)
            println("Sparse final vector: ", sparse_final_vector)
            println(typeof(sparse_final_vector))
            println("Difference: ", abs.(final_vector - sparse_final_vector))
            println()
        end
        @test maximum(abs.(final_vector - sparse_final_vector)) < 1e-8

    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_symmetry_sparsification(; verbose=true)
end
 
