using AtomicSymmetries
using LinearAlgebra
using Test

if abspath(PROGRAM_FILE) == @__FILE__
    include("define_cell.jl")
end


function test_filter(; verbose=false)
    positions, cell, types = get_pm3m_perovskite()
    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)

    n_before_filter = length(symmetry_group)
    if verbose
        println("Number of original symmetries: ", get_nsymmetries(symmetry_group))
    end

    # Now filter the symmetries
    filter_invariant_symmetries!(symmetry_group, [1.0, 0.0, 0.0])

    n_after_filter = length(symmetry_group)

    if verbose
        println("Number of symmetries after filter: ", get_nsymmetries(symmetry_group))
    end

    # Try extending the symmetry group
    complete_symmetry_group!(symmetry_group)

    n_after_complete = length(symmetry_group)

    if verbose
        println("Number of symmetries after complete: ", get_nsymmetries(symmetry_group))
    end

    @test n_before_filter > n_after_filter
    @test n_after_filter == n_after_complete
    @test n_after_complete == 8

    # Check that the inversion is not in the symmetry group
    for i in 1:length(symmetry_group)
        sym = symmetry_group.symmetries[i]
        if verbose
            println("Symmetry $i: ", sym)
        end
        @test max(abs.(sym .+ I(3))...) > 1e-12
    end
end

function test_filter_full_vector(; verbose=false)
    positions, cell, types = get_pm3m_perovskite()
    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)

    n_before_filter = length(symmetry_group)
    if verbose
        println("Number of original symmetries: ", get_nsymmetries(symmetry_group))
    end

    # Now filter the symmetries
    vector_filter = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    filter_invariant_symmetries!(symmetry_group, vector_filter)

    n_after_filter = length(symmetry_group)

    if verbose
        println("Number of symmetries after filter: ", get_nsymmetries(symmetry_group))
    end

    # Try extending the symmetry group
    complete_symmetry_group!(symmetry_group)

    n_after_complete = length(symmetry_group)

    if verbose
        println("Number of symmetries after complete: ", get_nsymmetries(symmetry_group))
    end

    @test n_before_filter > n_after_filter
    @test n_after_filter == n_after_complete

    @test n_after_complete == 8

    # Check that the inversion is not in the symmetry group
    for i in 1:length(symmetry_group)
        sym = symmetry_group.symmetries[i]
        if verbose
            println("Symmetry $i: ", sym)
        end
        @test max(abs.(sym .+ I(3))...) > 1e-12
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    test_filter(verbose=true)
    test_filter_full_vector(verbose=true)
end

