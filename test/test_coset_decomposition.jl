using AtomicSymmetries
using LinearAlgebra
using Test

if abspath(PROGRAM_FILE) == @__FILE__
    include("define_cell.jl")
end


function test_coset_decomposition(; verbose=false)
    positions, cell, types = get_pm3m_perovskite()
    original_group = get_symmetry_group_from_spglib(positions, cell, types)

    n_original = length(original_group)
    if verbose
        println("Number of original symmetries: ", n_original)
    end

    # Build CosetSymmetryGroup with [1,0,0] perturbation
    perturbation = [1.0, 0.0, 0.0]
    coset_group = CosetSymmetryGroup(original_group, perturbation, cell)

    n_subgroup = length(coset_group.subgroup)
    n_cosets = length(coset_group.cosets)

    if verbose
        println("Number of subgroup symmetries: ", n_subgroup)
        println("Number of coset representatives: ", n_cosets)
    end

    # Test 1: |subgroup| * |cosets| == |original_group|
    @test n_subgroup * n_cosets == n_original

    # Test 2: The subgroup should have 8 symmetries (C4v for cubic with [1,0,0] perturbation)
    @test n_subgroup == 8
    @test n_cosets == 6  # 48/8

    # Test 3: The subgroup is invariant under the perturbation
    # (filtering it again should not remove any element)
    subgroup_copy = deepcopy(coset_group.subgroup)
    filter_invariant_symmetries!(subgroup_copy, perturbation, cell)
    @test length(subgroup_copy) == n_subgroup

    # Test 4: Each coset representative combined with the subgroup
    # reconstructs elements of the original group (right cosets: h_j * g_k)
    n_found = 0
    for k in 1:n_cosets
        g_k = coset_group.cosets.symmetries[k]
        for j in 1:n_subgroup
            h_j = coset_group.subgroup.symmetries[j]
            composed = h_j * g_k

            # Check this is in the original group
            found = false
            for i in 1:n_original
                if maximum(abs.(composed .- original_group.symmetries[i])) < 1e-6
                    found = true
                    break
                end
            end
            @test found
            if found
                n_found += 1
            end
        end
    end
    @test n_found == n_original

    # Test 5: All cosets are disjoint (no element of G appears in two different cosets)
    all_composed = Matrix{Float64}[]
    for k in 1:n_cosets
        g_k = coset_group.cosets.symmetries[k]
        for j in 1:n_subgroup
            h_j = coset_group.subgroup.symmetries[j]
            push!(all_composed, h_j * g_k)
        end
    end
    for i in 1:length(all_composed)
        for j in i+1:length(all_composed)
            @test maximum(abs.(all_composed[i] .- all_composed[j])) > 1e-6
        end
    end

    # Test 6: irt is properly initialized for cosets
    nat = original_group.n_particles
    for k in 1:n_cosets
        @test length(coset_group.cosets.irt[k]) == nat
        # irt values should be valid atom indices
        for idx in coset_group.cosets.irt[k]
            @test 1 <= idx <= nat
        end
    end
end


function test_coset_decomposition_qspace(; verbose=false)
    positions, cell, types = get_pm3m_perovskite()
    original_group = get_symmetry_group_from_spglib(positions, cell, types)

    # Create a q-point grid (2x2x2 for a cubic cell)
    q_points = Float64[0 0.5 0 0 0.5 0.5 0 0.5;
                       0 0 0.5 0 0.5 0 0.5 0.5;
                       0 0 0 0.5 0 0.5 0.5 0.5]
    q_group = SymmetriesQSpace(original_group, q_points)

    n_original = length(q_group)
    if verbose
        println("Q-space: Number of original symmetries: ", n_original)
    end

    perturbation = [1.0, 0.0, 0.0]
    coset_group = CosetSymmetryGroup(q_group, perturbation, cell)

    n_subgroup = length(coset_group.subgroup)
    n_cosets = length(coset_group.cosets)

    if verbose
        println("Q-space: Number of subgroup symmetries: ", n_subgroup)
        println("Q-space: Number of coset representatives: ", n_cosets)
    end

    # Test 1: |subgroup| * |cosets| == |original_group|
    @test n_subgroup * n_cosets == n_original

    # Test 2: irt_q is properly initialized for both subgroup and cosets
    nq = size(q_points, 2)
    for k in 1:n_cosets
        @test length(coset_group.cosets.irt_q[k]) == nq
        for idx in coset_group.cosets.irt_q[k]
            @test 1 <= idx <= nq
        end
    end
    for k in 1:n_subgroup
        @test length(coset_group.subgroup.irt_q[k]) == nq
        for idx in coset_group.subgroup.irt_q[k]
            @test 1 <= idx <= nq
        end
    end

    # Test 3: The real-space part is consistent
    @test length(coset_group.subgroup.symmetries) == n_subgroup
    @test length(coset_group.cosets.symmetries) == n_cosets
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_coset_decomposition(verbose=true)
    test_coset_decomposition_qspace(verbose=true)
end
