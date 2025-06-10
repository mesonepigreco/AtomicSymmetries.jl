using Test
using AtomicSymmetries

function test_symmetrize_q_space(; verbose=false)
    a = 2.87
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    positions = [0.0 0.6
                 0.0 0.6
                 0.0 0.6]
    types = ones(Int, 2)

    # Create a 4x4x4 supercell
    supercell = [1, 1, 1]

    nat = size(positions, 2)
    ndims = size(positions, 1)

    n_sc = prod(supercell)
    nat_sc = nat*n_sc

    # Prepare the R and q space
    R_lat = zeros(Float64, ndims, nat*n_sc)
    q_vec = zeros(Float64, ndims, n_sc)

    super_cell = zeros(Float64, 3, 3)
    super_positions = zeros(Float64, ndims, nat * n_sc)
    super_types = ones(Int, nat*n_sc)
    super_itau = zeros(Int, nat*n_sc)
    
    for i in 1:ndims
        @views super_cell[:, i] .= cell[:, i] * supercell[i]
    end

    counter = 1
    for i in 1:supercell[1]
        for j in 1:supercell[2]
            for k in 1:supercell[3]
                latvec = [i-1, j-1, k-1] ./ supercell
                @views q_vec[:, counter] .= latvec

                for h in 1:nat
                    R_lat[:, nat * (counter - 1) + h] = latvec 
                    super_positions[:, nat * (counter - 1) + h] = positions[:, h] ./ supercell +  latvec

                    super_itau[nat * (counter - 1) + h] = h
                end

                counter += 1
            end
        end
    end
    
    if verbose
        println("Atoms:")
        for i in 1:nat_sc
            println(super_positions[:, i])
        end
    end

    # Supercell symmetry group
    uc_group = get_symmetry_group_from_spglib(positions, cell, types)
    sc_group = get_symmetry_group_from_spglib(super_positions, super_cell, super_types)

    if verbose
        println("The symmetries : $(length(uc_group))")
        for i in 1:length(uc_group)
            println("S:")
            @show uc_group.symmetries[i]'
            @show uc_group.translations[i]
        end
    end

    n_rand = 1
    u_coordinates = randn(Float64, n_rand, nat * ndims * n_sc)

    # Convert to q-space
    q_coordinates = zeros(Complex{Float64}, 1, nat*ndims, n_sc)
    u_next = zeros(Float64, nat*ndims*n_sc)
    q_next = zeros(Complex{Float64}, 1, nat*ndims, n_sc)
    u_next_back = zeros(Float64, 1, nat*ndims*n_sc)
    irt_q = zeros(Int, n_sc)

    #
    # Apply the symmetry in real space

    irt = zeros(Int, nat * n_sc)
    AtomicSymmetries.get_irt!(irt, super_positions, uc_group.symmetries[2], uc_group.translations[2] ./ supercell)
    AtomicSymmetries.apply_sym_centroid!(u_next, u_coordinates[1, :], uc_group.symmetries[2], ndims, irt)

    # Apply the symmetry in q-space
    AtomicSymmetries.vector_r2q!(q_coordinates, u_coordinates, q_vec, super_itau, R_lat)

    #@views AtomicSymmetries.apply_symmetry_vectorq!(q_next[1, :, :], q_coordinates[1, :, :], uc_group.symmetries[2], uc_group.irt[2], irt_q)

    AtomicSymmetries.vector_q2r!(u_next_back, q_coordinates, q_vec, super_itau, R_lat)
    
    # Test the fourier transform
    if verbose
        println("Test fourier transform back")

        println()
        @show u_coordinates[1, :]
        @show q_coordinates[1, :, :]'
        @show u_next_back[1, :]
    end
    for i in 1:length(u_next_back)
        @test u_coordinates[1, i] ≈ u_next_back[1, i] rtol = 1e-7
    end
    
    AtomicSymmetries.get_irt_q!(irt_q, q_vec, uc_group.symmetries[2])
    @views AtomicSymmetries.apply_symmetry_vectorq!(q_next[1, :, :], q_coordinates[1, :, :], uc_group.symmetries[2], uc_group.irt[2], irt_q)

    AtomicSymmetries.vector_q2r!(u_next_back, q_next, q_vec, super_itau, R_lat)
 

    # Compare the two vectors
    if verbose 
        println("Test symmetries fourier")
    end
    for i in 1:length(u_next_back)
        @test u_next_back[1, i] ≈ u_next[i] rtol = 1e-7
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_symmetrize_q_space(; verbose=true)
end
