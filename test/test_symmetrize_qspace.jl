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
    types = [1, 2]

    # Create a 4x4x4 supercell
    supercell = [4, 4, 4]

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
                    R_lat[:, nat * (counter - 1) + h] = latvec  .* supercell
                    super_positions[:, nat * (counter - 1) + h] = positions[:, h] ./ supercell +  latvec

                    super_itau[nat * (counter - 1) + h] = h
                end

                counter += 1
            end
        end
    end

    for iat in 1:nat_sc
        super_types[iat] = types[super_itau[iat]]
    end
    
    if verbose
        println("Atoms:")
        for i in 1:nat_sc
            println(super_positions[:, i])
        end

        println()
        println("R lat:")
        for i in 1:nat_sc
            println(R_lat[:,i])
        end
        println()
        println("Q vec:")
        for i in 1:n_sc
            println(q_vec[:, i])
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
    i_sym = 3

    irt = zeros(Int, nat * n_sc)
    #TODO: It may not be the same symmetry (that is why it crashes). Check the translation
    AtomicSymmetries.get_irt!(irt, super_positions, uc_group.symmetries[i_sym], uc_group.translations[i_sym] ./ supercell)
    AtomicSymmetries.apply_sym_centroid!(u_next, u_coordinates[1, :], uc_group.symmetries[i_sym], ndims, irt)

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
    
    AtomicSymmetries.get_irt_q!(irt_q, q_vec, uc_group.symmetries[i_sym])
    @views AtomicSymmetries.apply_symmetry_vectorq!(q_next[1, :, :], q_coordinates[1, :, :], uc_group.symmetries[i_sym], uc_group.irt[i_sym], irt_q)

    AtomicSymmetries.vector_q2r!(u_next_back, q_next, q_vec, super_itau, R_lat)
 

    # Compare the two vectors
    if verbose 
        println("Test symmetries fourier")
        println()

        @show irt_q
        println()

        println("Original  | Transformed | Transformed Fourier")
        for i_at in 1:nat_sc
            println("$(u_coordinates[3*(i_at - 1)+1 : 3i_at])   |    $(u_next[3*(i_at-1)+1:3i_at])   |  $(u_next_back[1, 3*(i_at - 1)+1: 3i_at])")
        end
    end
    for i in 1:length(u_next_back)
        @test u_next_back[1, i] ≈ u_next[i] rtol = 1e-7
    end



    # APPLY SYMMETRY ON FC MATRIX
    fc_trial = randn(Float64, ndims*nat_sc, ndims*nat_sc)
    fc_trial .+= fc_trial'
    dynq_trial = zeros(Complex{Float64}, ndims*nat, ndims*nat, n_sc)
    fc_backward1 = zeros(Float64, ndims*nat_sc, ndims*nat_sc) 
    fc_backward2 = zeros(Float64, ndims*nat_sc, ndims*nat_sc) 
    dynq_back2 = zeros(Complex{Float64}, ndims*nat, ndims*nat, n_sc) 


    # Perform the fourier transform of the dynamical matrix
    AtomicSymmetries.matrix_r2q!(dynq_trial, fc_trial, q_vec, super_itau,
                                 R_lat)
    AtomicSymmetries.matrix_q2r!(fc_backward1, dynq_trial, q_vec, super_itau, R_lat)
    AtomicSymmetries.matrix_r2q!(dynq_back2, fc_backward1, q_vec, super_itau,
                                 R_lat)

    # Test the fourier transform of the dynamical matrix
    for iq in 1:n_sc
        for i in 1:ndims*nat
            for j in 1:ndims*nat
                @test dynq_trial[j, i, iq] ≈ dynq_back2[j, i, iq]
            end
        end
    end


    # TODO: Now we must test the symmetry application in Fourier space
    # Apply the symmetry in the dynq matrix
    dynq_back2 .= 0.0
    AtomicSymmetries.apply_symmetry_matrixq!(dynq_back2,
                                           dynq_trial,
                                           uc_group.symmetries[i_sym],
                                           uc_group.irt[i_sym],
                                           irt_q)
    AtomicSymmetries.matrix_q2r!(fc_trial, dynq_back2, q_vec, super_itau, R_lat)
    # Now, fc_trial contains the dynamical matrix with the symmetry applied in q space
    # Let us apply the symmetry also to fc_backward1 -> fc_backward2
    AtomicSymmetries.apply_sym_fc!(fc_backward2, fc_backward1, uc_group.symmetries[i_sym], ndims, irt)

    # Now we can compare fc_backward2 and fc_backward
    print_next = false
    if verbose && print_next
        println("Testing the symmetry application of the force constant matrices")
        println("Matrix")
        for i in 1:ndims
            println(uc_group.symmetries[i_sym][:, i])
        end
        @show irt_q
        @show irt
        println("Before | After symmetry (real)")
        for i in 1:ndims*nat_sc
            for j in 1:ndims*nat_sc
                print(fc_backward1[j, i] > 0  ? " " : "")
                print("$(round(fc_backward1[j, i]; digits=3)) ")
            end
            print("        ")
            for j in 1:ndims*nat_sc
                print(fc_backward2[j, i] > 0  ? " " : "")
                print("$(round(fc_backward2[j, i]; digits=3)) ")
            end
            println()
        end

        println()
        @show uc_group.irt[i_sym]
        println("Before | After symmetry (qspace)")
        for iq in 1:n_sc
            println("IQ = $iq")
            for i in 1:ndims*nat
                for j in 1:ndims*nat
                    print(real(dynq_trial[j, i, iq]) > 0 ? " " : "")
                    print("$(round(real(dynq_trial[j, i, iq]); digits=3)) ")
                end
                print("        ")
                for j in 1:ndims*nat_sc
                    print(real(dynq_back2[j, i, iq]) > 0 ? " " : "")
                    print("$(round(real(dynq_back2[j, i, iq]); digits=3)) ")
                end
                println()
            end
        end
    end

    for i in 1:ndims*nat_sc
        for j in 1:ndims*nat_sc
            if abs(fc_backward2[j, i]) < 1e-10
                @test abs(fc_trial[j, i]) < 1e-10
            else
                @test fc_backward2[j, i] ≈ fc_trial[j, i]
            end
        end
    end


end

if abspath(PROGRAM_FILE) == @__FILE__
    test_symmetrize_q_space(; verbose=true)
end
