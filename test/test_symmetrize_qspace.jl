using Test
using AtomicSymmetries

function test_symmetrize_q_space(; verbose=false)
    a = 2.87
    b = 2.93
    c = 2.60
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    # positions = zeros(Float64, 3, 1)
    # positions .= [0.0 
    #              0.0 
    #              0.0 ]
    # types = [1]

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
    translations = AtomicSymmetries.get_translations(sc_group)

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


    rmat = randn(Float64, ndims*nat_sc, ndims*nat_sc)
    fc_trial = zeros(Float64, ndims*nat_sc, ndims*nat_sc)
    dynq_trial = zeros(Complex{Float64}, ndims*nat, ndims*nat, n_sc)
    fc_backward1 = zeros(Float64, ndims*nat_sc, ndims*nat_sc) 
    fc_backward2 = zeros(Float64, ndims*nat_sc, ndims*nat_sc) 
    dynq_back2 = zeros(Complex{Float64}, ndims*nat, ndims*nat, n_sc) 


    # APPLY SYMMETRY ON FC MATRIX
    print_next = true #n_sc == 1
    for i_sym in 1:length(uc_group)
        if verbose
            println()
            println("-----------------------")
            println("Testing the symmetry $i_sym:")
            println("-----------------------")
            for k in 1:ndims
                println(uc_group.symmetries[i_sym][k, :])
            end
            @show uc_group.translations[i_sym]
        end

        dynq_trial .= 0
        fc_backward1 .= 0
        fc_backward2 .= 0
        dynq_back2 .= 0
        fc_trial .= rmat
        fc_trial .+= fc_trial'

        # Update irt
        AtomicSymmetries.get_irt!(irt, super_positions, uc_group.symmetries[i_sym], uc_group.translations[i_sym] ./ supercell)

        
        # Perform the fourier transform of the dynamical matrix
        AtomicSymmetries.matrix_r2q!(dynq_trial, fc_trial, q_vec, super_itau,
                                     R_lat)
        AtomicSymmetries.matrix_q2r!(fc_backward1, dynq_trial, q_vec, super_itau, R_lat; translations=translations)
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
        AtomicSymmetries.matrix_q2r!(fc_trial, dynq_back2, q_vec, super_itau, R_lat; translations=translations)
        # Now, fc_trial contains the dynamical matrix with the symmetry applied in q space
        # Let us apply the symmetry also to fc_backward1 -> fc_backward2
        AtomicSymmetries.apply_sym_fc!(fc_backward2, fc_backward1, uc_group.symmetries[i_sym], ndims, irt)

        # Now we can compare fc_backward2 and fc_backward
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
                    for j in 1:ndims*nat
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
    # At this point, the application of individual symmetry operation works in fourier space.
    # However, we must still check if the full symmetrization works


    # Now, we can test the full symmetrization
    # TODO: THe errpr see,s tp be that fc_backward1 has a lot of zeros and 
    # it is not the real starting point!!!
    q_symmetries = AtomicSymmetries.SymmetriesQSpace(uc_group, q_vec)
    AtomicSymmetries.symmetrize_matrix_q!(dynq_back2, dynq_trial, q_symmetries)
    AtomicSymmetries.matrix_q2r!(fc_trial, dynq_back2, q_vec, super_itau, R_lat; translations=translations)

    # Now let us perform the symmetrization directly in cartesian space
    AtomicSymmetries.matrix_q2r!(fc_backward2, dynq_trial, q_vec, super_itau, R_lat; translations=translations)
    fc_backward1 = copy(fc_backward2)
    sc_group.symmetrize_fc!(fc_backward2)

    if verbose && print_next
        println("Testing the symmetrization of the force constant matrices")
        println("Before | After symmetriezation (real)")
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
        println("Before | After symmetrization (qspace)")
        for iq in 1:n_sc
            println("IQ = $iq")
            for i in 1:ndims*nat
                for j in 1:ndims*nat
                    print(real(dynq_trial[j, i, iq]) > 0 ? " " : "")
                    print("$(round(real(dynq_trial[j, i, iq]); digits=3)) ")
                end
                print("        ")
                for j in 1:ndims*nat
                    print(real(dynq_back2[j, i, iq]) > 0 ? " " : "")
                    print("$(round(real(dynq_back2[j, i, iq]); digits=3)) ")
                end
                println()
            end
        end
    end


    # Test the comparison between the two
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


function test_fourier_matrix(; verbose=false)
    # Define the values
    q_tot = [0.0 0.5; 0.0 0.0; 0.0 0.0]
    R_lat = [0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    itau = [1, 2, 1, 2]
    translations = [[1, 2, 3, 4], [3, 4, 1, 2]]

    
    # Why this does not pass???
    phi_r = [0.000438140025113431 -3.253043901411059e-19 2.8916214609205817e-20 -2.7715586166034142e-5 -3.419413084247844e-19 5.693724722070692e-22 0.00035295127968299085 3.993258130969505e-20 -2.927765786020131e-20 -2.771558616603443e-5 -6.150097785543722e-19 -2.969132245658229e-21; -3.253043901411059e-19 0.0003064298227686857 -1.2939857906710487e-20 -7.096425988974665e-20 9.023040681526278e-6 -4.118508929857172e-21 -3.2908603145680523e-19 -0.000305254564575895 1.790570273163483e-20 1.6161493113436632e-19 9.023040681526395e-6 -3.462746986889771e-21; 2.8916214609205817e-20 -1.2939857906710487e-20 5.487582086902722e-5 2.051602034199394e-20 -1.7816559862726685e-20 -5.867451475292805e-5 6.846883744848607e-20 2.0225650295352952e-20 5.4875820869027176e-5 -5.2682927198537854e-20 1.004279935668362e-19 -5.867451475292839e-5; -2.7715586166034142e-5 -7.096425988974665e-20 2.051602034199394e-20 0.00029620051951981507 -1.6207173971121348e-18 -6.75152186267878e-20 -2.7715586166034058e-5 5.655596364714862e-20 3.7337131152926256e-20 2.540966649977185e-5 -2.116169389430518e-18 5.812234226403181e-21; -3.419413084247844e-19 9.023040681526278e-6 -1.7816559862726685e-20 -1.6207173971121348e-18 0.00021631370598860212 -2.999512975092786e-20 -5.439640291776784e-19 9.023040681526278e-6 2.4053157488047915e-20 -1.5820480880547106e-18 6.078374762586847e-5 8.115284513570756e-21; 5.693724722070692e-22 -4.118508929857172e-21 -5.867451475292805e-5 -6.75152186267878e-20 -2.999512975092786e-20 0.00019620902014832575 9.03021689436743e-22 -3.872384176705582e-21 -5.867451475292801e-5 -8.763839013210105e-20 -7.228196745223495e-20 -7.073668544892055e-5; 0.00035295127968299085 -3.2908603145680523e-19 6.846883744848607e-20 -2.7715586166034058e-5 -5.439640291776784e-19 9.03021689436743e-22 0.0004381400251134308 6.081283811833534e-20 -7.833719700277762e-20 -2.771558616603424e-5 -2.198982237621582e-19 1.657123362265235e-21; 3.993258130969505e-20 -0.000305254564575895 2.0225650295352952e-20 5.655596364714862e-20 9.023040681526278e-6 -3.872384176705582e-21 6.081283811833534e-20 0.00030642982276868607 -1.6298718233991108e-20 -2.17446117758632e-19 9.023040681526146e-6 -2.5980418476214058e-21; -2.927765786020131e-20 1.790570273163483e-20 5.4875820869027176e-5 3.7337131152926256e-20 2.4053157488047915e-20 -5.867451475292801e-5 -7.833719700277762e-20 -1.6298718233991108e-20 5.487582086902714e-5 1.6609235481087098e-19 -2.863977938970003e-20 -5.8674514752928346e-5; -2.771558616603443e-5 1.6161493113436632e-19 -5.2682927198537854e-20 2.540966649977185e-5 -1.5820480880547106e-18 -8.763839013210105e-20 -2.771558616603424e-5 -2.17446117758632e-19 1.6609235481087098e-19 0.00029620051951981615 -2.025781547051897e-18 -1.9066109743533083e-20; -6.150097785543722e-19 9.023040681526395e-6 1.004279935668362e-19 -2.116169389430518e-18 6.078374762586847e-5 -7.228196745223495e-20 -2.198982237621582e-19 9.023040681526146e-6 -2.863977938970003e-20 -2.025781547051897e-18 0.00021631370598860207 -5.902673655019174e-20; -2.969132245658229e-21 -3.462746986889771e-21 -5.867451475292839e-5 5.812234226403181e-21 8.115284513570756e-21 -7.073668544892055e-5 1.657123362265235e-21 -2.5980418476214058e-21 -5.8674514752928346e-5 -1.9066109743533083e-20 -5.902673655019174e-20 0.00019620902014832553]
    
    phi_q = zeros(Complex{Float64}, 6, 6, 2) 
    phi_q2r = similar(phi_r)
    phi_q2r .= 0
    phi_q_bis = similar(phi_q)

    # Convert in q space
    matrix_r2q!(phi_q, phi_r, q_tot, itau, R_lat)

    # Convert back in r space
    matrix_q2r!(phi_q2r, phi_q, q_tot, itau, R_lat;
                translations)

    matrix_r2q!(phi_q_bis, phi_q2r, q_tot, itau, R_lat)

    for iq in 1:2
        for k in 1:6
            for h in 1:6
                @test phi_q_bis[h, k, iq] ≈ phi_q[h, k, iq]
            end
        end
    end
    if verbose
        println("Q pass the test")
    end

    # Check if they are equal
    for i in 1:12
        for j in 1:12
            @test phi_q2r[i, j] ≈ phi_r[i, j]
        end
    end
    if verbose
        println("R pass the test")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_fourier_matrix(; verbose=true)
    test_symmetrize_q_space(; verbose=false)
end
