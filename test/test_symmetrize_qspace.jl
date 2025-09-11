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
    phi_r = [6.944077713497171e-9 -3.7725100300087146e-24 -2.7916441566992746e-24 -2.6502766015545497e-9 -3.549208423473465e-24 8.664561479838025e-25 -1.64352451038805e-9 1.431307838807254e-25 2.4957276696667897e-24 -2.6502766015545567e-9 3.796838388233213e-24 -1.0515594427328438e-24; -3.7725100300087146e-24 1.0247441921697906e-8 -3.1270747218558373e-24 2.775336528332932e-24 -3.52850276853492e-9 -1.239793161881364e-23 -4.8300034362667475e-24 -3.1904363846278865e-9 -5.78965596508118e-24 1.7438715494390113e-24 -3.528502768535103e-9 1.6466259623049095e-23; -2.7916441566992746e-24 -3.1270747218558373e-24 1.1692015910447308e-8 -1.2896862337118913e-24 -3.66031626862765e-24 -2.5527304141603894e-9 1.9423782165784355e-24 -1.2451775904376244e-24 -6.5865550821265495e-9 -3.4492648217963074e-24 8.484309948448528e-24 -2.5527304141603745e-9; -2.6502766015545497e-9 2.775336528332932e-24 -1.2896862337118913e-24 6.745007606626527e-9 -3.365587283602352e-24 -1.3395244122046136e-24 -2.6502766015545576e-9 6.2068579361313445e-24 3.516698839572828e-24 -1.4444544035174181e-9 -4.6071021130034915e-24 -5.164822483403609e-24; -3.549208423473465e-24 -3.52850276853492e-9 -3.66031626862765e-24 -3.365587283602352e-24 1.0692874210984586e-8 1.8481131247656558e-23 1.6762961146945018e-24 -3.528502768534917e-9 9.575401897091406e-24 5.047160113039026e-24 -3.635868673914753e-9 -2.6335213548191474e-23; 8.664561479838025e-25 -1.239793161881364e-23 -2.5527304141603894e-9 -1.3395244122046136e-24 1.8481131247656558e-23 1.0537531217909593e-8 7.048074030485992e-25 -8.776504679457027e-24 -2.5527304141603803e-9 2.598641203715694e-24 -5.847207241187604e-25 -5.432070389588828e-9; -1.64352451038805e-9 -4.8300034362667475e-24 1.9423782165784355e-24 -2.6502766015545576e-9 1.6762961146945018e-24 7.048074030485992e-25 6.944077713497157e-9 3.652955011663864e-24 -2.6772896034501183e-24 -2.6502766015545563e-9 2.8595707022475537e-24 -3.7211657035649032e-25; 1.431307838807254e-25 -3.1904363846278865e-9 -1.2451775904376244e-24 6.2068579361313445e-24 -3.528502768534917e-9 -8.776504679457027e-24 3.652955011663864e-24 1.024744192169784e-8 -1.0498413543432898e-23 -1.1074661113428624e-23 -3.5285027685350364e-9 1.8238540043773188e-23; 2.4957276696667897e-24 -5.78965596508118e-24 -6.5865550821265495e-9 3.516698839572828e-24 9.575401897091406e-24 -2.5527304141603803e-9 -2.6772896034501183e-24 -1.0498413543432898e-23 1.1692015910447286e-8 1.4782376420899345e-25 6.156807344352184e-24 -2.5527304141603724e-9; -2.6502766015545567e-9 1.7438715494390113e-24 -3.4492648217963074e-24 -1.4444544035174181e-9 5.047160113039026e-24 2.598641203715694e-24 -2.6502766015545563e-9 -1.1074661113428624e-23 1.4782376420899345e-25 6.745007606626523e-9 3.2970552721212182e-24 5.8633744795800024e-24; 3.796838388233213e-24 -3.528502768535103e-9 8.484309948448528e-24 -4.6071021130034915e-24 -3.635868673914753e-9 -5.847207241187604e-25 2.8595707022475537e-24 -3.5285027685350364e-9 6.156807344352184e-24 3.2970552721212182e-24 1.06928742109849e-8 -4.98744144435483e-24; -1.0515594427328438e-24 1.6466259623049095e-23 -2.5527304141603745e-9 -5.164822483403609e-24 -2.6335213548191474e-23 -5.432070389588828e-9 -3.7211657035649032e-25 1.8238540043773188e-23 -2.5527304141603724e-9 5.8633744795800024e-24 -4.98744144435483e-24 1.0537531217909603e-8]
    
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
                @test phi_q_bis[h, k, iq] ≈ phi_q[h, k, iq] atol = 1e-14
            end
        end
    end
    if verbose
        println("Q pass the test")
    end

    # Check if they are equal
    for i in 1:12
        for j in 1:12
            @test phi_q2r[i, j] ≈ phi_r[i, j]  atol = 1e-14
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
