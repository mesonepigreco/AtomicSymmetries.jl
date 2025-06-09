using AtomicSymmetries


function test_symmetrize_q_space(; verbose=false)
    a = 2.87
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    positions = [0.0 0.6
                 0.0 0.6
                 0.0 0.6]

    # Create a 4x4x4 supercell
    supercell = [4, 4, 4]

    nat = size(positions, 2)
    ndims = size(positions, 1)

    n_sc = prod(supercell)

    # Prepare the R and q space
    R_lat = zeros(Float64, ndims, n_sc)
    q_vec = zeros(Float64, ndims, n_sc)

    super_cell = zeros(Float64, 3, 3)
    super_positions = zeros(Float64, ndims, nat * n_sc)
    super_types = ones(Int, nat*n_sc)

    for i in 1:ndims
        @views super_cell[:, i] .= cell[:, i] * supercell[i]
    end

    counter = 1
    for i in 1:supercell[1]
        for j in 1:supercell[2]
            for k in 1:supercell[3]
                R_lat[:, counter] .= [i-1, j-1, k-1]
                @views q_vec[:, counter] .= R_lat[:, counter]
                q_vec[:, counter] ./= supercell

                for h in 1:nat
                    super_positions[:, nat * (counter - 1) + h] = positions[:, h] + R_lat[:, counter]
                end

                counter += 1
            end
        end
    end

    # Supercell symmetry group
    sc_group = get_symmetry_group_from_spglib(super_positions, super_cell, super_types)

    u_coordinates = randn(Float64, nat * ndims * n_sc)




end

if abspath(PROGRAM_FILE) == @__FILE__
    test_symmetrize_q_space(; verboese=true)
end
