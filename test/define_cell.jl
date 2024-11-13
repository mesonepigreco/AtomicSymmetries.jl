function get_pm3m_perovskite()
    positions = transpose([0.0 0.0 0.0
                 0.5 0.5 0.5
                 0.5 0.5 0.0
                 0.0 0.5 0.5
                 0.5 0.0 0.5])
    cell = [1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0]
    types = [1, 2, 3, 3, 3]
    return (positions, cell, types)
end


function get_gold_unit_cell()
    cell = [2.8900186179635647     0.0000000000000000     0.0000000000000000
        1.4450093089817819     2.5028295405664411     0.0000000000000000
        1.4450093089817826     0.8342765135221475     2.3599213882279275]

    position = zeros(Float64, 3, 1)
    position[:, 1] = [0.0000000000000000, 0.0000000000000000, 0.0000000000000000]
    types = [1]
    return (position, cell', types)
end

function get_pm3m_supercell(; supercell=[2, 2, 2])
    uc_positions, cell, types = get_pm3m_perovskite()
    get_supercell(uc_positions, cell, types, supercell)
end

function get_supercell(uc_positions, cell, types, supercell)
    # Now create a supercell
    sup_dim = prod(supercell)
    nat_uc = size(uc_positions, 2)
    positions = zeros(3, sup_dim*nat_uc)
    types = zeros(Int, sup_dim*nat_uc)
    new_cell = similar(cell)
    for i_x in 1:supercell[1]
        for i_y in 1:supercell[2]
            for i_z in 1:supercell[3]
                for k in 1:nat_uc
                    cell_index = begin 
                        (i_x-1)*supercell[3]*supercell[2]*nat_uc + 
                        (i_y-1)*supercell[3]*nat_uc + 
                        (i_z - 1)*nat_uc
                    end
                    @views positions[:, cell_index + k] = uc_positions[:, k] + (i_x-1)*cell[:, 1] + (i_y-1)*cell[:, 2] + (i_z-1)*cell[:, 3]
                    types[cell_index + k] = types[k]
                end
            end
        end
    end

    
    # for i in 1:3
    #     @views positions[i, :] ./= supercell[i]
    # end

    for i in 1:3
        for j in 1:3
            new_cell[i, j] = cell[i, j]*supercell[j]
        end
    end

    crystal = similar(positions)
    get_crystal_coords!(crystal, positions, new_cell)

    return (crystal, new_cell, types)
end
    

