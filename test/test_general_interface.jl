using AtomicSymmetries
using LinearAlgebra
using Test


function test_rotate_vector_real(; verbose=false)
    # Use PbTe unit cell (non-orthogonal, 2 atoms)
    a = 12.21
    cell = collect([-a 0.0 a; 0.0 a a; -a a 0.0]')
    positions = collect([0.0 0.0 0.0; 0.4 0.4 0.4]')

    sym_group = get_symmetry_group_from_spglib(positions, cell, [1, 2])
    n_sym = get_nsymmetries(sym_group)
    n_dims = 3
    n_atoms = 2
    n_modes = n_dims * n_atoms

    reciprocal_vectors = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_vectors, cell)

    # Create a random vector and symmetrize via rotate_vector! + averaging
    old_vector = randn(Float64, n_modes)
    asr! = ASRConstraint!(n_dims)
    asr!(old_vector)

    avg_vector = zeros(Float64, n_modes)
    new_vector = zeros(Float64, n_modes)
    for i in 1:n_sym
        new_vector .= 0
        rotate_vector!(new_vector, old_vector, cell, reciprocal_vectors, sym_group, i)
        avg_vector .+= new_vector
    end
    avg_vector ./= n_sym

    # Compare with the existing symmetrize_vector!
    ref_vector = copy(old_vector)
    symmetrize_vector!(ref_vector, cell, sym_group)

    if verbose
        println("rotate_vector avg: ", avg_vector)
        println("symmetrize_vector: ", ref_vector)
        println("diff: ", maximum(abs.(avg_vector - ref_vector)))
    end

    @test isapprox(avg_vector, ref_vector; atol=1e-10)
end


function test_rotate_matrix_real(; verbose=false)
    # Use a cubic BCC cell
    a = 2.87
    cell = [a 0.0 0.0
            0.0 a 0.0
            0.0 0.0 a]
    positions = [0.0 0.5
                 0.0 0.5
                 0.0 0.5]
    types = [1, 1]

    sym_group = get_symmetry_group_from_spglib(positions, cell, types)
    n_sym = get_nsymmetries(sym_group)

    reciprocal_vectors = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_vectors, cell)

    # Create a random 3x3 matrix (stress tensor)
    old_matrix = randn(Float64, 3, 3)
    old_matrix .= (old_matrix .+ old_matrix') ./ 2  # Make symmetric

    avg_matrix = zeros(Float64, 3, 3)
    new_matrix = zeros(Float64, 3, 3)
    for i in 1:n_sym
        new_matrix .= 0
        rotate_matrix!(new_matrix, old_matrix, cell, reciprocal_vectors, sym_group, i)
        avg_matrix .+= new_matrix
    end
    avg_matrix ./= n_sym

    # For a cubic cell, the symmetrized stress tensor should be proportional to identity
    trace_val = tr(avg_matrix) / 3
    expected = trace_val * I(3)

    if verbose
        println("Symmetrized stress tensor:")
        println(avg_matrix)
        println("Expected (proportional to I): ", trace_val, " * I")
        println("diff: ", maximum(abs.(avg_matrix - expected)))
    end

    @test isapprox(avg_matrix, expected; atol=1e-10)
end


function test_rotate_dynamical_matrix_real(; verbose=false)
    # Use PbTe unit cell
    a = 12.21
    cell = collect([-a 0.0 a; 0.0 a a; -a a 0.0]')
    positions = collect([0.0 0.0 0.0; 0.4 0.4 0.4]')

    sym_group = get_symmetry_group_from_spglib(positions, cell, [1, 2])
    n_sym = get_nsymmetries(sym_group)
    n_dims = 3
    n_atoms = 2
    n_modes = n_dims * n_atoms

    reciprocal_vectors = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_vectors, cell)

    # Create a random force constant matrix
    old_fc = randn(Float64, n_modes, n_modes)
    old_fc .= (old_fc .+ old_fc') ./ 2  # Make symmetric

    avg_fc = zeros(Float64, n_modes, n_modes)
    new_fc = zeros(Float64, n_modes, n_modes)
    for i in 1:n_sym
        new_fc .= 0
        rotate_dynamical_matrix!(new_fc, old_fc, cell, reciprocal_vectors, sym_group, i)
        avg_fc .+= new_fc
    end
    avg_fc ./= n_sym

    # Compare with the existing symmetrize_fc!
    ref_fc = copy(old_fc)
    symmetrize_fc!(ref_fc, cell, sym_group)

    if verbose
        println("diff: ", maximum(abs.(avg_fc - ref_fc)))
    end

    @test isapprox(avg_fc, ref_fc; atol=1e-10)
end


function test_rotate_vector_qspace(; verbose=false)
    # Use a tetragonal cell with 2 atoms (lower symmetry → non-trivial symmetrized vector)
    a = 3.0
    c = 4.5
    unit_cell = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 c]
    positions = [0.0 0.5; 0.0 0.5; 0.0 0.3]  # atom 2 at general z → non-trivial result
    types = [1, 2]

    symmetry_group_uc = get_symmetry_group_from_spglib(positions, unit_cell, types)

    reciprocal_lattice = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_lattice, unit_cell)

    # Build q-points on a 2x2x2 grid
    nq_side = 2
    nq = nq_side^3
    q_points_cryst = zeros(Float64, 3, nq)
    iq = 0
    for ix in 0:nq_side-1
        for iy in 0:nq_side-1
            for iz in 0:nq_side-1
                iq += 1
                q_points_cryst[:, iq] = [ix / nq_side, iy / nq_side, iz / nq_side]
            end
        end
    end

    symmetry_group_q = SymmetriesQSpace(symmetry_group_uc, q_points_cryst)
    n_sym = length(symmetry_group_q)

    n_dims = 3
    n_atoms = 2
    n_modes = n_dims * n_atoms

    # Create a random complex vector in q-space (Cartesian)
    old_vector = randn(Complex{Float64}, n_modes, nq)

    # Compute rotate_vector! average over all symmetries
    avg_vector = zeros(Complex{Float64}, n_modes, nq)
    new_vector = zeros(Complex{Float64}, n_modes, nq)
    for i in 1:n_sym
        new_vector .= 0
        rotate_vector!(new_vector, old_vector, unit_cell, reciprocal_lattice, symmetry_group_q, i)
        avg_vector .+= new_vector
    end
    avg_vector ./= n_sym

    # Compare at gamma with symmetrize_vector_cartesian_q!
    # (it converts cart→cryst, calls symmetrize_vector_q!, converts back)
    ref_vector = copy(old_vector)
    symmetrize_vector_cartesian_q!(ref_vector, unit_cell, symmetry_group_q)

    if verbose
        println("n_sym = ", n_sym)
        println("rotate avg at gamma: ", real.(avg_vector[:, 1]))
        println("symmetrize at gamma: ", real.(ref_vector[:, 1]))
        println("diff at gamma: ", maximum(abs.(real.(avg_vector[:, 1]) - real.(ref_vector[:, 1]))))
    end

    # Compare real parts at gamma (symmetrize_vector_cartesian_q! returns real gamma)
    @test isapprox(real.(avg_vector[:, 1]), real.(ref_vector[:, 1]); atol=1e-10)
    # Check that result is non-trivial (not all zeros)
    @test maximum(abs.(real.(ref_vector[:, 1]))) > 1e-12
end


function test_rotate_dynamical_matrix_qspace(; verbose=false)
    # Use the same setup as test_symmetrize_cartesian_qspace
    q_tot = [0.0 -0.08970485720864457 -0.08970485720864457 0.0 -0.08970485720864457 -0.08970485720864457 -0.0 0.0; 0.0 -0.05179112345678227 0.05179112345678227 0.10358224691409373 0.05179112345678227 -0.05179112345678227 -0.10358224691409373 0.0; 0.0 0.07324370920331642 -0.07324370920331642 0.07324370920331642 0.03662185460165821 -0.03662185460165821 0.03662185460165821 -0.10986556380497464]

    unit_cell_structure = zeros(Float64, 3, 1)
    unit_cell = [2.94954603 1.4747730150000002 1.4747730150000002; 0.0 2.554381791611538 0.8514605972038461; 0.0 0.0 2.408294248783948]
    unit_cell .*= 1.889725989

    symmetry_group_uc = get_symmetry_group_from_spglib(unit_cell_structure, unit_cell, [1])

    reciprocal_lattice = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_lattice, unit_cell)

    q_points_cryst = zeros(Float64, size(q_tot)...)
    cryst_cart_conv!(q_points_cryst, q_tot, unit_cell, reciprocal_lattice, false; q_space=true)

    symmetry_group_q = SymmetriesQSpace(symmetry_group_uc, q_points_cryst)
    n_sym = length(symmetry_group_q)

    n_dims = 3
    n_atoms = 1
    n_modes = n_dims * n_atoms
    nq = size(q_tot, 2)

    # Create a random complex dynamical matrix in q-space
    old_matrix = randn(Complex{Float64}, n_modes, n_modes, nq)
    # Make Hermitian at each q
    for iq in 1:nq
        old_matrix[:, :, iq] .= (old_matrix[:, :, iq] .+ old_matrix[:, :, iq]') ./ 2
    end

    avg_matrix = zeros(Complex{Float64}, n_modes, n_modes, nq)
    new_matrix = zeros(Complex{Float64}, n_modes, n_modes, nq)
    for i in 1:n_sym
        new_matrix .= 0
        rotate_dynamical_matrix!(new_matrix, old_matrix, unit_cell, reciprocal_lattice, symmetry_group_q, i)
        avg_matrix .+= new_matrix
    end
    avg_matrix ./= n_sym

    # Compare with symmetrize_matrix_cartesian_q!
    ref_matrix = copy(old_matrix)
    symmetrize_matrix_cartesian_q!(ref_matrix, unit_cell, symmetry_group_q)

    if verbose
        println("diff: ", maximum(abs.(avg_matrix - ref_matrix)))
    end

    @test isapprox(avg_matrix, ref_matrix; atol=1e-10)
end


function test_rotate_centroid_real(; verbose=false)
    # Use PbTe unit cell (non-orthogonal, 2 atoms)
    a = 12.21
    cell = collect([-a 0.0 a; 0.0 a a; -a a 0.0]')
    crystal = collect([0.0 0.0 0.0; 0.5 0.5 0.5]')

    sym_group = get_symmetry_group_from_spglib(crystal, cell, [1, 2])
    n_dims = 3
    n_atoms = 2
    n_modes = n_dims * n_atoms
    n_sym = length(sym_group)

    reciprocal_vectors = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_vectors, cell)

    positions = zeros(Float64, n_dims, n_atoms)
    cryst_cart_conv!(positions, crystal, cell, reciprocal_vectors, true)

    # Create random atomic positions in Cartesian coordinates
    # Start with crystal coordinates, then convert to Cartesian
    old_centroid = randn(Float64, n_dims, n_atoms) * 0.01
    old_centroid .+= positions

    # Compute rotate_centroid! average over all symmetries
    avg_centroid = zeros(Float64, n_modes)
    new_centroid = zeros(Float64, n_modes)
    for i in 1:n_sym
        new_centroid .= 0
        rotate_centroid!(new_centroid, reshape(old_centroid, :), cell, reciprocal_vectors, sym_group, i)
        avg_centroid .+= new_centroid
    end
    avg_centroid ./= n_sym

    # Compare with symmetrize_positions!
    # symmetrize_positions! expects positions as matrix (n_dims × n_atoms)

    
    if verbose
        println("Number of symmetries: ", n_sym)
        println("rotate_centroid avg: ", avg_centroid)
        println("symmetrize_positions: ", reshape(positions, :))
        println("diff: ", maximum(abs.(avg_centroid - reshape(positions,:))))
    end

    @test isapprox(avg_centroid, reshape(positions, :); atol=1e-10)
end


function test_rotate_centroid_identity(; verbose=false)
    # Test with identity symmetry group (no translations)
    n_dims = 3
    n_atoms = 2
    n_modes = n_dims * n_atoms
    
    cell = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]
    sym_group = AtomicSymmetries.get_identity_symmetry_group(Float64; dims=n_dims, n_atoms=n_atoms, translations=false)
    
    reciprocal_vectors = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_vectors, cell)
    
    # Create random centroid
    old_centroid = randn(Float64, n_modes)
    
    # Apply rotate_centroid! with identity symmetry (index 1)
    new_centroid = zeros(Float64, n_modes)
    rotate_centroid!(new_centroid, old_centroid, cell, reciprocal_vectors, sym_group, 1)
    
    # With identity symmetry and no translations, result should be the same as input
    # (after coordinate conversion)
    if verbose
        println("Identity test - old: ", old_centroid)
        println("Identity test - new: ", new_centroid)
        println("diff: ", maximum(abs.(new_centroid - old_centroid)))
    end
    
    @test isapprox(new_centroid, old_centroid; atol=1e-10)
end


function test_rotate_centroid_supercell(; verbose=false)
    # PbTe 2x2x2 supercell (non-orthogonal, 16 atoms, 384 symmetries)
    a = 12.21
    cell_uc = collect([-a 0.0 a; 0.0 a a; -a a 0.0]')
    crystal_uc = collect([0.0 0.0 0.0; 0.5 0.5 0.5]')

    supercell = [2, 2, 2]
    nat_uc = size(crystal_uc, 2)
    nat_sc = prod(supercell) * nat_uc
    n_dims = 3
    n_modes = n_dims * nat_sc
    types_uc = [1, 2]

    # Build supercell positions in Cartesian
    uc_cart = cell_uc * crystal_uc
    cell = similar(cell_uc)
    for i in 1:3, j in 1:3
        cell[i, j] = cell_uc[i, j] * supercell[j]
    end

    positions_cart = zeros(Float64, 3, nat_sc)
    sc_types = zeros(Int, nat_sc)
    idx = 0
    for ix in 0:supercell[1]-1, iy in 0:supercell[2]-1, iz in 0:supercell[3]-1
        for k in 1:nat_uc
            idx += 1
            @views positions_cart[:, idx] = uc_cart[:, k] + ix*cell_uc[:, 1] + iy*cell_uc[:, 2] + iz*cell_uc[:, 3]
            sc_types[idx] = types_uc[k]
        end
    end

    # Convert to crystal coords and get symmetry group
    crystal = zeros(Float64, 3, nat_sc)
    get_crystal_coords!(crystal, positions_cart, cell)

    sym_group = get_symmetry_group_from_spglib(crystal, cell, sc_types)
    n_sym = length(sym_group)

    reciprocal_vectors = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_vectors, cell)

    # Equilibrium positions in Cartesian
    positions = zeros(Float64, n_dims, nat_sc)
    cryst_cart_conv!(positions, crystal, cell, reciprocal_vectors, true)

    # Perturb positions slightly around equilibrium
    old_centroid = copy(positions)
    old_centroid .+= randn(Float64, n_dims, nat_sc) * 0.01

    # Average rotate_centroid! over all symmetries
    avg_centroid = zeros(Float64, n_modes)
    new_centroid = zeros(Float64, n_modes)
    for i in 1:n_sym
        new_centroid .= 0
        rotate_centroid!(new_centroid, reshape(old_centroid, :), cell, reciprocal_vectors, sym_group, i)
        avg_centroid .+= new_centroid
    end
    avg_centroid ./= n_sym

    # The average should recover the Wyckoff (equilibrium) positions
    ref = reshape(positions, :)

    if verbose
        println("Number of symmetries: ", n_sym)
        println("diff: ", maximum(abs.(avg_centroid - ref)))
    end

    @test isapprox(avg_centroid, ref; atol=1e-10)
end


if abspath(PROGRAM_FILE) == @__FILE__
    test_rotate_vector_real()
    test_rotate_matrix_real()
    test_rotate_centroid_real()
    test_rotate_centroid_supercell()
    test_rotate_dynamical_matrix_qspace(; verbose=true)
end
