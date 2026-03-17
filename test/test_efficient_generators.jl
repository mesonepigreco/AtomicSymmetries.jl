using AtomicSymmetries
using Test
using LinearAlgebra

"""
Test that apply_symmetry_tensor! for rank-2 matches apply_sym_fc!,
and that rank-3 symmetrize is idempotent.
"""
function test_apply_symmetry_tensor()
    # PbTe unit cell
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'
    nat = 2
    atomic_numbers = [1, 2]

    symmetry_group = get_symmetry_group_from_spglib(positions, cell, atomic_numbers)
    dim = AtomicSymmetries.get_dimensions(symmetry_group)
    n_sym = get_nsymmetries(symmetry_group)
    n_modes = dim * nat

    # ── Rank-2: compare with apply_sym_fc! ──
    fc = randn(Float64, n_modes, n_modes)
    fc += fc'

    # Convert to crystal coords
    fc_cryst = similar(fc)
    AtomicSymmetries.cart_cryst_matrix_conversion!(fc_cryst, fc, cell; cart_to_cryst=true)

    # Method 1: apply_sym_fc! (existing code)
    result_fc = zeros(Float64, n_modes, n_modes)
    for s in 1:n_sym
        AtomicSymmetries.apply_sym_fc!(result_fc, fc_cryst,
            symmetry_group.symmetries[s], dim, symmetry_group.irt[s])
    end
    result_fc ./= n_sym

    # Method 2: apply_symmetry_tensor! (new code)
    result_tensor = zeros(Float64, n_modes, n_modes)
    for s in 1:n_sym
        apply_symmetry_tensor!(result_tensor, fc_cryst,
            symmetry_group.symmetries[s], dim, symmetry_group.irt[s])
    end
    result_tensor ./= n_sym

    @test result_fc ≈ result_tensor atol=1e-10

    # ── Rank-3: symmetrize idempotency ──
    tensor3 = randn(Float64, n_modes, n_modes, n_modes)
    # Make it permutation-symmetric
    for i in 1:n_modes, j in 1:n_modes, k in 1:n_modes
        v = (tensor3[i,j,k] + tensor3[i,k,j] + tensor3[j,i,k] +
             tensor3[j,k,i] + tensor3[k,i,j] + tensor3[k,j,i]) / 6
        tensor3[i,j,k] = v
        tensor3[i,k,j] = v
        tensor3[j,i,k] = v
        tensor3[j,k,i] = v
        tensor3[k,i,j] = v
        tensor3[k,j,i] = v
    end

    tensor3_sym = copy(tensor3)
    symmetrize_tensor!(tensor3_sym, cell, symmetry_group)

    # Symmetrize again: should be idempotent
    tensor3_sym2 = copy(tensor3_sym)
    symmetrize_tensor!(tensor3_sym2, cell, symmetry_group)

    @test tensor3_sym ≈ tensor3_sym2 atol=1e-10
end


"""
Test that rank-2 tensor generators match the existing get_matrix_generators.
"""
function test_rank2_consistency()
    # PbTe 2×2×2 supercell (same as test_pbte_generators.jl)
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'
    nat = 2

    # Generate supercell
    new_pos = zeros(Float64, 3, 16)
    for i in 1:2
        for j in 1:2
            for k in 1:2
                for iat in 1:nat
                    new_pos[:, (i-1)*8 + (j-1)*4 + (k-1)*2 + iat] = positions[:, iat] + [i-1, j-1, k-1]
                end
            end
        end
    end
    new_pos /= 2.0
    cell *= 2.0
    atomic_numbers = [(i - 1) % nat + 1 for i in 1:nat*8]

    symmetry_group = get_symmetry_group_from_spglib(new_pos, cell, atomic_numbers)

    # Old method: get_matrix_generators
    old_generators = AtomicSymmetries.get_matrix_generators(symmetry_group, cell)

    # New method: get_tensor_generators with rank=2
    new_generators = get_tensor_generators(symmetry_group, cell; rank=2)

    # Same number of independent generators
    @test length(old_generators) == length(new_generators)

    # Verify reconstruction: random fc → symmetrize → project → reconstruct
    n_modes = 3 * 16
    fc = randn(Float64, n_modes, n_modes)
    fc += fc'

    # Symmetrize with symmetry group
    fc_sym = copy(fc)
    symmetrize_fc!(fc_sym, cell, symmetry_group)

    # Project onto new generators and reconstruct
    coeffs = zeros(Float64, length(new_generators))
    get_coefficients_from_tensor!(coeffs, fc_sym, new_generators, cell)

    fc_reconstructed = zeros(Float64, n_modes, n_modes)
    reconstruct_tensor!(fc_reconstructed, new_generators, coeffs, cell)

    @test fc_sym ≈ fc_reconstructed atol=1e-8
end


"""
Test rank-3 generators: find generators, symmetrize → project → reconstruct round-trip.
Uses zinc-blende (F-43m, non-centrosymmetric) since centrosymmetric groups
force all fully-symmetric odd-rank tensors to zero.
"""
function test_rank3_generators()
    # Zinc-blende (F-43m, 2 atoms, non-centrosymmetric)
    positions, cell, types = get_zincblende()

    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)
    dim = AtomicSymmetries.get_dimensions(symmetry_group)
    nat = AtomicSymmetries.get_n_atoms(symmetry_group)
    n_modes = dim * nat

    # Find rank-3 generators
    generators = get_tensor_generators(symmetry_group, cell; rank=3)
    @test length(generators) > 0

    # Create a random permutation-symmetric rank-3 tensor and symmetrize it
    tensor = randn(Float64, n_modes, n_modes, n_modes)
    # Make permutation-symmetric
    tensor_sym = zeros(Float64, n_modes, n_modes, n_modes)
    for i in 1:n_modes, j in 1:n_modes, k in 1:n_modes
        v = (tensor[i,j,k] + tensor[i,k,j] + tensor[j,i,k] +
             tensor[j,k,i] + tensor[k,i,j] + tensor[k,j,i]) / 6
        tensor_sym[i,j,k] = v
    end

    symmetrize_tensor!(tensor_sym, cell, symmetry_group)

    # Project and reconstruct
    coeffs = zeros(Float64, length(generators))
    get_coefficients_from_tensor!(coeffs, tensor_sym, generators, cell)

    tensor_recon = zeros(Float64, n_modes, n_modes, n_modes)
    reconstruct_tensor!(tensor_recon, generators, coeffs, cell)

    @test tensor_sym ≈ tensor_recon atol=1e-8
end


"""
Test contraction: contract rank-2 generator with vector,
compare with explicit matrix-vector product.
"""
function test_contraction()
    # Pm-3m perovskite
    positions, cell, types = get_pm3m_perovskite()
    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)
    dim = AtomicSymmetries.get_dimensions(symmetry_group)
    nat = AtomicSymmetries.get_n_atoms(symmetry_group)
    n_modes = dim * nat

    # Get rank-2 generators
    generators = get_tensor_generators(symmetry_group, cell; rank=2)
    @test length(generators) > 0

    # Pick the first generator
    gen = generators[1]
    n_targets = size(gen.cartesian_blocks, 3)

    # Reconstruct the full matrix from this single generator
    full_matrix = zeros(Float64, n_modes, n_modes)
    coeff = [1.0]
    reconstruct_tensor!(full_matrix, [gen], coeff, cell)

    # Create a random Cartesian vector
    vector_cart = randn(Float64, n_modes)

    # Explicit matrix-vector product (in Cartesian)
    explicit_result = full_matrix * vector_cart

    # Contract using generator (contract on index 2: matrix times vector)
    result_blocks = zeros(Float64, dim, n_targets)
    result_target_atoms = zeros(Int, 1, n_targets)
    contract_generator_vector!(result_blocks, result_target_atoms, gen, vector_cart, 2)

    # Reconstruct contracted vector in Cartesian
    contracted_cart = zeros(Float64, n_modes)
    for t in 1:n_targets
        atom = result_target_atoms[1, t]
        contracted_cart[dim*(atom-1)+1:dim*atom] .+= gen.normalization .* result_blocks[:, t]
    end

    @test explicit_result ≈ contracted_cart atol=1e-10
end


"""
Test compact reconstruction: reconstruct full tensor from Generator struct,
verify matches symmetrize_tensor! applied to the seed.
"""
function test_compact_reconstruction()
    # Pm-3m perovskite
    positions, cell, types = get_pm3m_perovskite()
    symmetry_group = get_symmetry_group_from_spglib(positions, cell, types)
    dim = AtomicSymmetries.get_dimensions(symmetry_group)
    nat = AtomicSymmetries.get_n_atoms(symmetry_group)
    n_modes = dim * nat

    # Get rank-2 generators
    generators = get_tensor_generators(symmetry_group, cell; rank=2)

    # Pick the first generator, reconstruct it as a full matrix
    gen = generators[1]
    full_tensor = zeros(Float64, n_modes, n_modes)
    coeff = [1.0]
    reconstruct_tensor!(full_tensor, [gen], coeff, cell)

    # Symmetrizing the result should not change it (it's already symmetric)
    full_tensor_sym = copy(full_tensor)
    symmetrize_fc!(full_tensor_sym, cell, symmetry_group)

    @test full_tensor ≈ full_tensor_sym atol=1e-8

    # Check that the norm is approximately 1
    @test norm(full_tensor) ≈ 1.0 atol=1e-6
end


"""
Test index helpers: round-trip linearized <-> mode indices.
"""
function test_index_helpers()
    n_modes = 15  # e.g. 5 atoms * 3 dims
    rank = 3

    modes = zeros(Int, rank)
    for idx in 1:n_modes^rank
        AtomicSymmetries.linearized_to_mode_indices!(modes, idx, n_modes)
        idx2 = AtomicSymmetries.mode_indices_to_linear(modes, n_modes)
        @test idx == idx2
    end

    # Test atom/cartesian decomposition round-trip
    dim = 3
    nat = 5
    atom_idx = zeros(Int, rank)
    cart_idx = zeros(Int, rank)
    for idx in 1:n_modes^rank
        AtomicSymmetries.get_atomic_indices!(atom_idx, idx, dim, nat)
        AtomicSymmetries.get_cartesian_indices!(cart_idx, idx, dim, nat)

        # Reconstruct mode indices
        for j in 1:rank
            mode = AtomicSymmetries.atom_cart_to_mode(atom_idx[j], cart_idx[j], dim)
            modes_check = zeros(Int, rank)
            AtomicSymmetries.linearized_to_mode_indices!(modes_check, idx, n_modes)
            @test mode == modes_check[j]
        end
    end
end


"""
Test that get_tensor_generators_fast produces the same results as
get_tensor_generators on a non-orthogonal cell (PbTe unit cell).

Both methods must yield the same number of generators, and
projecting a random matrix onto the fast generators and reconstructing
must match direct symmetrization via symmetrize_fc!.
"""
function test_fast_vs_standard_generators()
    # PbTe unit cell (non-orthogonal cell, 2 atoms)
    a = 12.21
    cell = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions = [0.0 0.0 0.0; 0.5 0.5 0.5]'
    atomic_numbers = [1, 2]

    symmetry_group = get_symmetry_group_from_spglib(positions, cell, atomic_numbers)
    dim = AtomicSymmetries.get_dimensions(symmetry_group)
    nat = AtomicSymmetries.get_n_atoms(symmetry_group)
    n_modes = dim * nat

    # Standard method (known correct)
    generators_std = get_tensor_generators(symmetry_group, cell; rank=2)

    # Fast method (orbit decomposition)
    generators_fast = get_tensor_generators_fast(symmetry_group, cell; rank=2)

    # Same number of generators
    @test length(generators_std) == length(generators_fast)

    # Random symmetric matrix
    fc = randn(n_modes, n_modes)
    fc = (fc + fc') / 2

    # Direct symmetrization (ground truth)
    fc_sym = copy(fc)
    symmetrize_fc!(fc_sym, cell, symmetry_group)

    # Project onto fast generators and reconstruct
    coeffs_fast = zeros(length(generators_fast))
    get_coefficients_from_tensor!(coeffs_fast, fc, generators_fast, cell)
    fc_fast = zeros(n_modes, n_modes)
    reconstruct_tensor!(fc_fast, generators_fast, coeffs_fast, cell)

    # Must match direct symmetrization
    @test fc_fast ≈ fc_sym atol = 1e-8
end


"""
Test distance-based cutoff for tensor generator construction.

Uses the PbTe 2×2×2 supercell (16 atoms) from test_rank2_consistency.
"""
function test_cutoff_generators()
    # PbTe 2×2×2 supercell
    a = 12.21
    cell_uc = [-a 0.0 a; 0.0 a a; -a a 0.0]'
    positions_uc = [0.0 0.0 0.0; 0.5 0.5 0.5]'
    nat_uc = 2

    new_pos = zeros(Float64, 3, 16)
    for i in 1:2
        for j in 1:2
            for k in 1:2
                for iat in 1:nat_uc
                    new_pos[:, (i-1)*8 + (j-1)*4 + (k-1)*2 + iat] = positions_uc[:, iat] + [i-1, j-1, k-1]
                end
            end
        end
    end
    new_pos /= 2.0
    cell = cell_uc * 2.0
    atomic_numbers = [(i - 1) % nat_uc + 1 for i in 1:nat_uc*8]

    symmetry_group = get_symmetry_group_from_spglib(new_pos, cell, atomic_numbers)
    nat = AtomicSymmetries.get_n_atoms(symmetry_group)
    dim = AtomicSymmetries.get_dimensions(symmetry_group)
    n_modes = dim * nat

    # The positions are already in crystal coordinates (from the construction)
    positions_cryst = new_pos

    # --- Test 1: Inf cutoff should give the same generators as no cutoff ---
    generators_nocutoff = get_tensor_generators(symmetry_group, cell; rank=2)
    generators_inf = get_tensor_generators(symmetry_group, cell; rank=2,
        positions=positions_cryst, cutoff=Inf)
    @test length(generators_nocutoff) == length(generators_inf)

    # --- Test 2: Finite cutoff gives fewer generators ---
    # Use a cutoff that is shorter than the supercell diagonal but allows
    # nearest-neighbor interactions (PbTe nearest-neighbor distance ~ a/sqrt(2) ≈ 8.63 Å)
    cutoff_nn = a * 1.1  # slightly above nearest-neighbor distance
    generators_cutoff = get_tensor_generators(symmetry_group, cell; rank=2,
        positions=positions_cryst, cutoff=cutoff_nn)
    @test length(generators_cutoff) < length(generators_nocutoff)
    @test length(generators_cutoff) > 0

    # --- Test 3: Cutoff generators produce a valid symmetry-invariant tensor ---
    # Project a random symmetric tensor onto the cutoff generators, reconstruct,
    # then verify the result is symmetry-invariant
    fc = randn(Float64, n_modes, n_modes)
    fc = (fc + fc') / 2

    coeffs = zeros(Float64, length(generators_cutoff))
    get_coefficients_from_tensor!(coeffs, fc, generators_cutoff, cell)

    fc_recon = zeros(Float64, n_modes, n_modes)
    reconstruct_tensor!(fc_recon, generators_cutoff, coeffs, cell)

    # Symmetrize the reconstruction: should not change it
    fc_recon_sym = copy(fc_recon)
    symmetrize_fc!(fc_recon_sym, cell, symmetry_group)
    @test fc_recon ≈ fc_recon_sym atol=1e-8

    # --- Test 4: Both methods agree ---
    generators_fast_cutoff = get_tensor_generators_fast(symmetry_group, cell; rank=2,
        positions=positions_cryst, cutoff=cutoff_nn)
    @test length(generators_cutoff) == length(generators_fast_cutoff)

    # Also check fast method with Inf cutoff
    generators_fast_inf = get_tensor_generators_fast(symmetry_group, cell; rank=2,
        positions=positions_cryst, cutoff=Inf)
    generators_fast_nocutoff = get_tensor_generators_fast(symmetry_group, cell; rank=2)
    @test length(generators_fast_inf) == length(generators_fast_nocutoff)

    # --- Test 5: Error on missing positions ---
    @test_throws ArgumentError get_tensor_generators(symmetry_group, cell; rank=2,
        cutoff=5.0)
    @test_throws ArgumentError get_tensor_generators_fast(symmetry_group, cell; rank=2,
        cutoff=5.0)
end


if abspath(PROGRAM_FILE) == @__FILE__
    include("define_cell.jl")
    test_index_helpers()
    println("Index helpers passed")
    test_apply_symmetry_tensor()
    println("Apply symmetry tensor passed")
    test_rank2_consistency()
    println("Rank-2 consistency passed")
    test_rank3_generators()
    println("Rank-3 generators passed")
    test_contraction()
    println("Contraction passed")
    test_compact_reconstruction()
    println("Compact reconstruction passed")
    test_fast_vs_standard_generators()
    println("Fast vs standard generators passed")
    test_cutoff_generators()
    println("Cutoff generators passed")
    println("All efficient generator tests passed!")
end
