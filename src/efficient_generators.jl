@doc raw"""
    struct Generator{T, N, M}

An efficient memory representation of a symmetry generator for rank-N tensors.

Each generator stores only the small `dim^N` cartesian block per distinct target
atom tuple, never materializing the full `(dim*nat)^N` tensor.

The generator represents the symmetrization of a unit seed tensor:
starting from a tensor with 1 at position `(atom_indices, cartesian_indices)`
(plus permutation-symmetric copies), symmetrized over the full symmetry group.

# Fields
- `index::Int` — Linear seed index in `n_modes^N` space.
- `atom_indices::Vector{Int}` — Length N, sorted (canonical atom tuple of the seed).
- `cartesian_indices::Vector{Int}` — Length N, cartesian directions of the seed.
- `target_atoms::Matrix{Int}` — `(N, n_targets)`: distinct target atom tuples.
- `cartesian_blocks::Array{T, M}` — `(dim,...,dim, n_targets)`: accumulated blocks. M = N+1.
- `normalization::T` — Overall normalization factor.
- `dimension::Int` — Spatial dimension (typically 3).
- `n_atoms::Int` — Number of atoms.
"""
struct Generator{T,N,M}
    index::Int
    atom_indices::Vector{Int}
    cartesian_indices::Vector{Int}
    target_atoms::Matrix{Int}
    cartesian_blocks::Array{T,M}
    normalization::T
    dimension::Int
    n_atoms::Int
end

# ─────────────────────────────────────────────────────────────────────────────
# Index helpers
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    linearized_to_mode_indices!(modes, index, n_modes)

Convert a 1-based linear index to a k-tuple of 1-based mode indices.
Mode indices are in column-major order: the first index varies fastest.
"""
function linearized_to_mode_indices!(modes::AbstractVector{Int}, index::Int, n_modes::Int)
    rank = length(modes)
    idx = index - 1
    for i in 1:rank
        modes[i] = idx % n_modes + 1
        idx = idx ÷ n_modes
    end
end

@doc raw"""
    mode_indices_to_linear(modes, n_modes)

Convert a k-tuple of 1-based mode indices to a 1-based linear index.
"""
function mode_indices_to_linear(modes::AbstractVector{Int}, n_modes::Int)
    rank = length(modes)
    idx = 0
    for i in rank:-1:1
        idx = idx * n_modes + (modes[i] - 1)
    end
    return idx + 1
end

@doc raw"""
    get_atomic_indices!(indices, index, dimension, nat)

Convert a 1-based linear index in `n_modes^rank` space to the atom indices.
Mode `m` (0-based) decomposes as `m = cart + dim * atom`, so
`atom = m ÷ dim` and `cart = m % dim`.
"""
function get_atomic_indices!(indices::AbstractVector{Int}, index::Int, dimension::Int, nat::Int)
    rank = length(indices)
    n_modes = dimension * nat
    idx = index - 1
    for i in 1:rank
        mode = idx % n_modes
        indices[i] = mode ÷ dimension + 1
        idx = idx ÷ n_modes
    end
end

@doc raw"""
    get_cartesian_indices!(indices, index, dimension, nat)

Convert a 1-based linear index in `n_modes^rank` space to the cartesian indices.
"""
function get_cartesian_indices!(indices::AbstractVector{Int}, index::Int, dimension::Int, nat::Int)
    rank = length(indices)
    n_modes = dimension * nat
    idx = index - 1
    for i in 1:rank
        mode = idx % n_modes
        indices[i] = mode % dimension + 1
        idx = idx ÷ n_modes
    end
end

@doc raw"""
    atom_cart_to_mode(atom, cart, dim)

Convert (1-based atom index, 1-based cartesian index) to 1-based mode index.
"""
function atom_cart_to_mode(atom::Int, cart::Int, dim::Int)
    return (atom - 1) * dim + cart
end

# ─────────────────────────────────────────────────────────────────────────────
# Permutation helpers
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    sorted_permutations(k)

Return all k! permutations of (1,...,k) as a vector of vectors.
"""
function sorted_permutations(k::Int)
    if k == 0
        return [Int[]]
    end
    if k == 1
        return [[1]]
    end
    result = Vector{Vector{Int}}()
    _permute!(result, collect(1:k), 1)
    return result
end

function _permute!(result, arr, start)
    if start > length(arr)
        push!(result, copy(arr))
        return
    end
    for i in start:length(arr)
        arr[start], arr[i] = arr[i], arr[start]
        _permute!(result, arr, start + 1)
        arr[start], arr[i] = arr[i], arr[start]
    end
end

@doc raw"""
    is_sorted_tuple(t)

Return true if the vector t is sorted in non-decreasing order.
"""
function is_sorted_tuple(t::AbstractVector{Int})
    for i in 1:length(t)-1
        if t[i] > t[i+1]
            return false
        end
    end
    return true
end

@doc raw"""
    is_canonical_under_irt(atom_tuple, symmetry_group)

Return true if no symmetry operation maps this sorted atom tuple to a
lexicographically smaller sorted tuple (i.e., this is the canonical representative
of its orbit under the symmetry group's atom permutations).
"""
function is_canonical_under_irt(atom_tuple::AbstractVector{Int}, symmetry_group::Symmetries)
    n_sym = get_nsymmetries(symmetry_group)
    n_atoms = get_n_atoms(symmetry_group)
    rank = length(atom_tuple)

    mapped = zeros(Int, rank)
    irt_inv = zeros(Int, n_atoms)

    for s in 1:n_sym
        # Build irt_inv for this symmetry
        for k in 1:n_atoms
            irt_inv[symmetry_group.irt[s][k]] = k
        end

        # Map atoms through inverse permutation and sort
        for j in 1:rank
            mapped[j] = irt_inv[atom_tuple[j]]
        end
        sort!(mapped)

        # Compare lexicographically with atom_tuple
        for j in 1:rank
            if mapped[j] < atom_tuple[j]
                return false
            elseif mapped[j] > atom_tuple[j]
                break
            end
        end
    end
    return true
end

@doc raw"""
    get_little_group(atom_tuple, symmetry_group; is_permutational=true)

Return the indices of the symmetries that belong to the Little Group (stabilizer)
of the given atom tuple.

A symmetry $s$ is in the Little Group if the transformed tuple $irt_s(T)$
is a permutation of $T$ (if `is_permutational=true`) or exactly equal to $T$
(if `is_permutational=false`).
"""
function get_little_group(atom_tuple::AbstractVector{Int}, symmetry_group::Symmetries; is_permutational=true)
    n_sym = get_nsymmetries(symmetry_group)
    rank = length(atom_tuple)
    little_group_indices = Int[]

    # Pre-sort for fast comparison if permutational
    target = is_permutational ? sort(atom_tuple) : atom_tuple
    mapped = zeros(Int, rank)

    for s in 1:n_sym
        irt = symmetry_group.irt[s]
        for j in 1:rank
            mapped[j] = irt[atom_tuple[j]]
        end

        if is_permutational
            sort!(mapped)
        end

        if mapped == target
            push!(little_group_indices, s)
        end
    end

    return little_group_indices
end

@doc raw"""
    get_atom_orbit_representatives(sym_group)

Return the indices of atoms that are representatives of their orbits under the symmetry group.
"""
function get_atom_orbit_representatives(sym_group::Symmetries)
    nat = get_n_atoms(sym_group)
    nsym = get_nsymmetries(sym_group)
    reps = Int[]
    seen = fill(false, nat)
    for i in 1:nat
        if !seen[i]
            push!(reps, i)
            for s in 1:nsym
                seen[sym_group.irt[s][i]] = true
            end
        end
    end
    return reps
end

@doc raw"""
    _kron_power(S, k)

Compute the k-fold Kronecker product S ⊗ S ⊗ ... ⊗ S.
For dim=3, rank=4 this produces an 81×81 matrix.
"""
function _kron_power(S::AbstractMatrix{T}, k::Int) where T
    result = S
    for _ in 2:k
        result = kron(result, S)
    end
    return result
end

@doc raw"""
    _build_perm_map(sigma, dim, rank)

Build a linear-index permutation map for a Cartesian block permutation.
`sigma` specifies that `result[c₁,...,cₖ] = block[c_{σ(1)},...,c_{σ(k)}]`.

Returns a vector `perm_map` such that `result_flat[i] = block_flat[perm_map[i]]`.
"""
function _build_perm_map(sigma::AbstractVector{Int}, dim::Int, rank::Int)
    block_shape = ntuple(_ -> dim, rank)
    lin = LinearIndices(block_shape)
    ci_all = CartesianIndices(block_shape)
    perm_map = Vector{Int}(undef, dim^rank)
    for ci in ci_all
        t = Tuple(ci)
        permuted = ntuple(j -> t[sigma[j]], rank)
        perm_map[lin[ci]] = lin[permuted...]
    end
    return perm_map
end

@doc raw"""
    _build_sym_lin_indices(sym_indices, dim, rank)

Convert NTuple-based sym_indices to linear index vectors for fast flat-array access.
Returns a vector of vectors of linear indices.
"""
function _build_sym_lin_indices(sym_indices::Vector{Vector{NTuple{N,Int}}}, dim::Int) where N
    block_shape = ntuple(_ -> dim, N)
    lin = LinearIndices(block_shape)
    sym_lin = Vector{Vector{Int}}(undef, length(sym_indices))
    for i in eachindex(sym_indices)
        sym_lin[i] = [lin[idx...] for idx in sym_indices[i]]
    end
    return sym_lin
end

@doc raw"""
    _apply_rotation_to_block!(result, block, S)

In-place application of rotation $S$ to all indices of the Cartesian block.
Convention matches `_contract_rotation_dim!` (2-arg) and `apply_sym_fc!`:
$result_{c_1 \dots c_k} = \sum_{d_1 \dots d_k} S_{d_1 c_1} \dots S_{d_k c_k} block_{d_1 \dots d_k}$.

For rank 2 this gives $\text{result} = S^T \cdot \text{block} \cdot S$.
"""
function _apply_rotation_to_block!(result::AbstractArray{T,rank}, block::AbstractArray{T,rank},
    S::AbstractMatrix{T}, tmp::AbstractArray{T,rank}, current::AbstractArray{T,rank}) where {T,rank}
    dim = size(S, 1)
    # Use pre-allocated workspace buffers instead of allocating new ones
    tmp .= block
    current .= block

    for d in 1:rank
        # Contract dimension d with S
        # result[..., β, ...] = sum_α S[α, β] * current[..., α, ...]
        fill!(tmp, zero(T))
        _contract_rotation_dim!(tmp, current, S, d, dim, rank)
        current .= tmp
    end
    result .= current
end

@doc raw"""
    _contract_rotation_dim!(dest, src, S, which_dim, dim, rank)

Specialized contraction for a single Cartesian dimension (3-argument form).
Convention: `dest[..., β, ...] = Σ_α S[α, β] * src[..., α, ...]`

This matches the 2-argument `_contract_rotation_dim!` used by `apply_symmetry_tensor!`.
"""
function _contract_rotation_dim!(dest, src, S, which_dim, dim, rank)
    cart_idx = ones(Int, rank)
    src_idx = zeros(Int, rank)
    total = dim^rank
    for _ in 1:total
        val = zero(eltype(dest))
        for alpha in 1:dim
            # Build src_idx manually: replace which_dim with alpha
            for j in 1:rank
                src_idx[j] = j == which_dim ? alpha : cart_idx[j]
            end
            val += S[alpha, cart_idx[which_dim]] * src[src_idx...]
        end
        dest[cart_idx...] = val

        # Increment multi-index
        for j in 1:rank
            cart_idx[j] += 1
            if cart_idx[j] <= dim
                break
            end
            cart_idx[j] = 1
        end
    end
end

@doc raw"""
    _permute_block_indices!(result, block, sigma)

In-place permutation of the Cartesian indices of the block according to sigma.
$result_{c_1 \dots c_k} = block_{c_{\sigma(1)} \dots c_{\sigma(k)}}$.
"""
function _permute_block_indices!(result::AbstractArray{T,rank}, block::AbstractArray{T,rank},
    sigma::AbstractVector{Int}) where {T,rank}
    cart_idx = ones(Int, rank)
    permuted_idx = zeros(Int, rank)
    total = size(block, 1)^rank
    for _ in 1:total
        # Permute index according to sigma
        for j in 1:rank
            permuted_idx[j] = cart_idx[sigma[j]]
        end
        result[cart_idx...] = block[permuted_idx...]

        # Increment multi-index
        for j in 1:rank
            cart_idx[j] += 1
            if cart_idx[j] <= size(block, 1)
                break
            end
            cart_idx[j] = 1
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Generator constructor
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    Generator{T, rank}(index, symmetry_group, cell; buffer=default_buffer())

Construct a Generator of the given rank from a linear seed index and symmetry group.

The generator is the symmetrization of a unit seed tensor. For rank-2, this
matches `get_matrix_generator!`; for rank-1, `get_vector_generator!`.

The constructor:
1. Decomposes the index into atom/cartesian indices
2. Builds the full seed tensor (with permutation symmetry)
3. Symmetrizes it: cart→cryst, average over all symmetries, cryst→cart
4. Extracts the compact block representation from the symmetrized tensor
"""
function Generator{T,rank}(index::Int, symmetry_group::Symmetries,
    cell::AbstractMatrix; buffer=default_buffer()) where {T,rank}
    dim = symmetry_group.dimension
    n_sym = get_nsymmetries(symmetry_group)
    nat = get_n_atoms(symmetry_group)
    n_modes = dim * nat

    # Decompose index
    atom_indices = zeros(Int, rank)
    cartesian_indices = zeros(Int, rank)
    get_atomic_indices!(atom_indices, index, dim, nat)
    get_cartesian_indices!(cartesian_indices, index, dim, nat)

    # Build seed tensor: set 1 at the seed position (+ permutation-symmetric copies)
    tensor = zeros(T, ntuple(_ -> n_modes, rank))
    mode_indices = [atom_cart_to_mode(atom_indices[j], cartesian_indices[j], dim) for j in 1:rank]

    # Set all permutations of the mode indices to 1
    perms = sorted_permutations(rank)
    perm_modes = zeros(Int, rank)
    for perm in perms
        for j in 1:rank
            perm_modes[j] = mode_indices[perm[j]]
        end
        tensor[perm_modes...] = one(T)
    end

    # Symmetrize the tensor
    symmetrize_tensor!(tensor, cell, symmetry_group)

    # Extract compact representation: find all nonzero atom-tuple blocks
    # Collect unique atom tuples with nonzero blocks
    target_list = Vector{Vector{Int}}()
    block_list = Vector{Array{T}}()

    atom_idx = ones(Int, rank)
    total_tuples = nat^rank
    block_dims = ntuple(_ -> dim, rank)
    thr = 1e-15

    for _ in 1:total_tuples
        # Extract block
        ranges = ntuple(j -> (dim*(atom_idx[j]-1)+1):(dim*atom_idx[j]), rank)
        block = tensor[ranges...]

        if maximum(abs, block) > thr
            push!(target_list, copy(atom_idx))
            push!(block_list, copy(block))
        end

        # Increment atom multi-index
        for j in 1:rank
            atom_idx[j] += 1
            if atom_idx[j] <= nat
                break
            end
            atom_idx[j] = 1
        end
    end

    n_targets = length(target_list)
    if n_targets == 0
        # Zero generator
        target_atoms = zeros(Int, rank, 0)
        cartesian_blocks = zeros(T, block_dims..., 0)
        return Generator{T,rank,rank + 1}(index, atom_indices, cartesian_indices,
            target_atoms, cartesian_blocks,
            zero(T), dim, nat)
    end

    target_atoms = zeros(Int, rank, n_targets)
    cartesian_blocks = zeros(T, block_dims..., n_targets)
    for i in 1:n_targets
        target_atoms[:, i] .= target_list[i]
        selectdim(cartesian_blocks, rank + 1, i) .= block_list[i]
    end

    # Normalization: 1 (already normalized by symmetrize_tensor! which divides by n_sym)
    normalization = one(T)

    return Generator{T,rank,rank + 1}(index, atom_indices, cartesian_indices,
        target_atoms, cartesian_blocks,
        normalization, dim, nat)
end

# ─────────────────────────────────────────────────────────────────────────────
# Tensor operations (operate on full tensors for verification)
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    apply_symmetry_tensor!(result, tensor, S, dim, irt; buffer=default_buffer())

Apply one symmetry operation to a rank-N tensor in crystal coordinates.

Generalizes `apply_sym_fc!` (rank 2) and `apply_sym_centroid!` (rank 1):
loops over atom N-tuples, extracts dim^k block at source atoms,
contracts S on each index, and writes to destination.

The result is **added** to `result` (not overwritten), to match the
convention of `apply_sym_fc!`.
"""
function apply_symmetry_tensor!(result::AbstractArray{T}, tensor::AbstractArray{T},
    S::AbstractMatrix, dim::Int, irt::AbstractVector{Int};
    buffer=default_buffer()) where T
    rank = ndims(tensor)
    nat = size(tensor, 1) ÷ dim

    # For each source atom tuple, compute the rotated block and add to result
    atom_idx = ones(Int, rank)
    target_idx = zeros(Int, rank)
    total_tuples = nat^rank

    block_dims = ntuple(_ -> dim, rank)

    for _ in 1:total_tuples
        # Get target atoms via irt
        for j in 1:rank
            target_idx[j] = irt[atom_idx[j]]
        end

        # Extract source block
        src_ranges = ntuple(j -> (dim*(atom_idx[j]-1)+1):(dim*atom_idx[j]), rank)
        src_block = copy(tensor[src_ranges...])

        # Apply rotation on each index: S' * block * S (generalized)
        # For consistency with apply_sym_fc!: result_{irt[a],irt[b]}[β₁,β₂] += Σ S[α₁,β₁] * tensor[a,b][α₁,α₂] * S[α₂,β₂]
        # Which is: S' * tensor_block * S for rank 2
        # Generalized: contract S^T on each dimension
        tmp_block = copy(src_block)
        for j in 1:rank
            _contract_rotation_dim!(tmp_block, S, j, dim, rank)
        end

        # Add to result at target position
        tgt_ranges = ntuple(j -> (dim*(target_idx[j]-1)+1):(dim*target_idx[j]), rank)
        view(result, tgt_ranges...) .+= tmp_block

        # Increment atom multi-index
        for j in 1:rank
            atom_idx[j] += 1
            if atom_idx[j] <= nat
                break
            end
            atom_idx[j] = 1
        end
    end
end

@doc raw"""
    _contract_rotation_dim!(block, S, which_dim, dim, rank)

Contract the rotation matrix S on dimension `which_dim` of the block in-place.
Convention: `new_block[..., β, ...] = Σ_α S[α, β] * old_block[..., α, ...]`

This matches apply_sym_fc! where: result = S' * source * S,
i.e. result[β₁,β₂] = Σ_{α₁,α₂} S[α₁,β₁] * source[α₁,α₂] * S[α₂,β₂]
"""
function _contract_rotation_dim!(block::AbstractArray{T}, S, which_dim, dim, rank) where T
    block_copy = copy(block)

    cart_idx = ones(Int, rank)
    total = dim^rank
    for _ in 1:total
        val = zero(T)
        for alpha in 1:dim
            src_idx = ntuple(j -> j == which_dim ? alpha : cart_idx[j], rank)
            val += S[alpha, cart_idx[which_dim]] * block_copy[src_idx...]
        end
        block[cart_idx...] = val

        # Increment
        for j in 1:rank
            cart_idx[j] += 1
            if cart_idx[j] <= dim
                break
            end
            cart_idx[j] = 1
        end
    end
end

@doc raw"""
    tensor_cryst_cart_conversion!(dest, src, cell, dim, n_atoms; cart_to_cryst=true, buffer=default_buffer())

Rank-N coordinate conversion generalizing `cart_cryst_matrix_conversion!`.

For each atom N-tuple, the dim^N block is transformed by contracting
the coordinate transformation matrix on each index.
"""
function tensor_cryst_cart_conversion!(dest::AbstractArray{T}, src::AbstractArray{T},
    cell::AbstractMatrix, dim::Int, n_atoms::Int;
    cart_to_cryst::Bool=true, buffer=default_buffer()) where T
    rank = ndims(dest)

    metric_tensor = cell' * cell
    inv_metric_tensor = inv(metric_tensor)
    transform_matrix = inv_metric_tensor * cell'

    if cart_to_cryst
        transform_matrix = inv(transform_matrix)
    end

    dest .= zero(T)

    atom_idx = ones(Int, rank)
    total_tuples = n_atoms^rank

    for _ in 1:total_tuples
        ranges = ntuple(j -> (dim*(atom_idx[j]-1)+1):(dim*atom_idx[j]), rank)
        block = copy(src[ranges...])

        # Contract transform_matrix on each dimension
        for j in 1:rank
            _contract_transform_dim!(block, transform_matrix, j, dim, rank)
        end

        dest[ranges...] .= block

        # Increment
        for j in 1:rank
            atom_idx[j] += 1
            if atom_idx[j] <= n_atoms
                break
            end
            atom_idx[j] = 1
        end
    end
end

@doc raw"""
    _contract_transform_dim!(block, M, which_dim, dim, rank)

Contract transform matrix M on dimension `which_dim`:
`new[...,β,...] = Σ_α M'[β,α] * old[...,α,...]`
"""
function _contract_transform_dim!(block::AbstractArray{T}, M, which_dim, dim, rank) where T
    block_copy = copy(block)
    cart_idx = ones(Int, rank)
    src_idx = zeros(Int, rank)
    total = dim^rank
    for _ in 1:total
        val = zero(T)
        for alpha in 1:dim
            # Build src_idx manually: replace which_dim with alpha
            for j in 1:rank
                src_idx[j] = j == which_dim ? alpha : cart_idx[j]
            end
            val += M[alpha, cart_idx[which_dim]] * block_copy[src_idx...]
        end
        block[cart_idx...] = val

        for j in 1:rank
            cart_idx[j] += 1
            if cart_idx[j] <= dim
                break
            end
            cart_idx[j] = 1
        end
    end
end

@doc raw"""
    symmetrize_tensor!(tensor, cell, symmetry_group; buffer=default_buffer())

Symmetrize a rank-N tensor in Cartesian coordinates.
Converts cart→cryst, averages over all symmetries, converts cryst→cart.

This generalizes `symmetrize_fc!` (rank 2) and `symmetrize_vector!` (rank 1).
"""
function symmetrize_tensor!(tensor::AbstractArray{T}, cell::AbstractMatrix,
    symmetry_group::Symmetries; buffer=default_buffer()) where T
    rank = ndims(tensor)

    # Delegate to existing optimized implementations for rank ≤ 2
    if rank == 1
        symmetrize_vector!(tensor, cell, symmetry_group; buffer=buffer)
        return
    elseif rank == 2
        symmetrize_fc!(tensor, cell, symmetry_group; buffer=buffer)
        return
    end

    # Generic implementation for rank ≥ 3
    dim = symmetry_group.dimension
    nat = get_n_atoms(symmetry_group)
    n_sym = get_nsymmetries(symmetry_group)

    tensor_cryst = similar(tensor)
    tensor_cryst_cart_conversion!(tensor_cryst, tensor, cell, dim, nat; cart_to_cryst=true)

    result = zeros(T, size(tensor))
    for s in 1:n_sym
        apply_symmetry_tensor!(result, tensor_cryst, symmetry_group.symmetries[s],
            dim, symmetry_group.irt[s])
    end
    result ./= n_sym

    tensor_cryst_cart_conversion!(tensor, result, cell, dim, nat; cart_to_cryst=false)
end

# ─────────────────────────────────────────────────────────────────────────────
# Fast Generator Finding via Orbit Decomposition
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    get_tensor_generators_fast(symmetry_group, cell; rank, type=Float64)

Faster version of `get_tensor_generators` using orbit decomposition.
"""
function get_tensor_generators_fast(symmetry_group::Symmetries{U}, cell::AbstractMatrix{T};
    rank::Int, type::Type=Float64) where {U,T}

    dim = get_dimensions(symmetry_group)
    nat = get_n_atoms(symmetry_group)
    n_sym = get_nsymmetries(symmetry_group)

    generators = Vector{Generator{type,rank,rank + 1}}()

    # Find orbit representatives of atoms
    atom_reps = get_atom_orbit_representatives(symmetry_group)

    # Temporary storage for canonicality check
    atom_tuple = ones(Int, rank)

    # Iterate over sorted atom tuples starting with a primitive representative
    for i1 in atom_reps
        atom_tuple[1] = i1
        # Recursively fill the rest of the tuple: i1 <= i2 <= ... <= ik
        _iterate_sorted_from!(generators, atom_tuple, 2, symmetry_group, cell, rank, type)
    end

    return generators
end

function _iterate_sorted_from!(generators, tuple, pos, sym_group, cell, rank, type)
    nat = get_n_atoms(sym_group)
    if pos > rank
        # Tuple completely built. Check if canonical.
        if is_canonical_under_irt(tuple, sym_group)
            # Find independent generators for this orbit
            _find_generators_for_orbit!(generators, tuple, sym_group, cell, rank, type)
        end
        return
    end

    for i in tuple[pos-1]:nat
        tuple[pos] = i
        _iterate_sorted_from!(generators, tuple, pos + 1, sym_group, cell, rank, type)
    end
end

function _find_generators_for_orbit!(generators, tuple, sym_group, cell, rank, T)
    # 1. Little Group
    H = get_little_group(tuple, sym_group; is_permutational=true)

    # 2. Invariant basis of Cartesian blocks
    dim = sym_group.dimension
    nat = get_n_atoms(sym_group)
    basis = _get_invariant_cartesian_basis(H, tuple, sym_group, rank, T)

    if isempty(basis)
        return
    end

    block_size = dim^rank
    block_shape = ntuple(_ -> dim, rank)
    n_sym = get_nsymmetries(sym_group)

    # Precompute crystal-to-Cartesian transform matrix and its Kronecker power
    metric_tensor = cell' * cell
    inv_metric_tensor = inv(metric_tensor)
    cryst_to_cart_M = inv_metric_tensor * cell'
    # Transpose: convention is new[β] = Σ_α M[α,β] * old[α], so need M' in Kronecker
    kron_M = Matrix(_kron_power(cryst_to_cart_M, rank)')

    # Precompute Kronecker rotation matrices for all symmetries (transposed for S[α,β] convention)
    kron_S_all = Vector{Matrix{T}}(undef, n_sym)
    for s in 1:n_sym
        kron_S_all[s] = Matrix(_kron_power(sym_group.symmetries[s], rank)')
    end

    # Precompute permutation maps and identity permutation for expansion step
    all_perms = sorted_permutations(rank)
    identity_perm = collect(1:rank)
    perm_maps_expand = Dict{Vector{Int}, Vector{Int}}()
    for perm in all_perms
        if perm != identity_perm
            perm_maps_expand[perm] = _build_perm_map(perm, dim, rank)
        end
    end

    # Workspace flat vectors for orbit rotation + coordinate conversion
    rot_flat = Vector{T}(undef, block_size)
    perm_flat = Vector{T}(undef, block_size)
    conv_flat = Vector{T}(undef, block_size)
    mapped = Vector{Int}(undef, rank)
    target_buf = Vector{Int}(undef, rank)

    # 3. For each basis vector, build a Generator
    for (i_gen, B_rep) in enumerate(basis)
        B_rep_flat = vec(B_rep)

        # Find all tuples in the orbit and their blocks (in crystal coordinates)
        orbit_dict = Dict{Vector{Int},Array{T,rank}}()

        for s in 1:n_sym
            irt = sym_group.irt[s]

            # Target tuple (unsorted then sorted)
            for j in 1:rank
                mapped[j] = irt[tuple[j]]
                target_buf[j] = mapped[j]
            end
            sort!(target_buf)
            target = copy(target_buf)

            if !haskey(orbit_dict, target)
                sigma_fixed = _find_correct_permutation(target, mapped)
                pm = _build_perm_map(sigma_fixed, dim, rank)

                # Rotate: rot = kron(S,...,S) * B_rep
                mul!(rot_flat, kron_S_all[s], B_rep_flat)

                # Permute: perm[i] = rot[pm[i]]
                @inbounds for i in 1:block_size
                    perm_flat[i] = rot_flat[pm[i]]
                end

                orbit_dict[target] = reshape(copy(perm_flat), block_shape)
            end
        end

        # Convert all blocks from crystal to Cartesian coordinates using kron_M
        for (_, block) in orbit_dict
            block_flat = vec(block)
            mul!(conv_flat, kron_M, block_flat)
            copyto!(block_flat, conv_flat)
        end

        # Expand sorted tuples to include all permutations (permutation symmetry):
        # the Generator stores blocks at ALL atom tuples, not just sorted ones.
        # For unsorted tuple π(T), the block has its Cartesian indices permuted by π.
        expanded = Dict{Vector{Int},Array{T,rank}}()
        for (sorted_target, block) in orbit_dict
            block_flat = vec(block)
            for perm in all_perms
                perm_target = [sorted_target[perm[j]] for j in 1:rank]
                if !haskey(expanded, perm_target)
                    if perm == identity_perm
                        expanded[perm_target] = block
                    else
                        pm = perm_maps_expand[perm]
                        perm_block_flat = Vector{T}(undef, block_size)
                        @inbounds for i in 1:block_size
                            perm_block_flat[i] = block_flat[pm[i]]
                        end
                        expanded[perm_target] = reshape(perm_block_flat, block_shape)
                    end
                end
            end
        end

        # Convert dict to Generator
        n_modes = dim * nat
        mode_indices = [atom_cart_to_mode(tuple[j], 1, dim) for j in 1:rank]
        lin_idx = mode_indices_to_linear(mode_indices, n_modes)

        gen = _blocks_to_generator(T, expanded, lin_idx, tuple, ones(Int, rank), dim, nat)

        # Normalize
        gnorm = _generator_frobenius_norm(gen)
        if gnorm > 1e-10
            push!(generators, _make_normalized_generator(gen, gnorm))
        end
    end
end

function _find_correct_permutation(sorted_T, mapped_T)
    rank = length(sorted_T)
    sigma = zeros(Int, rank)
    used = fill(false, rank)
    for j in 1:rank
        val = mapped_T[j]
        for k in 1:rank
            if !used[k] && sorted_T[k] == val
                sigma[j] = k
                used[k] = true
                break
            end
        end
    end
    return sigma
end

function _get_invariant_cartesian_basis(H, tuple, sym_group, rank, T)
    dim = sym_group.dimension
    sym_indices = _get_tuple_symmetric_indices(tuple, dim)
    D_sym = length(sym_indices)
    block_size = dim^rank

    P = zeros(T, D_sym, D_sym)

    # Precompute linear indices for each symmetry group (for fast flat-array access)
    sym_lin = _build_sym_lin_indices(sym_indices, dim)

    # Precompute inverse sqrt of group sizes
    inv_sqrt_len = [one(T) / sqrt(T(length(sym_indices[j]))) for j in 1:D_sym]

    # Precompute Kronecker rotation matrices and permutation maps for each symmetry in H
    mapped = Vector{Int}(undef, rank)
    kron_perm_list = Vector{Tuple{Matrix{T}, Vector{Int}}}(undef, length(H))
    for (k, s_idx) in enumerate(H)
        S = sym_group.symmetries[s_idx]
        irt = sym_group.irt[s_idx]
        for j in 1:rank
            mapped[j] = irt[tuple[j]]
        end
        sigma = _find_correct_permutation(tuple, mapped)
        # Transpose: convention is dest[β] = Σ_α S[α,β] * src[α], so need S' in Kronecker
        kron_perm_list[k] = (Matrix(_kron_power(S, rank)'), _build_perm_map(sigma, dim, rank))
    end

    # Workspace: 3 flat vectors of length block_size, reused across iterations
    B_work_flat = zeros(T, block_size)
    B_rot_flat = zeros(T, block_size)
    B_perm_flat = zeros(T, block_size)

    # Hot loop: accumulate projection matrix P
    for (kS, pm) in kron_perm_list
        for j in 1:D_sym
            # Fill B_work_flat: set entries corresponding to sym_indices[j] to normalized value
            fill!(B_work_flat, zero(T))
            isqrt = inv_sqrt_len[j]
            @inbounds for lin in sym_lin[j]
                B_work_flat[lin] = isqrt
            end

            # BLAS rotation: B_rot = kron(S,...,S) * B_work
            mul!(B_rot_flat, kS, B_work_flat)

            # Gather permutation: B_perm[i] = B_rot[perm_map[i]]
            @inbounds for i in 1:block_size
                B_perm_flat[i] = B_rot_flat[pm[i]]
            end

            # Accumulate into P matrix
            for i in 1:D_sym
                val = zero(T)
                @inbounds for lin in sym_lin[i]
                    val += B_perm_flat[lin]
                end
                P[i, j] += val * inv_sqrt_len[i]
            end
        end
    end

    P ./= length(H)

    if all(abs.(P) .< 1e-12)
        return []
    end

    U, Svals, V = svd(P)
    basis_vectors = Vector{Array{T,rank}}()
    block_shape = ntuple(_ -> dim, rank)
    for i in 1:D_sym
        if Svals[i] > 0.5
            B_flat = zeros(T, block_size)
            for j in 1:D_sym
                coeff = U[j, i] * inv_sqrt_len[j]
                @inbounds for lin in sym_lin[j]
                    B_flat[lin] += coeff
                end
            end
            push!(basis_vectors, reshape(B_flat, block_shape))
        end
    end

    return basis_vectors
end

function _get_symmetric_indices(rank, dim)
    indices = Vector{Vector{NTuple{rank,Int}}}()
    seen = Set{NTuple{rank,Int}}()
    curr = ones(Int, rank)
    while true
        t = ntuple(i -> curr[i], rank)
        if !(t in seen)
            perms = _get_all_permutations(t)
            push!(indices, perms)
            for p in perms
                push!(seen, p)
            end
        end
        done = true
        for j in 1:rank
            curr[j] += 1
            if curr[j] <= dim
                done = false
                break
            end
            curr[j] = 1
        end
        if done
            break
        end
    end
    return indices
end

@doc raw"""
    _get_tuple_symmetric_indices(atom_tuple, dim)

Return a basis for Cartesian blocks that are symmetric under internal permutations
of the atom tuple (i.e., permutations that map the tuple to itself).
"""
function _get_tuple_symmetric_indices(atom_tuple::AbstractVector{Int}, dim::Int)
    rank = length(atom_tuple)

    # Identify groups of indices that share the same atom
    groups = Vector{Vector{Int}}()
    seen_atoms = Int[]
    for (idx, a) in enumerate(atom_tuple)
        if a in seen_atoms
            p = findfirst(==(a), seen_atoms)
            push!(groups[p], idx)
        else
            push!(seen_atoms, a)
            push!(groups, [idx])
        end
    end

    indices = Vector{Vector{NTuple{rank,Int}}}()
    seen = Set{NTuple{rank,Int}}()
    curr = ones(Int, rank)
    while true
        t = ntuple(i -> curr[i], rank)
        if !(t in seen)
            # Find all permutations that only swap indices WITHIN atom groups
            perms = _get_tuple_permutations(t, groups)
            push!(indices, perms)
            for p in perms
                push!(seen, p)
            end
        end
        done = true
        for j in 1:rank
            curr[j] += 1
            if curr[j] <= dim
                done = false
                break
            end
            curr[j] = 1
        end
        if done
            break
        end
    end
    return indices
end

function _get_tuple_permutations(t::NTuple{N,Int}, groups::Vector{Vector{Int}}) where N
    res = Set{NTuple{N,Int}}()
    push!(res, t)

    # Work on a copy to allow permutations
    t_v = collect(t)

    for g in groups
        if length(g) > 1
            new_res = Set{NTuple{N,Int}}()
            for existing_t in res
                sub_t = collect(existing_t[g])
                sub_perms = Set{NTuple{length(g),Int}}()
                _permute_recursive!(sub_perms, sub_t, 1, length(g))

                for p in sub_perms
                    candidate = collect(existing_t)
                    candidate[g] .= p
                    push!(new_res, ntuple(i -> candidate[i], N))
                end
            end
            union!(res, new_res)
        end
    end
    return collect(res)
end

function _get_all_permutations(t::NTuple{N,T}) where {N,T}
    res = Set{NTuple{N,T}}()
    _permute_recursive!(res, collect(t), 1, N)
    return collect(res)
end

function _permute_recursive!(res, a, k, n)
    if k == n
        push!(res, ntuple(i -> a[i], n))
    else
        for i in k:n
            a[k], a[i] = a[i], a[k]
            _permute_recursive!(res, a, k + 1, n)
            a[k], a[i] = a[i], a[k]
        end
    end
end
# Block-level arithmetic for compact Gram-Schmidt
# ─────────────────────────────────────────────────────────────────────────────

"""
    _blocks_from_generator(gen) -> Dict{Vector{Int}, Array{T, rank}}

Extract a mutable dictionary of (atom_tuple => Cartesian block) from a Generator.
Blocks are scaled by gen.normalization so the dict directly represents the tensor values.
"""
function _blocks_from_generator(gen::Generator{T,N,M}) where {T,N,M}
    rank = N
    n_targets = size(gen.target_atoms, 2)
    blocks = Dict{Vector{Int},Array{T,N}}()
    for t in 1:n_targets
        key = gen.target_atoms[:, t]
        block = Array{T}(undef, ntuple(_ -> gen.dimension, rank))
        block .= selectdim(gen.cartesian_blocks, rank + 1, t) .* gen.normalization
        blocks[key] = block
    end
    return blocks
end

"""
    _blocks_dot(blocks1, blocks2) -> T

Compute the Frobenius inner product between two block-dictionaries.
Only matching atom tuples contribute.
"""
function _blocks_dot(blocks1::Dict{Vector{Int},Array{T,N}},
    blocks2::Dict{Vector{Int},Array{T,N}}) where {T,N}
    val = zero(T)
    for (key, b1) in blocks1
        b2 = get(blocks2, key, nothing)
        if b2 !== nothing
            val += sum(b1 .* b2)
        end
    end
    return val
end

"""
    _blocks_axpy!(blocks, α, other_blocks)

Compute blocks += α * other_blocks in-place.
Creates new entries for atom tuples present only in other_blocks.
"""
function _blocks_axpy!(blocks::Dict{Vector{Int},Array{T,N}},
    α::T, other_blocks::Dict{Vector{Int},Array{T,N}}) where {T,N}
    for (key, b_other) in other_blocks
        b = get(blocks, key, nothing)
        if b !== nothing
            b .+= α .* b_other
        else
            blocks[key] = α .* b_other
        end
    end
end

"""
    _blocks_norm(blocks) -> T

Compute the Frobenius norm of a block-dictionary.
"""
function _blocks_norm(blocks::Dict{Vector{Int},Array{T,N}}) where {T,N}
    val = zero(T)
    for (_, b) in blocks
        val += sum(abs2, b)
    end
    return sqrt(val)
end

"""
    _blocks_to_generator(blocks, index, atom_indices, cart_indices, dim, nat, rank) -> Generator

Convert a block-dictionary to a compact Generator struct, removing near-zero blocks.
"""
function _blocks_to_generator(::Type{T}, blocks::Dict{Vector{Int},Array{T,N}},
    index::Int, atom_indices::Vector{Int}, cart_indices::Vector{Int},
    dim::Int, nat::Int) where {T,N}
    rank = N
    thr = 1e-15

    target_list = Vector{Vector{Int}}()
    block_list = Vector{Array{T,N}}()

    for (key, block) in blocks
        if maximum(abs, block) > thr
            push!(target_list, copy(key))
            push!(block_list, copy(block))
        end
    end

    block_dims = ntuple(_ -> dim, rank)
    n_targets = length(target_list)
    if n_targets == 0
        target_atoms = zeros(Int, rank, 0)
        cartesian_blocks = zeros(T, block_dims..., 0)
    else
        target_atoms = zeros(Int, rank, n_targets)
        cartesian_blocks = zeros(T, block_dims..., n_targets)
        for i in 1:n_targets
            target_atoms[:, i] .= target_list[i]
            selectdim(cartesian_blocks, rank + 1, i) .= block_list[i]
        end
    end

    return Generator{T,rank,rank + 1}(index, copy(atom_indices), copy(cart_indices),
        target_atoms, cartesian_blocks,
        one(T), dim, nat)
end

# ─────────────────────────────────────────────────────────────────────────────
# get_tensor_generators — Find Independent Generators
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    get_tensor_generators(symmetry_group, cell; rank, type=Float64)

Find a set of independent generators for the invariant subspace of rank-`rank`
tensors under the given symmetry group, assuming full permutation symmetry
between tensor indices.

Returns a `Vector{Generator}` — the generators are orthonormal compact structs.

# Algorithm
1. Iterate over sorted atom k-tuples (only sorted due to permutation symmetry)
2. For each cartesian seed:
   - Build the symmetrized tensor from the unit seed (compact representation)
   - Check if zero → skip
   - Apply Modified Gram-Schmidt against existing generators using compact
     block-level arithmetic (never materializing the full tensor)
   - If the residual has significant norm → normalize and accept as new generator
"""
function get_tensor_generators(symmetry_group::Symmetries{U}, cell::AbstractMatrix{T};
    rank::Int, type::Type=Float64) where {U,T}

    dim = get_dimensions(symmetry_group)
    nat = get_n_atoms(symmetry_group)
    n_modes = dim * nat

    generators = Vector{Generator{type,rank,rank + 1}}()
    # Cache block-dicts of accepted generators for fast Gram-Schmidt projections
    gen_blocks_cache = Vector{Dict{Vector{Int},Array{type,rank}}}()

    # Iterate over all sorted atom tuples
    atom_tuple = ones(Int, rank)
    n_cart = dim^rank

    while true
        # Fast skip: if this atom tuple is not the canonical representative
        # of its orbit under the symmetry group, all its generators are
        # linear combinations of generators at the canonical representative.
        if !is_canonical_under_irt(atom_tuple, symmetry_group)
            if !_increment_sorted_tuple!(atom_tuple, nat)
                break
            end
            continue
        end

        cart_idx = ones(Int, rank)
        for _ in 1:n_cart
            # For permutation symmetry with equal atoms: only use sorted cart indices
            if _should_skip_cartesian(atom_tuple, cart_idx, rank)
                _increment_multi_index!(cart_idx, dim, rank)
                continue
            end

            # Build seed tensor and symmetrize (expensive)
            mode_indices = [atom_cart_to_mode(atom_tuple[j], cart_idx[j], dim) for j in 1:rank]
            lin_idx = mode_indices_to_linear(mode_indices, n_modes)

            gen = Generator{type,rank}(lin_idx, symmetry_group, cell)

            # Check if zero
            gen_norm = _generator_frobenius_norm(gen)
            if gen_norm < 1e-10
                _increment_multi_index!(cart_idx, dim, rank)
                continue
            end

            # Extract blocks and normalize
            candidate = _blocks_from_generator(gen)
            inv_norm = one(type) / gen_norm
            for (_, block) in candidate
                block .*= inv_norm
            end

            # Modified Gram-Schmidt: subtract projections onto existing generators
            for q_blocks in gen_blocks_cache
                overlap = _blocks_dot(candidate, q_blocks)
                if abs(overlap) > 1e-14
                    _blocks_axpy!(candidate, -overlap, q_blocks)
                end
            end

            # Check residual norm
            residual_norm = _blocks_norm(candidate)
            if residual_norm < 1e-8
                _increment_multi_index!(cart_idx, dim, rank)
                continue
            end

            # Normalize residual
            inv_res = one(type) / residual_norm
            for (_, block) in candidate
                block .*= inv_res
            end

            # Convert to Generator and accept
            atom_indices = zeros(Int, rank)
            cart_indices_vec = zeros(Int, rank)
            get_atomic_indices!(atom_indices, lin_idx, dim, nat)
            get_cartesian_indices!(cart_indices_vec, lin_idx, dim, nat)

            new_gen = _blocks_to_generator(type, candidate, lin_idx,
                atom_indices, cart_indices_vec, dim, nat)
            push!(generators, new_gen)
            push!(gen_blocks_cache, candidate)

            _increment_multi_index!(cart_idx, dim, rank)
        end

        # Increment sorted atom tuple
        if !_increment_sorted_tuple!(atom_tuple, nat)
            break
        end
    end

    return generators
end

@doc raw"""
    _extract_generator(type, rank, index, atom_indices, cart_indices, tensor, dim, nat)

Extract a compact Generator struct from a full normalized tensor.
Finds all nonzero atom-tuple blocks and stores them.
"""
function _extract_generator(::Type{T}, rank, index, atom_indices, cart_indices,
    tensor::AbstractArray{T}, dim, nat) where T
    target_list = Vector{Vector{Int}}()
    block_list = Vector{Array{T}}()
    block_dims = ntuple(_ -> dim, rank)
    thr = 1e-15

    atom_idx = ones(Int, rank)
    total_tuples = nat^rank
    for _ in 1:total_tuples
        ranges = ntuple(j -> (dim*(atom_idx[j]-1)+1):(dim*atom_idx[j]), rank)
        block = tensor[ranges...]
        if maximum(abs, block) > thr
            push!(target_list, copy(atom_idx))
            push!(block_list, copy(block))
        end

        for j in 1:rank
            atom_idx[j] += 1
            if atom_idx[j] <= nat
                break
            end
            atom_idx[j] = 1
        end
    end

    n_targets = length(target_list)
    if n_targets == 0
        target_atoms = zeros(Int, rank, 0)
        cartesian_blocks = zeros(T, block_dims..., 0)
    else
        target_atoms = zeros(Int, rank, n_targets)
        cartesian_blocks = zeros(T, block_dims..., n_targets)
        for i in 1:n_targets
            target_atoms[:, i] .= target_list[i]
            selectdim(cartesian_blocks, rank + 1, i) .= block_list[i]
        end
    end

    return Generator{T,rank,rank + 1}(index, copy(atom_indices), copy(cart_indices),
        target_atoms, cartesian_blocks,
        one(T), dim, nat)
end

@doc raw"""
    _should_skip_cartesian(atom_tuple, cart_idx, rank)

For permutation symmetry: when two atom indices are equal,
the cartesian indices for those positions should be sorted to avoid
generating the same generator twice.
"""
function _should_skip_cartesian(atom_tuple, cart_idx, rank)
    for i in 1:rank-1
        if atom_tuple[i] == atom_tuple[i+1] && cart_idx[i] > cart_idx[i+1]
            return true
        end
    end
    return false
end

@doc raw"""
    _generator_value_at_seed(gen, atom_tuple, cart_idx)

Check if a generator has a nonzero block value at the given atom tuple
and cartesian indices.
"""
function _generator_value_at_seed(gen::Generator{T,N,M}, atom_tuple, cart_idx) where {T,N,M}
    rank = N
    n_targets = size(gen.target_atoms, 2)
    for t in 1:n_targets
        match = true
        for j in 1:rank
            if gen.target_atoms[j, t] != atom_tuple[j]
                match = false
                break
            end
        end
        if match
            block = selectdim(gen.cartesian_blocks, rank + 1, t)
            return block[cart_idx...]
        end
    end
    return zero(T)
end

@doc raw"""
    _generator_frobenius_norm(gen)

Compute the Frobenius norm of the full tensor represented by a generator.
"""
function _generator_frobenius_norm(gen::Generator{T,N,M}) where {T,N,M}
    rank = N
    n_targets = size(gen.cartesian_blocks, rank + 1)
    total = zero(T)
    for t in 1:n_targets
        block = selectdim(gen.cartesian_blocks, rank + 1, t)
        total += sum(abs2, block)
    end
    return sqrt(total) * abs(gen.normalization)
end

@doc raw"""
    _make_normalized_generator(gen, current_norm)

Create a new generator with blocks scaled so the full tensor has unit Frobenius norm.
"""
function _make_normalized_generator(gen::Generator{T,N,M}, current_norm::T) where {T,N,M}
    scale = one(T) / current_norm
    new_blocks = gen.cartesian_blocks .* (scale * gen.normalization)
    return Generator{T,N,M}(gen.index, gen.atom_indices, gen.cartesian_indices,
        gen.target_atoms, new_blocks,
        one(T), gen.dimension, gen.n_atoms)
end

@doc raw"""
    _is_independent_generator(gen, existing_generators)

Check if a generator is independent from existing generators.
Uses dot products of the compact block representations.
"""
function _is_independent_generator(gen::Generator{T,N,M},
    existing::Vector{Generator{T,N,M}}) where {T,N,M}
    for other in existing
        overlap = _generator_dot(gen, other)
        if abs(overlap) > 1e-6
            return false
        end
    end
    return true
end

@doc raw"""
    _generator_dot(gen1, gen2)

Compute the dot product of two generators using their compact representations.

The full tensor element at atom tuple (a₁,...,aₖ) and cartesian indices (c₁,...,cₖ)
is Σ_t [target_atoms[:,t] == (a₁,...,aₖ)] * normalization * block[c₁,...,cₖ,t].

The dot product is:
Σ_{a₁,...,aₖ} Σ_{c₁,...,cₖ} gen1_value * gen2_value
= Σ_{t1,t2: target1[:,t1] == target2[:,t2]} norm1*norm2 * dot(block1[:,...,:,t1], block2[:,...,:,t2])
"""
function _generator_dot(gen1::Generator{T,N,M},
    gen2::Generator{T,N,M}) where {T,N,M}
    rank = N
    n_t1 = size(gen1.cartesian_blocks, rank + 1)
    n_t2 = size(gen2.cartesian_blocks, rank + 1)

    dot_val = zero(T)

    for t1 in 1:n_t1
        target1 = view(gen1.target_atoms, :, t1)
        for t2 in 1:n_t2
            target2 = view(gen2.target_atoms, :, t2)

            if target1 == target2
                block1 = selectdim(gen1.cartesian_blocks, rank + 1, t1)
                block2 = selectdim(gen2.cartesian_blocks, rank + 1, t2)
                dot_val += sum(block1 .* block2)
            end
        end
    end

    return dot_val * gen1.normalization * gen2.normalization
end

@doc raw"""
    _increment_sorted_tuple!(tuple, max_val)

Increment a sorted tuple to the next sorted tuple in lexicographic order.
Returns false if the tuple has reached its maximum (all values = max_val).
"""
function _increment_sorted_tuple!(tuple::AbstractVector{Int}, max_val::Int)
    k = length(tuple)
    # Find the rightmost position that can be incremented
    for i in k:-1:1
        if tuple[i] < max_val
            tuple[i] += 1
            # Reset all positions to the right to the minimum allowed value
            for j in i+1:k
                tuple[j] = tuple[i]
            end
            return true
        end
    end
    return false
end

function _increment_multi_index!(idx::AbstractVector{Int}, max_val::Int, rank::Int)
    for j in 1:rank
        idx[j] += 1
        if idx[j] <= max_val
            return
        end
        idx[j] = 1
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Contraction functions
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    contract_generator_vector!(result_blocks, result_target_atoms, gen, vector_cart, contract_index)

Contract a generator with a Cartesian vector along `contract_index`, reducing rank by 1.

Returns contracted blocks and reduced target atom tuples.
The result is in Cartesian coordinates (same as the generator blocks and input vector).
"""
function contract_generator_vector!(result_blocks::AbstractArray{T},
    result_target_atoms::AbstractMatrix{Int},
    gen::Generator{T,N,M},
    vector_cart::AbstractVector{T},
    contract_index::Int) where {T,N,M}

    rank = N
    dim = gen.dimension
    n_targets = size(gen.target_atoms, 2)

    for t in 1:n_targets
        target = view(gen.target_atoms, :, t)
        atom_c = target[contract_index]

        # Extract the vector block for the contracted atom (Cartesian)
        v = view(vector_cart, dim*(atom_c-1)+1:dim*atom_c)

        # Fill reduced target atoms
        idx = 0
        for j in 1:rank
            if j != contract_index
                idx += 1
                result_target_atoms[idx, t] = target[j]
            end
        end

        # Contract the block
        block = selectdim(gen.cartesian_blocks, rank + 1, t)
        new_rank = rank - 1
        if new_rank == 0
            result_view = selectdim(result_blocks, 1, t)
        else
            result_view = selectdim(result_blocks, new_rank + 1, t)
        end
        _contract_block_vector!(result_view, block, v, contract_index, dim, rank)
    end
end

function _contract_block_vector!(result, block, v, contract_dim, dim, rank)
    new_rank = rank - 1
    if new_rank == 0
        val = zero(eltype(result))
        for alpha in 1:dim
            idx = ntuple(j -> j == contract_dim ? alpha : 1, rank)
            val += block[idx...] * v[alpha]
        end
        result[] = val
        return
    end

    cart_idx = ones(Int, new_rank)
    total = dim^new_rank
    for _ in 1:total
        val = zero(eltype(result))
        for alpha in 1:dim
            full_idx = ntuple(rank) do j
                if j < contract_dim
                    cart_idx[j]
                elseif j == contract_dim
                    alpha
                else
                    cart_idx[j-1]
                end
            end
            val += block[full_idx...] * v[alpha]
        end
        result[cart_idx...] = val

        for j in 1:new_rank
            cart_idx[j] += 1
            if cart_idx[j] <= dim
                break
            end
            cart_idx[j] = 1
        end
    end
end

@doc raw"""
    contract_generator_matrix!(result_blocks, result_target_atoms, gen, matrix_cart, idx1, idx2)

Contract a generator with a Cartesian matrix along indices `idx1` and `idx2`,
reducing rank by 2.
"""
function contract_generator_matrix!(result_blocks::AbstractArray{T},
    result_target_atoms::AbstractMatrix{Int},
    gen::Generator{T,N,M},
    matrix_cart::AbstractMatrix{T},
    idx1::Int, idx2::Int) where {T,N,M}

    rank = N
    dim = gen.dimension
    n_targets = size(gen.target_atoms, 2)
    new_rank = rank - 2

    @assert idx1 != idx2
    if idx1 > idx2
        idx1, idx2 = idx2, idx1
    end

    for t in 1:n_targets
        target = view(gen.target_atoms, :, t)
        atom1 = target[idx1]
        atom2 = target[idx2]

        mat_block = view(matrix_cart, dim*(atom1-1)+1:dim*atom1, dim*(atom2-1)+1:dim*atom2)

        idx = 0
        for j in 1:rank
            if j != idx1 && j != idx2
                idx += 1
                result_target_atoms[idx, t] = target[j]
            end
        end

        block = selectdim(gen.cartesian_blocks, rank + 1, t)
        if new_rank == 0
            result_view = selectdim(result_blocks, 1, t)
            _contract_block_matrix_scalar!(result_view, block, mat_block, idx1, idx2, dim, rank)
        else
            result_view = selectdim(result_blocks, new_rank + 1, t)
            _contract_block_matrix!(result_view, block, mat_block, idx1, idx2, dim, rank)
        end
    end
end

function _contract_block_matrix_scalar!(result, block, mat, idx1, idx2, dim, rank)
    val = zero(eltype(result))
    for a1 in 1:dim
        for a2 in 1:dim
            full_idx = ntuple(rank) do j
                if j == idx1
                    a1
                elseif j == idx2
                    a2
                else
                    1
                end
            end
            val += block[full_idx...] * mat[a1, a2]
        end
    end
    result[] = val
end

function _contract_block_matrix!(result, block, mat, idx1, idx2, dim, rank)
    new_rank = rank - 2
    cart_idx = ones(Int, new_rank)
    total = dim^new_rank

    for _ in 1:total
        val = zero(eltype(result))
        for a1 in 1:dim
            for a2 in 1:dim
                full_idx = _build_contracted_index_2(cart_idx, a1, a2, idx1, idx2, rank)
                val += block[full_idx...] * mat[a1, a2]
            end
        end
        result[cart_idx...] = val

        for j in 1:new_rank
            cart_idx[j] += 1
            if cart_idx[j] <= dim
                break
            end
            cart_idx[j] = 1
        end
    end
end

function _build_contracted_index_2(cart_idx, a1, a2, idx1, idx2, rank)
    reduced_pos = 0
    ntuple(rank) do j
        if j == idx1
            a1
        elseif j == idx2
            a2
        else
            reduced_pos += 1
            cart_idx[reduced_pos]
        end
    end
end

@doc raw"""
    contract_generators(gen1, gen2)

Compute the inner product of two generators.
"""
function contract_generators(gen1::Generator{T,N,M},
    gen2::Generator{T,N,M}) where {T,N,M}
    return _generator_dot(gen1, gen2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Reconstruction & projection
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    reconstruct_tensor!(tensor, generators, coefficients, cell)

Reconstruct a full rank-N tensor in Cartesian coordinates from generators and coefficients.

For each generator g with coefficient α:
  scatter-add `α * g.normalization * g.cartesian_blocks[:,...,:, t]`
  at the target atom positions for each target t.

Since the generators already store blocks in Cartesian coordinates
(they were computed from `symmetrize_tensor!` which returns Cartesian),
no coordinate conversion is needed.
"""
function reconstruct_tensor!(tensor::AbstractArray{T},
    generators::Vector{Generator{T,N,M}},
    coefficients::AbstractVector{T},
    cell::AbstractMatrix) where {T,N,M}

    rank = N
    @assert length(generators) == length(coefficients)
    tensor .= zero(T)

    if isempty(generators)
        return
    end

    for (i, gen) in enumerate(generators)
        coeff = coefficients[i]
        n_targets = size(gen.target_atoms, 2)

        for t in 1:n_targets
            target = view(gen.target_atoms, :, t)
            block = selectdim(gen.cartesian_blocks, rank + 1, t)

            dim = gen.dimension
            ranges = ntuple(j -> (dim*(target[j]-1)+1):(dim*target[j]), rank)
            view(tensor, ranges...) .+= coeff .* gen.normalization .* block
        end
    end
end

@doc raw"""
    get_coefficients_from_tensor!(coeffs, tensor_cart, generators, cell)

Project an input Cartesian tensor onto the generator basis.

Since generators store blocks in Cartesian coordinates, we compute
the dot product directly in Cartesian space.
"""
function get_coefficients_from_tensor!(coeffs::AbstractVector{T},
    tensor_cart::AbstractArray{T},
    generators::Vector{Generator{T,N,M}},
    cell::AbstractMatrix) where {T,N,M}

    rank = N
    @assert length(coeffs) == length(generators)
    if isempty(generators)
        return
    end

    for (i, gen) in enumerate(generators)
        n_targets = size(gen.target_atoms, 2)
        val = zero(T)
        dim = gen.dimension

        for t in 1:n_targets
            target = view(gen.target_atoms, :, t)
            block = selectdim(gen.cartesian_blocks, rank + 1, t)

            ranges = ntuple(j -> (dim*(target[j]-1)+1):(dim*target[j]), rank)
            tensor_block = view(tensor_cart, ranges...)
            val += sum(block .* tensor_block)
        end

        coeffs[i] = val * gen.normalization
    end
end
