# ─────────────────────────────────────────────────────────────────────────────
# Ensemble average: compute generator coefficients from stochastic data
# ─────────────────────────────────────────────────────────────────────────────
#
# Following Bianco et al. (arXiv:1703.03212), high-rank force constant tensors
# can be estimated as ensemble averages:
#   Rank 2: Φ_ab     = -⟨v_a f_b⟩
#   Rank 3: Φ_abc    = -⟨v_a v_b f_c⟩
#   Rank 4: Φ_abcd   = -⟨v_a v_b v_c f_d⟩
#
# We compute generator coefficients directly from ensemble data, never forming
# the full (dim*nat)^N tensor.

# ─────────────────────────────────────────────────────────────────────────────
# Rank-specialized contraction kernels (private)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _contract_block_rank2(block, v_at, f_at, dim)

Contract a `(dim, dim)` block with displacement vector `v_at` and force vector `f_at`:
`Σ_{a,b} block[a,b] * v[a] * f[b]`
"""
@inline function _contract_block_rank2(block::AbstractMatrix{T},
    v_at::AbstractVector{T}, f_at::AbstractVector{T}, dim::Int) where T
    val = zero(T)
    @inbounds for b in 1:dim
        fb = f_at[b]
        @simd for a in 1:dim
            val += block[a, b] * v_at[a] * fb
        end
    end
    return val
end

"""
    _contract_block_rank3(block, v_at1, v_at2, f_at, dim)

Contract a `(dim, dim, dim)` block with two displacement vectors and one force vector:
`Σ_{a,b,c} block[a,b,c] * v1[a] * v2[b] * f[c]`
"""
@inline function _contract_block_rank3(block::AbstractArray{T,3},
    v_at1::AbstractVector{T}, v_at2::AbstractVector{T},
    f_at::AbstractVector{T}, dim::Int) where T
    val = zero(T)
    @inbounds for c in 1:dim
        fc = f_at[c]
        for b in 1:dim
            vb_fc = v_at2[b] * fc
            @simd for a in 1:dim
                val += block[a, b, c] * v_at1[a] * vb_fc
            end
        end
    end
    return val
end

"""
    _contract_block_rank4(block, v_at1, v_at2, v_at3, f_at, dim)

Contract a `(dim, dim, dim, dim)` block with three displacement vectors and one force vector:
`Σ_{a,b,c,d} block[a,b,c,d] * v1[a] * v2[b] * v3[c] * f[d]`
"""
@inline function _contract_block_rank4(block::AbstractArray{T,4},
    v_at1::AbstractVector{T}, v_at2::AbstractVector{T},
    v_at3::AbstractVector{T}, f_at::AbstractVector{T}, dim::Int) where T
    val = zero(T)
    @inbounds for d in 1:dim
        fd = f_at[d]
        for c in 1:dim
            vc_fd = v_at3[c] * fd
            for b in 1:dim
                vb_vc_fd = v_at2[b] * vc_fd
                @simd for a in 1:dim
                    val += block[a, b, c, d] * v_at1[a] * vb_vc_fd
                end
            end
        end
    end
    return val
end

# ─────────────────────────────────────────────────────────────────────────────
# Single-config accumulation (streaming API)
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    accumulate_ensemble_config!(coeffs, v_config, f_config, generators)

Accumulate the contribution of a single configuration to ensemble-averaged
generator coefficients. This is the streaming variant: the caller is responsible
for zeroing `coeffs` before the first call, and dividing by `-n_configs` after
all configurations have been processed.

# Arguments
- `coeffs::AbstractVector{T}` — length `n_generators`, accumulated in-place
- `v_config::AbstractMatrix{T}` — shape `(dim, nat)`, Cartesian displacements for one config
- `f_config::AbstractMatrix{T}` — shape `(dim, nat)`, Cartesian forces for one config
- `generators::Vector{Generator{T,N,M}}` — compact generators (rank encoded in N)

# Notes
Due to permutation symmetry of the generators, the last tensor index is always
contracted with forces. For rank-2: `v * f`, rank-3: `v * v * f`, rank-4: `v * v * v * f`.
"""
function accumulate_ensemble_config!(coeffs::AbstractVector{T},
    v_config::AbstractMatrix{T}, f_config::AbstractMatrix{T},
    generators::Vector{Generator{T,2,3}}) where T

    _accumulate_rank2!(coeffs, v_config, f_config, generators)
end

function accumulate_ensemble_config!(coeffs::AbstractVector{T},
    v_config::AbstractMatrix{T}, f_config::AbstractMatrix{T},
    generators::Vector{Generator{T,3,4}}) where T

    _accumulate_rank3!(coeffs, v_config, f_config, generators)
end

function accumulate_ensemble_config!(coeffs::AbstractVector{T},
    v_config::AbstractMatrix{T}, f_config::AbstractMatrix{T},
    generators::Vector{Generator{T,4,5}}) where T

    _accumulate_rank4!(coeffs, v_config, f_config, generators)
end

# ── Rank-2 kernel ──

function _accumulate_rank2!(coeffs::AbstractVector{T},
    v::AbstractMatrix{T}, f::AbstractMatrix{T},
    generators::Vector{Generator{T,2,3}}) where T

    @inbounds for (i_gen, gen) in enumerate(generators)
        dim = gen.dimension
        n_targets = size(gen.target_atoms, 2)
        gen_val = zero(T)

        for t in 1:n_targets
            a1 = gen.target_atoms[1, t]
            a2 = gen.target_atoms[2, t]
            block = @view gen.cartesian_blocks[:, :, t]
            v_at = @view v[:, a1]
            f_at = @view f[:, a2]
            gen_val += _contract_block_rank2(block, v_at, f_at, dim)
        end

        coeffs[i_gen] += gen_val * gen.normalization
    end
end

# ── Rank-3 kernel ──

function _accumulate_rank3!(coeffs::AbstractVector{T},
    v::AbstractMatrix{T}, f::AbstractMatrix{T},
    generators::Vector{Generator{T,3,4}}) where T

    @inbounds for (i_gen, gen) in enumerate(generators)
        dim = gen.dimension
        n_targets = size(gen.target_atoms, 2)
        gen_val = zero(T)

        for t in 1:n_targets
            a1 = gen.target_atoms[1, t]
            a2 = gen.target_atoms[2, t]
            a3 = gen.target_atoms[3, t]
            block = @view gen.cartesian_blocks[:, :, :, t]
            v_at1 = @view v[:, a1]
            v_at2 = @view v[:, a2]
            f_at  = @view f[:, a3]
            gen_val += _contract_block_rank3(block, v_at1, v_at2, f_at, dim)
        end

        coeffs[i_gen] += gen_val * gen.normalization
    end
end

# ── Rank-4 kernel ──

function _accumulate_rank4!(coeffs::AbstractVector{T},
    v::AbstractMatrix{T}, f::AbstractMatrix{T},
    generators::Vector{Generator{T,4,5}}) where T

    @inbounds for (i_gen, gen) in enumerate(generators)
        dim = gen.dimension
        n_targets = size(gen.target_atoms, 2)
        gen_val = zero(T)

        for t in 1:n_targets
            a1 = gen.target_atoms[1, t]
            a2 = gen.target_atoms[2, t]
            a3 = gen.target_atoms[3, t]
            a4 = gen.target_atoms[4, t]
            block = @view gen.cartesian_blocks[:, :, :, :, t]
            v_at1 = @view v[:, a1]
            v_at2 = @view v[:, a2]
            v_at3 = @view v[:, a3]
            f_at  = @view f[:, a4]
            gen_val += _contract_block_rank4(block, v_at1, v_at2, v_at3, f_at, dim)
        end

        coeffs[i_gen] += gen_val * gen.normalization
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Batch API: process all configurations at once
# ─────────────────────────────────────────────────────────────────────────────

@doc raw"""
    get_coefficients_from_ensemble!(coeffs, v, f, generators)

Compute generator coefficients from ensemble-averaged displacement-force
correlations, following Bianco et al. (arXiv:1703.03212).

For rank-N force constants:
```math
\Phi_{a_1 \cdots a_N} = -\frac{1}{N_c} \sum_I v_{a_1}^I \cdots v_{a_{N-1}}^I f_{a_N}^I
```

The coefficients are computed directly without ever forming the full tensor.
Due to permutation symmetry of the generators, the last index is always
contracted with forces.

# Arguments
- `coeffs::AbstractVector{T}` — length `n_generators`, filled in-place
- `v::AbstractArray{T,3}` — shape `(dim, nat, n_configs)`, Cartesian displacements
- `f::AbstractArray{T,3}` — shape `(dim, nat, n_configs)`, Cartesian forces
- `generators::Vector{Generator{T,N,M}}` — compact generators (rank encoded in N)

# Example
```julia
generators = get_tensor_generators(symmetry_group, cell; rank=3)
coeffs = zeros(Float64, length(generators))
get_coefficients_from_ensemble!(coeffs, displacements, forces, generators)
# Reconstruct the full tensor (if needed):
# reconstruct_tensor!(tensor, generators, coeffs, cell)
```
"""
function get_coefficients_from_ensemble!(coeffs::AbstractVector{T},
    v::AbstractArray{T,3}, f::AbstractArray{T,3},
    generators::Vector{Generator{T,N,M}}) where {T,N,M}

    n_configs = size(v, 3)
    @assert size(v, 3) == size(f, 3) "v and f must have the same number of configurations"
    @assert length(coeffs) == length(generators)

    if isempty(generators)
        coeffs .= zero(T)
        return
    end

    coeffs .= zero(T)

    for I in 1:n_configs
        v_config = @view v[:, :, I]
        f_config = @view f[:, :, I]
        accumulate_ensemble_config!(coeffs, v_config, f_config, generators)
    end

    coeffs .*= -one(T) / n_configs
end
