abstract type ASRRule <: Function end

@doc raw"""
ASRConstraint!

Apply the ASR constraint to a rank-`rank` tensor of dimension `dimension`.
The application works as 

```julia
my_asr! = ASRconstraint!(dimension, 2)
my_tensor = rand(nat * dimension, nat * dimension)
my_asr!(my_tensor)
```

This will apply the ASR constraint to the tensor `my_tensor` in place.
"""
struct ASRConstraint! <: ASRRule
    dimension::Int
end
function (asr::ASRConstraint!)(vector::AbstractVector{T}) where T
    # Get the dimension of the tensor
    dimension = asr.dimension

    # Get the number of atoms
    nat = length(vector) ÷ dimension

    mytensor = reshape(vector, dimension, :)
    marginal = sum(mytensor, dims=2) ./ nat

    # Apply the ASR constraint
    for i in 1:nat
        @simd for k in 1:dimension
            index = dimension*(i-1) + k
            vector[index] -= marginal[k]
        end
    end
end

@doc raw"""
    (asr::ASRConstraint!)(matrix :: AbstractMatrix{T}; buffer=default_buffer(), differentiable=false) where T

Apply the ASR constraint to a rank-2 tensor of dimension `asr.dimension`.

The ASR is applied using the following formula to the $\Phi$ matrix:

```math
\Phi' = (I - \sum_t \left| t\right>\left< t\right|) \Phi(I - \sum_t \left| t\right>\left< t\right|)
```

where $\left |t\right>$ is the $t$-th global translation vector (1 for each dimension). 

The implementation follows the equation

```math
\Phi_{ij}^{\alpha\beta}' = \Phi_{ij}^{\alpha\beta} - \frac{1}{N_{\text{at}}}\sum_{tk} \Phi_{ik}^{\alpha t}\delta_{\beta t}
- \frac{1}{N_{\text{at}}}\sum_{tk} \Phi_{kj}^{t\beta}\delta_{\alpha t}
+ \frac{1}{N_{\text{at}}^2}\sum_{t_1t_2hk} \Phi_{hk}^{t_1t_2}\delta_{\alpha t_1}\delta_{\beta t_2}
```

which ensures that the translational invariance is mathematically preserved.

If differentiable is set to true, the function will not use bumper to allow memory allocation 
and differentiability of the function.
"""
function (asr::ASRConstraint!)(matrix :: AbstractMatrix{T}; buffer=default_buffer(), differentiable=false) where T
    nat = size(matrix, 1) ÷ asr.dimension
    nmodes = size(matrix, 1)
    
    # Get the last part
    @no_escape buffer begin
        if !differentiable
            trans = @alloc(T, nmodes, nmodes)
            tmp_mat = @alloc(T, nmodes, nmodes)
        else
            trans= zeros(T, nmodes, nmodes)
            tmp_mat = zeros(T, nmodes, nmodes)
        end
        trans .= 0
        for i in 1:nmodes
            trans[i, i] = 1
        end


        # Add the translational projector
        for i in 1:asr.dimension
            for h in 1:nmodes
               h_dim = (h-1) % asr.dimension + 1
               if h_dim != i
                   continue
               end

               for k in 1:nmodes
                   k_dim = (k-1) % asr.dimension + 1
                   if k_dim != i
                       continue
                   end

                   trans[h, k] -= 1 / nat
               end
           end
        end

        #TODO : This is not the most efficient way to do this
        #Probably we can exploit the sparsity of trans to speedup
        mul!(tmp_mat, matrix, trans)
        mul!(matrix, trans, tmp_mat)
        nothing
    end
end


@doc raw"""
    translation_mask!(mask::Vector{Bool}, pols::Matrix{T}, masses::Vector{T}; buffer=default_buffer()) where T

Identifies non-translational modes (e.g., optical modes) by checking
if the crystal's center of mass moves.

The function updates the `mask` in-place, setting `mask[i] = true` for
any mode `i` that does **not** move the center of mass.

## Mathematical and Physical Interpretation

This function correctly computes the (unnormalized) displacement of the
center of mass for each mode `i`.

1.  **Assumption:** The `pols` matrix contains the standard
    **mass-weighted** polarization eigenvectors, ``e_{j,\alpha}``.
    These relate to the real-space displacement vectors ``u_{j,\alpha}``
    by the formula:
    ``e_{j,\alpha} = \sqrt{m_j} \cdot u_{j,\alpha}``

2.  **Calculation:** The function computes a `ndims`-sized vector `dispv`
    for each mode. The calculation is:
    ``
    \texttt{dispv}_\alpha = \sum_j \left( \texttt{pols}_{j,\alpha}^{(i)} \cdot \sqrt{\texttt{masses}_j} \right)
    ``

3.  **Substitution:** By substituting the definition from (1), we see
    what `dispv` physically represents:
    ``
    \texttt{dispv}_\alpha = \sum_j \left( (\sqrt{m_j} \cdot u_{j,\alpha}^{(i)}) \cdot \sqrt{m_j} \right)
                         = \sum_j \left( m_j \cdot u_{j,\alpha}^{(i)} \right)
    ``
    This vector `dispv` is the component-wise sum of mass-weighted
    real-space displacements. This is equal to the total mass of the
    system times the displacement of its center of mass
    (``M_{\text{tot}} \cdot U_{\text{cm}}``).

4.  **Condition:** The check `max(abs.(dispv)...) <= 1e-6` is `true`
    if and only if `dispv` is the zero vector.
    * **Non-translational modes** (like optical modes, or acoustic
        modes at $\mathbf{q} \neq 0$) are *defined* by the
        property that they do not move the center of mass (``U_{\text{cm}} = 0``).
        For these modes, `dispv` is zero, the check passes,
        and `mask[i]` is set to `true`.
    * **Translational modes** (the 3 acoustic modes at $\mathbf{q} = 0$)
        *do* move the center of mass (``U_{\text{cm}} \neq 0``).
        For these modes, `dispv` is non-zero, the check fails,
        and `mask[i]` remains `false`.

## Arguments

- `mask::Vector{Bool}`: The boolean mask to be updated in-place.
  `mask[i]` will be set to `true` for non-translational modes.
- `pols::Matrix{T}`: The matrix of **mass-weighted** polarization
  eigenvectors, with modes as columns. Must have size
  `(Ndims * Natoms, Nmodes)`.
- `masses::Vector{T}`: A vector of the **atomic masses** for
  each atom, with size `(Natoms)`.
- `buffer`: (Optional) A `Bumper.jl` buffer for allocating
  the temporary `dispv` vector to avoid allocations.

!!! note "Finding Translational Modes"
    This function finds modes where `mask[i] = true` is `true` for
    **non**-translational modes. If you need a mask that is `true`
    for *only* the 3 translational modes, simply invert the result
    of this function:
    ```julia
    translation_mask!(non_trans_mask, pols, masses)
    trans_mask = .!non_trans_mask
    ```
"""
function translation_mask!(mask :: AbstractVector{Bool}, pols :: AbstractMatrix{T}, masses :: Vector{U}; buffer=default_buffer()) where {T, U}

    mask .= false
    nmod = length(mask)
    nat = length(masses)
    ndims = size(pols, 1) ÷ nat

    @no_escape buffer begin
        dispv = @alloc(T, ndims)

        for i in 1:nmod
            dispv .= 0
            for j in 1:nat 
                factor = sqrt(masses[j])
                @simd for k in 1:ndims
                    dispv[k] += pols[ndims * (j-1) + k, i] * factor
                end
            end

            if maximum(abs, dispv) <= 1e-6
                mask[i] = true
            end
        end
        nothing
    end
end


