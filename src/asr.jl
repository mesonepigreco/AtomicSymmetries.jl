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
    translation_mask!(mask::Vector{Bool}, pols::Matrix{T}, masses::Vector{T}, n_dims :: Int; multiply :: Bool =true; buffer=default_buffer()) where T

Identifies translational modes by checking
if the crystal's center of mass moves.

The function updates the `mask` in-place, setting `mask[i] = true` for
any mode `i` that does move the center of mass.

## Arguments

- `mask::Vector{Bool}`: The boolean mask to be updated in-place.
  `mask[i]` will be set to `true` for translational modes.
- `pols::Matrix{T}`: The matrix of **mass-weighted** polarization
  eigenvectors, with modes as columns. Must have size
  `(Ndims * Natoms, Nmodes)`.
- `masses::Vector{T}`: A vector of the **atomic masses** for
  each atom, with size `(Natoms)`.
- `n_dims` :: Int : The number of spatial dimensions (e.g., 3 for 3D systems).
- `multiply`: (Optional) If `true` (default), the function computes
  `dispv` using the mass-weighted polarization vectors. If `false`,
  it computes `dispv` using the unweighted polarization vectors.
- `buffer`: (Optional) A `Bumper.jl` buffer for allocating
  the temporary `dispv` vector to avoid allocations.
"""
function translation_mask!(mask :: AbstractVector{Bool}, pols :: AbstractMatrix{T}, masses :: Vector{U}, n_dims :: Int; multiply=true, buffer=default_buffer()) where {T, U}

    mask .= true
    nmod = length(mask)
    n_modes = size(pols, 1)


    @no_escape buffer begin
        dispv = @alloc(T, n_dims)

        for i in 1:nmod
            dispv .= 0
            for j in 1:n_modes
                j_dim = (j-1) % n_dims + 1
                factor = multiply ? sqrt(masses[j]) : 1/sqrt(masses[j])
                dispv[j_dim] += pols[j, i] * factor
            end


            # Center of mass displacement of the mode
            if maximum(abs, dispv) <= 1e-6
                mask[i] = false
            end
        end
        nothing
    end
end


