@doc raw"""
    vector_r2q!(
        v_q :: Array{Complex{T}, 3},
        v_sc :: Matrix{T},
        q_tot :: Matrix{T})


Fourier transform a vector from real space and q space.

$$
v_k(\vec q) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{-i 2\pi \vec R\cdot \vec q} v_k(\vec R)
$$


## Parameters

- v_q : (n_configs, n_modes, nq) 
    The target vector in Fourier space.
- v_sc : (n_configs, n_modes)
    The original vector in real space
- q_tot : (n_dims, nq)
    The list of q vectors
- itau : (nat_sc)
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- R_lat : (3, nat_sc)
    The origin coordinates of the supercell in which the atom is
"""
function vector_r2q!(
        v_q :: Array{Complex{T}, 3},
        v_sc :: Matrix{T},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}
    ) where {T <: AbstractFloat, I <: Integer}

    nq = size(q, 2)
    n_random = size(v_sc, 1)
    n_dims = size(q, 1)
    nat_sc = size(v_sc, 2) ÷ n_dims
    nat = size(v_q, 2)

    v_q .= 0

    for jq ∈ 1:nq
        for k ∈ 1:nat_sc
            @views q_dot_R = q[:, jq]' * R_lat[:, k]
            exp_value = exp(- 1im * 2π * q_dot_R)

            for α in 1:n_dims
                index_sc = n_dims * (k - 1) + α
                index_uc = n_dims * (itau[k] - 1) + α
                @simd for i ∈ 1:n_random
                    v_q[i, index_uc, jq] += exp_value * v_sc[i, index_sc]
                end
            end
        end
    end

    v_q ./= √(nq)
end

@doc raw"""
    vector_q2r!(
        v_sc :: Matrix{T},
        v_q :: Array{Complex{T}, 3},
        q_tot :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}) where {T <: AbstractFloat, I <: Integer}


Fourier transform a vector from q space to real space.

$$
v_k(\vec R) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{+i 2\pi \vec R\cdot \vec q} v_k(\vec q)
$$


## Parameters


- v_sc : (n_configs, 3*nat_sc)
    The target vector in real space
- v_q : (n_configs, nq, 3*nat) 
    The original vector in Fourier space. 
- q_tot : (3, nq)
    The list of q vectors
- itau : (nat_sc)
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- R_lat : (3, nat_sc)
    The origin coordinates of the supercell in which the atom is
"""
function vector_q2r!(
        v_sc :: Matrix{T},
        v_q :: Array{Complex{T}, 3},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}
    ) where {T <: AbstractFloat, I <: Integer}

    nq = size(q, 2)
    n_random = size(v_sc, 1)
    n_dims = size(R_lat, 1)
    nat_sc = size(v_sc, 2) ÷ n_dims
    tmp_vector = zeros(Complex{T}, (n_random, n_dims*nat_sc))

    v_sc .= 0
    for jq ∈ 1:nq
        for k ∈ 1:nat_sc
            @views q_dot_R = q[:, jq]' * R_lat[:, k]
            exp_value = exp(1im * 2π * q_dot_R)

            for α in 1:n_dims
                index_sc = n_dims * (k - 1) + α
                index_uc = n_dims * (itau[k] - 1) + α
                @simd for i ∈ 1:n_random
                    tmp_vector[i, index_sc] += exp_value * v_q[i, index_uc, jq]
                end
            end
        end
    end

    v_sc .= real(tmp_vector)
    v_sc ./= √(nq)
end

