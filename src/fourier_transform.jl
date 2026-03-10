
@doc raw"""
    vector_r2q!(
        v_q :: AbstractArray{Complex{T}, 3},
        v_sc :: AbstractMatrix{T},
        q_tot :: Matrix{T})
    vector_r2q!(v_q :: AbstractArray{Complex{T}, 2},
        v_sc :: AbstractVector{T},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}
    ) where {T <: AbstractFloat, I <: Integer}


Fourier transform a vector from real space and q space.

``
\displaystyle
v_k(\vec q) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{-i 2\pi \vec R\cdot \vec q} v_k(\vec R)
``

It works both on a single vector and on a series of vector. 
NOTE: In the latter case, the number
of configurations must be in the first column. 
This is not standard, 
but implemented in this way for performance reasons as
it is the most convenient memory rapresentation for vectorizing the
average calculation.

Notably, this convention introduces two main properties:

The ``\Gamma`` value of the fourier transform is not the average over the supercell of the same
quantity. If you want to obtain the average, you must divide by √nq (the number of q-points).

If the `R_lat` is not centered around zero, and the coordinates passed as `v_sc` are absolute values of positions,
then the ``\Gamma`` value of the fourier transform will be shifted by a total translation which is the average of the translations of the supercell lattice
vectors.

To avoid this behaviour (which is wrong), you can use the option `absolute_positions`,
which automatically rescales the `v_sc` to be coordinates relative
to the respective cell origin identified by `R_lat`.
Indeed, in this case, `R_lat` and `v_sc` must be of the same 
units, and coordinate types (you cannot mix crystalline and cartesian).


## Parameters

- `v_q` : (n_configs, 3nat, nq) 
    The target vector in Fourier space. Optionally, n_configs could be omitted if transforming only 1 vector
- `v_sc` : (n_configs, 3*nat_sc)
    The original vector in real space. Optionally, n_configs could be omitted if transforming only 1 vector
- `q_tot` : (3, nq)
    The list of q vectors
- `itau` : (nat_sc)
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- `R_lat` : (3, nat_sc)
    The origin coordinates of the supercell in which the atom is
- `absolute_positions` : Bool
    If true [default false], removes from v_sc the value of R_lat.
"""
function vector_r2q!(
        v_q :: AbstractArray{Complex{T}, 3},
        v_sc :: AbstractMatrix{T},
        q :: AbstractMatrix{T},
        itau :: AbstractVector{I},
        R_lat :: AbstractMatrix{T};
        absolute_positions :: Bool = false
    ) where {T <: AbstractFloat, I <: Integer}

    nq = size(q, 2)
    n_random = size(v_sc, 1)
    nat_sc = size(v_sc, 2) ÷ 3
    nat = size(v_q, 2)

    v_q .= 0


    for jq ∈ 1:nq
        for k ∈ 1:nat_sc
            @views q_dot_R = q[:, jq]' * R_lat[:, k]
            exp_value = exp(- 1im * 2π * q_dot_R)

            for α in 1:3
                index_sc = 3 * (k - 1) + α
                index_uc = 3 * (itau[k] - 1) + α
                @simd for i ∈ 1:n_random
                    # Remove the absolute position if needed
                    δ_value = v_sc[i, index_sc]
                    if absolute_positions
                        δ_value -= R_lat[α, k]
                    end

                    v_q[i, index_uc, jq] += exp_value * δ_value
                end
            end
        end
    end

    v_q ./= √nq
end
function vector_r2q!(
        v_q :: AbstractMatrix{Complex{T}},
        v_sc :: AbstractVector{T},
        q :: AbstractMatrix{T},
        itau :: AbstractVector{I},
        R_lat :: AbstractMatrix{T};
        kwargs...
    ) where {T <: AbstractFloat, I <: Integer}

    vector_r2q!(reshape(v_q, 1, size(v_q)...), reshape(v_sc, 1, size(v_sc)...), q, itau, R_lat; kwargs...)
end


@doc raw"""
    vector_q2r!(
        v_sc :: AbstractMatrix{T},
        v_q :: AbstractArray{Complex{T}, 3},
        q_tot :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T};
        absolute_positions :: Bool = false
        ) where {T <: AbstractFloat, I <: Integer}
    function vector_q2r!(
        v_sc :: AbstractVector{T},
        v_q :: AbstractMatrix{Complex{T}},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T};
        absolute_positions :: Bool = false
    ) where {T <: AbstractFloat, I <: Integer}


Fourier transform a vector from q space to real space.

``
\displaystyle
v_k(\vec R) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{+i 2\pi \vec R\cdot \vec q} v_k(\vec q)
``

It can be applied both to a single vector and in an ensemble.
NOTE: In the latter case, the configurations must be stored as the first index.
This choice is made for performance reason in computing averages (exploiting vectorization).


## Parameters


- `v_sc` : (n_configs, 3*nat_sc)
    The target vector in real space. Optionally, n_configs can be omitted
- `v_q` : (n_configs, nq, 3*nat) 
    The original vector in Fourier space. Optionally, n_configs can be omitted
- `q_tot` : (3, nq)
    The list of q vectors
- `itau : (nat_sc)`
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- `R_lat : (3, nat_sc)`
    The origin coordinates of the supercell in which the atom is
- `absolute_positions` : Bool
    If true, add the absolute position of the cell to the transformed v_sc
"""
function vector_q2r!(
        v_sc :: AbstractMatrix{T},
        v_q :: AbstractArray{Complex{T}, 3},
        q :: AbstractMatrix{T},
        itau :: AbstractVector{I},
        R_lat :: AbstractMatrix{T};
        absolute_positions :: Bool = false
    ) where {T <: AbstractFloat, I <: Integer}

    nq = size(q, 2)
    n_random = size(v_sc, 1)
    nat_sc = size(v_sc, 2) ÷ 3
    tmp_vector = zeros(Complex{T}, (n_random, 3*nat_sc))

    v_sc .= 0
    if absolute_positions
        for h in 1:nat_sc
            for k in 1:3
                δvalue = R_lat[k, h]
                index = 3*(h-1) + k
                v_sc[:, index] .+= δvalue
            end
        end
    end

    for jq ∈ 1:nq
        for k ∈ 1:nat_sc
            @views q_dot_R = q[:, jq]' * R_lat[:, k]
            exp_value = exp(1im * 2π * q_dot_R)

            for α in 1:3
                index_sc = 3 * (k - 1) + α
                index_uc = 3 * (itau[k] - 1) + α
                @simd for i ∈ 1:n_random
                    tmp_vector[i, index_sc] += exp_value * v_q[i, index_uc, jq]
                end
            end
        end
    end

    v_sc .+= real(tmp_vector)
    v_sc ./= √(nq)
end

function vector_q2r!(
        v_sc :: AbstractVector{T},
        v_q :: AbstractMatrix{Complex{T}},
        q :: AbstractMatrix{T},
        itau :: AbstractVector{I},
        R_lat :: AbstractMatrix{T};
        kwargs...
    ) where {T <: AbstractFloat, I <: Integer}

    vector_q2r!(reshape(v_sc, 1, size(v_sc)...), reshape(v_q, 1, size(v_q)...), q, itau, R_lat; kwargs...)
end


@doc raw"""
    matrix_r2q!(
        matrix_q :: Array{Complex{T}, 3},
        matrix_r :: AbstractMatrix{T},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T})

Fourier transform a matrix from real to q space

```math
M_{ab}(\vec q) = \sum_{\vec R} e^{2\pi i \vec q\cdot \vec R}\Phi_{a;b + \vec R}
```

Where ``\Phi_{ab}`` is the real space matrix, the ``b+\vec R`` indicates the corresponding atom in the supercell displaced by ``\vec R``. 


## Parameters

- `matrix_q` : (3nat, 3nat, nq) 
    The target matrix in Fourier space.
- `matrix_r` : (3*nat_sc, 3*nat)
    The original matrix in real space (supercell)
- `q_tot` : (3, nq)
    The list of q vectors
- `itau` : (nat_sc)
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- `R_lat` : (3, nat_sc)
    The origin coordinates of the supercell in which the corresponding atom is
"""
function matrix_r2q!(
        matrix_q :: AbstractArray{Complex{T}, 3},
        matrix_r :: AbstractMatrix{T},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}; buffer = default_buffer()) where {T, I<: Integer}
    nq = size(q, 2)
    ndims = size(q, 1)
    nat_sc = size(matrix_r, 1) ÷ ndims
    nat = size(matrix_q, 1) ÷ ndims

    matrix_q .= T(0.0) 

    @no_escape buffer begin
        ΔR⃗ = @alloc(T, ndims)
        tmp_mat = @alloc(Complex{T}, ndims, ndims)

        phase_i = Complex{T}(-2π * 1im)

        for iq in 1:nq
            for k_i in 1:nat
                @simd for h_i in 1:nat_sc
                    @views ΔR⃗ .= R_lat[:, k_i]
                    @views ΔR⃗ .-= R_lat[:, h_i]
                    @views q_dot_R = ΔR⃗' * q[:, iq]

                    h_i_uc = itau[h_i]

                    exp_factor = exp(phase_i * q_dot_R)
                    @views tmp_mat .= matrix_r[(ndims*(h_i - 1) + 1 : ndims * h_i), (ndims*(k_i - 1) +1 : ndims*k_i)]
                    tmp_mat .*= exp_factor 

                    @views matrix_q[(ndims*(h_i_uc - 1) + 1 : ndims * h_i_uc), (ndims*(k_i - 1) +1 : ndims*k_i), iq] .+= tmp_mat 
                    #matrix_q[(ndims*(h_i_uc - 1) +1 : ndims * h_i_uc), (ndims*(k_i - 1) + 1 : ndims*k_i), iq] .*= exp_factor
                end
            end
        end
        nothing
    end
end

@doc raw"""
    matrix_q2r!(
        matrix_r :: AbstractMatrix{T},
        matrix_q :: Array{Complex{T}, 3},
        q :: Matrix{T},
        itau :: Vector{Int},
        R_lat :: Matrix{T})

Fourier transform a matrix from q space into r space

``
\displaystyle
\Phi_{ab} = \frac{1}{N_q} \sum_{\vec q}
M_{ab}(\vec  q) e^{2i\pi \vec q\cdot[\vec R(a) - \vec R(b)]}
``

Where ``\Phi_{ab}`` is the real space matrix, ``M_{ab}(\vec q)`` is the q space matrix.


## Parameters


- matrix_r : (3*nat_sc, 3*nat)
    The target matrix in real space (supercell). If the second dimension is 3nat_sc, we also apply the translations
- matrix_q : (3nat, 3nat, nq) 
    The original matrix in Fourier space.
- q_tot : (3, nq)
    The list of q vectors
- itau : (nat_sc)
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- R_lat : (3, nat_sc)
    The origin coordinates of the supercell in which the corresponding atom is
- translations : Vector{Vector{Int}}
    The itau correspondance for each translational vector. Its size must be equal to the number of q point and
    contain all possible translations. This can be obtained from the `get_translations` subroutine. 
"""
function matrix_q2r!(
        matrix_r :: AbstractMatrix{T},
        matrix_q :: Array{Complex{T}, 3},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}; translations :: Union{Nothing, AbstractVector} = nothing, buffer = default_buffer()) where {T, I <: Integer}
    nq = size(q, 2)
    ndims = size(q, 1)
    nat_sc = size(matrix_r, 1) ÷ ndims
    nat = size(matrix_q, 1) ÷ ndims

    matrix_r .= T(0.0) 

    apply_translations = false
    if size(matrix_r, 2) > nat*ndims
        apply_translations = true
        if size(matrix_r, 2) != nat_sc*ndims
            error("Error, dimension mismatch in matrix_r: $(size(matrix_r))")
        end

        # Check if the translations are provided
        if translations == nothing
            error("Error in matrix_q2r! : if the size of the matrix_r is a square, then it is required to provide the translations.")
        end

        # Check if the translations have the correct lenght
        @assert length(translations) == nq "Error, the number of translations ($(length(translations))) must be equal with the number of q-points ($nq)"
    end

    @no_escape buffer begin
        ΔR⃗ = @alloc(T, ndims)

        phase_i = Complex{T}(2π * 1im)
        tmp_matrix = @alloc(Complex{T}, ndims, ndims)

        for iq in 1:nq

            for k_i in 1:nat
                @simd for h_i in 1:nat_sc
                    @views ΔR⃗ .= R_lat[:, k_i]
                    @views ΔR⃗ .-= R_lat[:, h_i]
                    @views q_dot_R = ΔR⃗' * q[:, iq]

                    h_i_uc = itau[h_i]

                    exp_factor = exp(phase_i * q_dot_R)
                    #TODO: createa temporaney structure before adding the exponential otherwise itis not real
                    @views tmp_matrix .= matrix_q[(ndims*(h_i_uc - 1) +1 : ndims * h_i_uc), (ndims*(k_i - 1)+1 : ndims*k_i), iq]
                    tmp_matrix .*= exp_factor / nq
                    @views matrix_r[(ndims*(h_i - 1) +1 : ndims * h_i), (ndims*(k_i - 1) + 1 : ndims*k_i)] .+= real(tmp_matrix)
                end
            end
        end

        # Check if we need to apply the translations
        if apply_translations
            new_tmp_matrix = @alloc(T, nat_sc*ndims, nat_sc*ndims)
            new_tmp_matrix .= T(0)
            for trans in translations
                for i in 1:nat_sc
                    i_t = trans[i]
                    for j in 1:nat_sc
                        j_t = trans[j]
                        @views new_tmp_matrix[(ndims*(j-1)+1:ndims*j), (ndims*(i-1)+1:ndims*i)] .+=  
                            matrix_r[(ndims*(j_t-1)+1:ndims*j_t), (ndims*(i_t-1)+1:ndims*i_t)]
                    end
                end
            end

            matrix_r .= new_tmp_matrix
        end
        nothing
    end
end


@doc raw"""
    shift_position_origin!(r_vector::AbstractVector{Complex{T}}, cell::AbstractMatrix{T}, R_lat::AbstractMatrix{T}; buffer) where T
    shift_position_origin!(r_vectors::AbstractMatrix{Complex{T}}, cell::AbstractMatrix{T}, R_lat::AbstractMatrix{T}; buffer) where T

Shifts a set of position vectors (`r_vector` or `r_vectors`) such that the average position (center of mass) of the reference lattice (`R_lat`) is moved to the origin. The function operates in-place.
This can be employed after performing a Fourier Transform of absolute positions, to correct for non centered `R_lat`.

This is typically used to remove the lattice translation from displacement vectors before calculating quantities that are invariant to rigid body translations.

The average position vector of the reference lattice $\mathbf{r}_{\text{avg}}$ is calculated in Cartesian coordinates, and then subtracted from the input vectors $\mathbf{r}$.

```math
\mathbf{r}_{\text{avg}} = \frac{1}{N_{\text{atoms}}} \sum_{k=1}^{N_{\text{atoms}}} \mathbf{r}_{\text{lat}, k}
```

The position vectors are then shifted:
```math
\mathbf{r}' = \mathbf{r} - \mathbf{r}_{\text{avg}}
```

## Arguments

- `r_vector`/`r_vectors` : The position vector(s) (of type `Complex{T}`) to be shifted in-place. 
- `cell` : The $N_{dims} \times N_{dims}$ matrix defining the unit cell (lattice vectors). 
- `R_lat` : The reference lattice positions (fractional/lattice coordinates, ``N_{dims} \times N_{atoms}``) whose average position determines the shift. 
- `buffer` : An optional buffer for temporary memory allocation (Bumer.jl)

## Details on `r_vectors` Matrix Method

The matrix method assumes that the input `r_vectors` has the dimensions ``(N_{\text{configs}}, N_{\text{dims}} \times N_{\text{atoms,v}})``, where ``N_{\text{configs}}`` is the number of configurations, and ``N_{\text{atoms,v}}`` is the number of atoms in the vectors being shifted.

The shift is applied simultaneously to all configurations using broadcasting (`@views ... .-=`).


## Example

An example of usage after the fourier transform

```julia
# Define positions_r as a n_configs, 3n_atoms_sc vector
# Perform the fourier transform in q space
# Here, we assume that R_lat and q_points are expressed in crystal coordinates.
# Otherwise, just pass the identity to the cell below.
vector_r2q!(positions_q, positions_r, q_points, itau, R_lat)

# Remove the translations
shift_position_origin!(positions_q, cell, R_lat)
```
"""
function shift_position_origin!(r_vector :: AbstractVector{U}, cell :: AbstractMatrix{T}, R_lat :: AbstractMatrix{T}; buffer = default_buffer()) where {T, U <: Union{T, Complex{T}}}

    n_dims = size(R_lat, 1)
    n_atoms = size(R_lat, 2)
    n_atoms_v = length(r_vector) ÷ n_dims

    @no_escape buffer begin
        r_lat_cart = @alloc(T, size(R_lat)...)
        r_lat_avg = @alloc(T, n_dims)

        get_cartesian_coords!(r_lat_cart, R_lat, cell)

        # Average the displacement
        r_lat_avg .= zero(T)
        for i in 1:n_dims
            @simd for k in 1:n_atoms
                r_lat_avg[i] += r_lat_cart[i, k]
            end
        end
        r_lat_avg ./= n_atoms

        for k in 1:n_atoms_v
            @simd for i in 1:n_dims
                r_vector[n_dims * (k-1) + i] -= r_lat_avg[i]
            end
        end
    end
    r_vector
end
function shift_position_origin!(r_vectors :: AbstractMatrix{U}, cell :: AbstractMatrix{T}, R_lat :: AbstractMatrix{T}; buffer = default_buffer()) where {T, U <: Union{T, Complex{T}}}

    n_dims = size(R_lat, 1)
    n_atoms = size(R_lat, 2)
    n_configs = size(r_vectors, 1)
    n_atoms_v = size(r_vectors) ÷ n_dims

    @no_escape buffer begin
        r_lat_cart = @alloc(T, size(R_lat)...)
        r_lat_avg = @alloc(T, n_dims)

        get_cartesian_coords!(r_lat_cart, R_lat, cell)

        # Average the displacement
        r_lat_avg .= zero(T)
        for i in 1:n_dims
            @simd for k in 1:n_atoms
                r_lat_avg[i] += r_lat_cart[i, k]
            end
        end
        r_lat_avg ./= n_atoms

        for k in 1:n_atoms_v
            @simd for i in 1:n_dims
                r_remove = r_lat_avg[i]
                index =  n_dims * (k-1) + i
                @views r_vectors[:, index] .-= r_remove
            end
        end
    end
    r_vectors
end


