@doc raw"""
    SymmetriesQSpace(symmetries :: Symmetries{T}, q_points :: AbstractMatrix{T}) :: SymmetriesQSpace{T} where T

    struct SymmetriesQSpace{T} <: GenericSymmetries where T
        symmetries :: Symmetries{T}
        irt_q :: Vector{Vector{Int}}
        minus_q_index :: Vector{Int}
    end

This structure contains the information to perform the symmetrization of a dynamical matrix directly in q space.

**Note that the `q_points` needs to be in crystal coordinates**,

and the symmetries must be of the primitive cell.
"""
struct SymmetriesQSpace{T} <: GenericSymmetries 
    symmetries :: Symmetries{T}
    q_points :: Matrix{T}
    irt_q :: Vector{Vector{Int}}
    minus_q_index :: Vector{Int}
end
function SymmetriesQSpace(symmetries :: Symmetries{T}, q_points :: AbstractMatrix{T}; buffer = default_buffer()) :: SymmetriesQSpace{T} where T
    n_symmetries = length(symmetries)
    n_q = size(q_points, 2)
    n_dims = size(q_points, 1)

    my_q_points = zeros(T, n_dims, n_q)
    my_q_points .= q_points

    irt_q = Vector{Vector{Int}}(undef, n_symmetries)
    for i in 1:n_symmetries
        irt_q[i] = Vector{Int}(undef, n_q)
        get_irt_q!(irt_q[i], q_points, symmetries.symmetries[i]; buffer = buffer)
    end

    minus_q_index = zeros(Int, n_q)
    get_minus_q!(minus_q_index, q_points)

    SymmetriesQSpace(symmetries, my_q_points, irt_q, minus_q_index)
end

Base.isempty(x :: SymmetriesQSpace) = isempty(x.symmetries)
Base.length(x :: SymmetriesQSpace) = length(x.symmetries)
function Base.getindex(x :: SymmetriesQSpace, k) 
    return x.symmetries.symmetries[k]
end


@doc raw"""
    check_symmetries(q_symmetries :: SymmetriesQSpace{T}, n_atoms :: Int) :: Bool

Check if the q_symmetries has been correctly initialized in the primitive cell.

Essentially, this subroutine checks the atomic correspondance by symmetry
and spots if there are atoms outside the primitive cell (whose index is
above `n_atoms`).

## Parameters

- `q_symmetries` : The symmetries in q space
- `n_atoms` : The number of atoms in the primitive cell

## Returns

`true` if no contraddiction have been detected,
`false` otherwise.
"""
function check_symmetries(q_symmetries :: SymmetriesQSpace{T}, n_atoms :: Int) :: Bool where T
    for i in 1:length(q_symmetries)
        n_length = length(q_symmetries.symmetries.irt[i])

        if n_length > n_atoms
            return false
        end
        for k in 1:n_length
            if q_symmetries.symmetries.irt[i][k] > n_atoms
                return false
            end
        end
    end
    return true
end



@doc raw"""
    apply_symmetry_vectorq!(target_vector :: AbstractMatrix{Complex{T}}, original_vector :: AbstractMatrix{Complex{T}}, symmetry_operation :: AbstractMatrix{U}, irt :: Vector{Int}, irt_q:: AbstractVector{Int})


Apply the symmetry on the original vector in q space

## Parameters

- `target_vector` : The result (modified inplace) (3n x nq)
- `original_vector` : The original vector (3n x nq)
- `symmetry_operation` : The 3x3 symmetry 
- `irt` : The atom-atom association by symmetry
- `irt_q` : The q-q association by symmetry
"""
function apply_symmetry_vectorq!(target_vector :: AbstractMatrix{Complex{T}}, original_vector :: AbstractMatrix{Complex{T}},
        symmetry_operation :: AbstractMatrix{U}, irt :: AbstractVector{Int}, irt_q :: AbstractVector{Int}) where {T, U}
    # Apply symmetries 
    nq = length(irt_q)
    n_dims = size(symmetry_operation, 1)
    n_atoms = size(target_vector, 1) ÷ n_dims


    for iq in 1:nq
        jq = irt_q[iq]

        for i in 1:n_atoms
            j = irt[i]

            @views mul!(target_vector[n_dims * (j - 1) + 1: n_dims * j, jq], 
                        symmetry_operation, 
                        original_vector[n_dims * (i - 1) + 1: n_dims * i, iq],
                        T(1.0), T(1.0))
        end
    end
end

@doc raw"""
    apply_symmetry_matrixq!(target_matrix :: AbstractArray{Complex{T}, 3},
        original_matrix :: AbstractArray{Complex{T}, 3},
        sym :: AbstractMatrix{U},
        irt :: AbstractVector{Int},
        irt_q :: AbstractVector{Int},
        unit_cell_translations :: AbstractMatrix{T},
        ; buffer = default_buffer()) where {T, U}


Apply the symmetry on the matrix in q space
This subroutine assumes the convention that the phase factor is for each supercell, not atoms.
In other words, all the atoms coordinates are computed from the same origin of the supercell they are associated with.

## Parameters

- `target_vector` : The result (modified inplace) (3n x nq)
- `original_vector` : The original vector (3n x nq)
- `symmetry_operation` : The 3x3 symmetry 
- `irt` : The atom-atom association by symmetry
- `irt_q` : The q-q association by symmetry
- `unit_cell_translations` : The translation vectors to move the transformed atom in the primitive cell
- `buffer` : The Bumper.jl buffer for caching memory allocations [Optional]
"""
function apply_symmetry_matrixq!(target_matrix :: AbstractArray{Complex{T}, 3},
        original_matrix :: AbstractArray{Complex{T}, 3},
        sym :: AbstractMatrix{U},
        irt :: AbstractVector{Int},
        irt_q :: AbstractVector{Int},
        unit_cell_translations :: AbstractMatrix{T},
        q_points :: AbstractMatrix{T}
        ; buffer = default_buffer()) where {T, U}

    nq = size(target_matrix, 3)
    n_dims = size(sym, 1)
    n_atoms = size(target_matrix, 1) ÷ n_dims


    @no_escape buffer begin
        work = @alloc(Complex{T}, n_dims, n_dims)
        δt = @alloc(T, n_dims)
        for iq in 1:nq
            iq_s = irt_q[iq]
            for i ∈ 1:n_atoms
                i_s = irt[i]
                for j in 1:n_atoms 
                    j_s = irt[j]
                    @views δt .= unit_cell_translations[:, i_s] .- unit_cell_translations[:, j_s]
                    @views q_dot_t = dot(q_points[:, iq], δt)
                    @views phase_factor = exp(1im * 2π * q_dot_t)

                    @views mul!(work, 
                                original_matrix[n_dims*(i_s-1) + 1: n_dims*i_s, n_dims*(j_s-1)+1 : n_dims*j_s, iq_s], 
                                sym, T(1.0), T(0.0))
                    @views mul!(target_matrix[n_dims*(i-1) + 1: n_dims*i, n_dims*(j - 1) + 1: n_dims*j, iq], 
                        sym', work, phase_factor, 1.0)
                end
            end
        end
        nothing
    end
end

@doc raw"""
    get_irt_q!(irt_q :: AbstractVector{Int}, q_points :: AbstractVector{T}, sym_mat :: AbstractMatrix)

Get the correspondance ``q' = S_\text{recip} q`` on the provided q grid.
Always assume everything is in crystal coordinates.

Note that in reciprocal space (crystal coordinates) the symmetry operation is the inverse transpose.

```math
S_\text{recip} = (S_\text{direct})^{-T}
```
The provided `sym_mat` is assumed to be in direct space.

This is needed for the correct application of the symmetries
"""
function get_irt_q!(irt_q :: AbstractVector{Int}, q_points :: AbstractMatrix{T}, sym_mat :: AbstractMatrix{U}; buffer = default_buffer()) where {T, U}
    nq = size(q_points, 2)
    ndims = size(q_points, 1)
    @no_escape buffer begin
        tmpvector = @alloc(T, ndims)
        tmp2 = @alloc(T, ndims)
        sym_rec = @alloc(U, ndims, ndims)
        sym_rec .= inv(sym_mat)'
        
        for i in 1:nq
            @views mul!(tmpvector, sym_rec, q_points[:, i])

            # Check the closest q point
            min_distance = T(Inf)
            min_index = 0
            for j in 1:nq
                @views tmp2 .= q_points[:, j] - tmpvector
                tmp2 .-= round.(tmp2)

                distance = sum(abs2, tmp2) 
                if distance < min_distance
                    min_index = j
                    min_distance = distance
                end
            end

            irt_q[i] = min_index
        end
        nothing
    end

    # Check if there are errors
    if !allunique(irt_q)
        error_msg = """
Error while inizialising the symmetries in q space.
In particular, the symmetry operation 
S = $sym_mat

Brings two distinct q points into the same vector.
This either occurs if S has a nonzero kernel 
(then it is a nonvalid symmetry operation), 
or if the q points are not properly organized in a grid
in crystal coordinates (most likely).

Please, check if the q points are correctly in a uniform
grid in fractional coordinates of the Brilluin zone.

q_points = $(q_points')

irt_q = $irt_q
(Here rows are the vectors)
"""
        error(error_msg)
    end
end

@doc raw"""
    get_minus_q!(minus_q_index :: AbstractVector{Int}, q_points :: AbstractMatrix{T}; buffer = default_buffer()) where T

Identify for each q point what is the corresponding -q:

``
\vec q \longrightarrow -\vec q + \vec G
``

where ``\vec G`` is a reciprocal vector. Since this is done in crystal coordinates``\vec G`` are all possible integers.
"""
function get_minus_q!(minus_q_index :: AbstractVector{Int}, q_points :: AbstractMatrix{T}; buffer = default_buffer()) where T
    @no_escape buffer begin
        ndims = size(q_points, 1)
        nq = size(q_points, 2)
        tmpvector = @alloc(T, ndims)
        tmp2 = @alloc(T, ndims)
        
        for i in 1:nq
            @views tmpvector .= -q_points[:, i]

            # Check the closest q point
            min_distance = T(Inf)
            min_index = 0
            for j in 1:nq
                @views tmp2 .= q_points[:, j] - tmpvector
                tmp2 .-= round.(tmp2)
                distance = sum(abs2, tmp2) 
                if distance < min_distance
                    min_index = j
                    min_distance = distance
                end
            end

            # TODO: Check if it is this or the opposite
            minus_q_index[i] = min_index
        end
        nothing
    end
end

@doc raw"""
    symmetrize_vector_q!(target_gamma :: AbstractVector{T}, original_q :: AbstractArray{Complex{T}, 2}, symmetries :: Symmetries, irt_q :: Vector{Vector{Int}}; buffer = default_buffer() where T

Impose the symmetrization of a vector in q space.
Since the symmetrization also imposes translational symmetries, the result is always a vector only at gamma.

The symmetrized vector is supposed to be a displacement (so no translations are applied)

NOTE: The provided vector must be in crystal coordinates
To symmetrize a vector incartesian coordinates, see the routine `symmetrize_vector_cartesian_q!`.


## Parameters

- `target_gamma` : The `n_at * n_dims` output symmetrized vector at ``\Gamma``
- `original_q` : The original vector in q-space of size `nat*n_dims, nq`
- `symmetries` : The symmetry group
- `irt_q` : A vector (one for each symmetry) of the correspondances of q points. For each symmetry can be obtained from `get_irt_q!`
- `gamma_index` : Specify which q vector is ``\Gamma``. If not specified, it is assumed to be the first one
"""
function symmetrize_vector_q!(target_gamma :: AbstractVector{T}, original_q :: AbstractArray{Complex{T}, 2}, sym :: Symmetries, irt_q :: Vector{Vector{Int}}; gamma_index=1, buffer = default_buffer()) where {T}

    n_modes = size(original_q, 1)
    n_q = size(original_q, 2)

    @assert gamma_index <= n_q "Error, the number of q points ($n_q) cannot be lower than the index of Γ ($gamma_index)"

    @no_escape buffer begin
        tmp_vector = @alloc(Complex{T}, n_modes, n_q)
        tmp_vector .= zero(T)

        for i in 1:length(sym)
            sym_mat = sym.symmetries[i]
            irt = sym.irt[i]
            q_irt = irt_q[i]

            apply_symmetry_vectorq!(tmp_vector, original_q, sym_mat, irt, q_irt)
        end

        tmp_vector ./= length(sym)
        @views broadcast!(real, target_gamma, tmp_vector[:, gamma_index])
        nothing
    end
end

@doc raw"""
    symmetrize_vector_cartesian_q!(vector_q_cart:: AbstractArray{Complex{T}, 2}, cell :: Matrix{T}, symmetries :: SymmetriesQSpace; buffer = default_buffer()) where {T}

Perform the symmetrization of a vector (overwriting it) in cartesian coordinates.
This is the go-to subroutine for performing symmetrization of vectors in q space.


## Parameters

- `vector_q_cart` : in-place symmetrize vector (q-space, cartesian coordinates)
- `cell` : 3x3 matrix of the primitive cell (column-based)
- `symmetries` : Symmetries in Q space
- `buffer` : Optional, Bumper stack buffer (caching)

"""
function symmetrize_vector_cartesian_q!(vector_q_cart:: AbstractArray{Complex{T}, 2}, cell :: Matrix{T}, symmetries :: SymmetriesQSpace; buffer = default_buffer()) where {T}
    # 
    n_modes = size(vector_q_cart, 1)
    n_q = size(vector_q_cart, 2)
    n_dims = size(cell, 1)
    n_atoms = n_modes ÷ n_dims

    # Check if it is coherent
    if !check_symmetries(symmetries, n_atoms)
        error("Error, the symmetries in q space must be initialized on the primitive cell!")
    end

    @no_escape buffer begin
        vector_cryst = @alloc(T, n_modes)
        tmp_vect = @alloc(T, n_modes)
        @views broadcast!(real, tmp_vect, vector_q_cart[:, 1])

        get_crystal_coords!(reshape(vector_cryst, n_dims, :),
                            reshape(tmp_vect, n_dims, :),
                            cell;
                            buffer=buffer)
        vector_q_cart[:, 1] .= vector_cryst
        tmp_vect .= 0


        symmetrize_vector_q!(tmp_vect, vector_q_cart, 
                             symmetries.symmetries,
                             symmetries.irt_q;
                             buffer=buffer)

        # Convert back in cartesian space
        get_cartesian_coords!(reshape(vector_cryst, n_dims, :),
                              reshape(tmp_vect, n_dims, :),
                              cell)

        vector_q_cart[:, 1] .= vector_cryst
        nothing
    end
end


@doc raw"""
    symmetrize_matrix_q!(target_q :: AbstractArray{Complex{T}, 3}, original_q :: AbstractArray{Complex{T}, 3}, symmetries :: Symmetries, irt_q :: Vector{Vector{Int}}; buffer = default_buffer() where T
    symmetrize_matrix_q!(matrix_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T
    symmetrize_matrix_q!(target_q :: AbstractArray{T, 3}, original_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T


Impose the symmetrization of a dynamical matrix in q space.
The matrix must be in crystal coordinates.


## Parameters

- `target_q` : The symmetrized matrix of size `n_modes, n_modes, nq` (modified in-place).
- `original_q` : The original matrix in q-space of size `n_modes, n_modes, nq`. It could be the same as target_q
- `symmetries` : The symmetry group
- `irt_q` : A vector (one for each symmetry) of the correspondances of q points. For each symmetry can be obtained from `get_irt_q!`
- `unit_cell_translations` :: Vector{Matrix{T}} : The translations of the unit cell to bring back the atoms in the primitive cell after the symmetry operation. Each vector elements corresponds to one symmetry operation, then the matrix is a n_dims x n_atoms translation. This is usually the same as the content of `symmetries.unit_cell_translations`.
- `minus_q_index` : A vector containing for each `q` the corresponding ``\vec {q'} = -\vec q + \vec G``, where ``\vec G`` is a generic reciprocal lattice vector.
- `q_points` : The vector containing the actual q points
"""
function symmetrize_matrix_q!(target_q :: AbstractArray{Complex{T}, 3}, original_q :: AbstractArray{Complex{T}, 3}, symmetries :: Symmetries, irt_q :: Vector{Vector{Int}}, unit_cell_translations :: Vector{Matrix{T}}, minus_q_index::Vector{Int}, q_points :: AbstractMatrix{T}; buffer = default_buffer())  where T

    n_modes = size(original_q, 1)
    n_q = size(original_q, 3)

    @no_escape buffer begin
        tmp_matrix = @alloc(Complex{T}, n_modes, n_modes, n_q)
        tmp_matrix .= 0.0
        target_q .= Complex{T}(0.0)

        for i in 1:length(symmetries)
            sym_mat = symmetries.symmetries[i]
            irt = symmetries.irt[i]
            q_irt = irt_q[i]

            #println("Applying symmetry $i")
            #println("Symmetry matrix:")
            #println(sym_mat)
            #@show q_irt
            #@show unit_cell_translations[i]

            apply_symmetry_matrixq!(tmp_matrix, original_q, sym_mat, irt, q_irt, unit_cell_translations[i], q_points; buffer=buffer)
        end

        tmp_matrix ./= length(symmetries)

        # Apply the hermitianity
        for iq in 1:n_q
            for h in 1:n_modes
                for k in 1:n_modes
                    target_q[k,h, iq] = tmp_matrix[k, h, iq]
                    target_q[k,h, iq] += conj(tmp_matrix[h, k, iq])
                end
            end
        end
        target_q ./= T(2)


        ## Apply the time-reversal symmetry
        #tmp_matrix .= target_q
        #for iq in 1:n_q
        #    @views target_q[:, :, iq] .+= conj.(tmp_matrix[:, :, minus_q_index[iq]]')
        #end
        #target_q ./= T(2)
        nothing
    end
end
function symmetrize_matrix_q!(target_q :: AbstractArray{Complex{T}, 3}, original_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T
    symmetrize_matrix_q!(target_q, original_q, q_symmetries.symmetries, q_symmetries.irt_q, q_symmetries.symmetries.unit_cell_translations, q_symmetries.minus_q_index, q_symmetries.q_points; buffer=buffer)
end
function symmetrize_matrix_q!(matrix_q :: AbstractArray{Complex{T}, 3}, q_symmetries :: SymmetriesQSpace; buffer = default_buffer())  where T
    @no_escape buffer begin
        target = @alloc(Complex{T}, size(matrix_q)...)
        target .= zero(T)
        symmetrize_matrix_q!(target, matrix_q, q_symmetries; buffer=buffer)
        matrix_q .= target
        nothing
    end
end


@doc raw"""
    symmetrize_matrix_cartesian_q!(matrix_q :: AbstractArray{Complex{T}, 3}, cell :: Matrix{T}, q_symmetries :: SymmetriesQSpace; buffer=default_buffer()) where T

Enforce the symmetries on the provided matrix (q-space), modifying it in-place.
The provided matrix must be in Cartesian Coordinates.


## Parameters

- `matrix_q` : The matrix to be symmetrized. Size (nmodes, nmodes, nq)
- `cell` : The 3x3 primitive cell (column ordering)
- `q_symmetries` : The symmetry group (q-space)
- `buffer` : Optional, stack for Bumper to cache allocations.
"""
function symmetrize_matrix_cartesian_q!(matrix_q :: AbstractArray{Complex{T}, 3}, cell :: Matrix{T}, q_symmetries :: SymmetriesQSpace; buffer=default_buffer()) where T

    n_q = size(matrix_q, 3)

    @no_escape buffer begin
        matrix_cryst_q = @alloc(Complex{T}, size(matrix_q)...)
        matrix_cryst_q .= 0

        # Convert in crystal coordinates
        cart_cryst_matrix_conversion!(matrix_cryst_q,
                                      matrix_q,
                                      cell;
                                      cart_to_cryst = true,
                                      buffer=buffer)


        # Perform the symmetrization
        symmetrize_matrix_q!(matrix_cryst_q, q_symmetries; buffer)

        # Convert back in cartesian coordinates
        cart_cryst_matrix_conversion!(matrix_q,
                                      matrix_cryst_q,
                                      cell;
                                      cart_to_cryst = false,
                                      buffer=buffer)

        nothing
    end
end


@doc raw"""
    get_R_lat!(R_lat :: Matrix{T}, primitive_coords :: Matrix{T}, supercell_coords :: Matrix{T})

Get the `R_lat` parameter to perform the fourier transform.

## Parameters

- `R_lat` the result lattice vectors, modified inplace
- `primitive_coords` : The coordinates in the primitive cell
- `supercell_coords` : The coordinates in the supercell
- `itau` : The correspondence for each atom in the supercell with the respective atom in the primitive cell

"""
function get_R_lat!(R_lat :: AbstractMatrix{T}, primitive_coords :: AbstractMatrix{T}, supercell_coords :: AbstractMatrix{T}, itau :: Vector{I}) where {T, I <: Integer}
    nat_uc = size(primitive_coords, 2)
    nat_sc = size(supercell_coords, 2)

    @simd for i in 1:nat_sc
        @views R_lat[:, i] .= supercell_coords[:, i]
        @views R_lat[:, i] .-= primitive_coords[:, itau[i]]
    end
end



@doc raw"""
    get_supercell(q_points::AbstractMatrix{T}, cell :: AbstractMatrix{T}) :: Vector{Int}
    get_supercell(q_points::AbstractMatrix{T}) :: Vector{Int}
    get_supercell!(supercell::AbstractVector{I}, q_points::AbstractMatrix{T}, cell :: AbstractMatrix{T}) where {T, I<:Integer}
    get_supercell!(supercell::AbstractVector{I}, q_points::AbstractMatrix{T}) where {T, I<:Integer}

Calculates the minimum supercell dimensions required to fold a set of commensurate
wave vectors (`q_points`) back to the Gamma point (Γ) of the Brillouin zone.

The resulting supercell dimension \$S_i\$ for each spatial direction is determined
by the reciprocal of the smallest non-zero q-point component in that direction.
For commensurate grids, this is mathematically equivalent to:

```math
S_i = \text{round} \left( \frac{1}{\min(|q_{i}|)} \right)
```

This function is primarily used when performing calculations in a real-space supercell
that is commensurate with the input k-point (or q-point) grid.

## Arguments

- `q_points`: A 2D matrix where **each column** represents a q-point vector, and
  **each row** corresponds to a dimension (x, y, z).
- `supercell`: An pre-allocated integer vector to store the result (used by the
  `get_supercell!` in-place version).
- `cell` : The primitive cell. If not provided, the code assumes that q_points are 
    in fractional coordinates.

## Important Note on Coordinates

The `q_points` **must** be provided in **crystal coordinates** (fractional coordinates)
if `cell` is not provided. 

# Example

```julia
# 3 dimensions, 2 q-points
q_points = [0.5 0.0;
            0.0 0.5;
            0.0 0.25]

# The supercell dimensions required are based on (1/0.5, 1/0.5, 1/0.25).
supercell = get_supercell(q_points)
# Result: [2, 2, 4]
```
"""
function get_supercell!(supercell :: AbstractVector{I}, q_points :: AbstractMatrix{T}) where {T, I<:Integer}
    function cellvalue(x) :: I
        if abs(x) < 1e-10
            return 1
        end
        round(1 / abs(x))
    end
    maximum!(cellvalue, supercell, q_points)
end
function get_supercell!(supercell :: AbstractVector{I}, q_points :: AbstractMatrix{T}, cell :: AbstractMatrix{T}, reciprocal_vectors:: AbstractMatrix{T}; buffer=default_buffer()) where {T, I<:Integer}
    # Convert the q points in crystal coordinates
    @no_escape buffer begin
        q_points_fract = @alloc(T, size(q_points)...)
        cryst_cart_conv!(q_points_fract, q_points, cell, reciprocal_vectors, false; q_space=true)
        get_supercell!(supercell, q_points_fract)
        nothing
    end
end
function get_supercell(q_points :: AbstractMatrix{T}, args...; kwargs...) :: Vector{Int} where T
    ndims = size(q_points, 1)
    supercell = zeros(Int, ndims)
    get_supercell!(supercell, q_points, args...; kwargs...)
    supercell
end



