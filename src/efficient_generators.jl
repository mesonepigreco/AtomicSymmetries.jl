@doc raw"""
    struct Generator{T, N} 
        index :: Int
        atom_index :: Vector
        matrix :: Array{T, N+1}
    end

An efficient memory representation of the generators.

The generators can be written so that they store only the ndim^rank matrix of the generator times the
number of symmetries in the primitive cell. 
This is an extremely memory efficient way to store the generators,
and it is also very easy to implement.
"""
struct Generator{T, N, M} 
    matrix :: Array{T, M}
    index :: Int
    atom_indices :: Vector{Int}
    cartesian_indices :: Vector{Int}
    transposed_indices :: Matrix{Int}
    dimension :: Int
    n_atoms :: Int
end

@doc raw"""
    Generator{T, M}(index :: Int, symmetry_group :: Symmetries) 

Construct the generator of rank M using the given index and symmetry group.
The type of the generator is identified by T
"""
function Generator{T, rank}(index :: Int, symmetry_group :: Symmetries; buffer=default_buffer()) 
    # Assert that M is an integer
    @assert rank isa Int

    dim = symmetry_group.dimension
    n_symmetries_primitive = get_n_symmetries_primitive_cell(symmetry_group)
    n_symmetries = get_nsymmetries(symmetry_group)
    n_atoms = get_n_atoms(symmetry_group)

    # Get the atom indices
    atom_indices = zeros(Int, rank)
    transposed_indices = zeros(Int, rank, n_symmetries)
    cartesian_indices = zeros(Int, rank)

    get_atomic_indices!(atom_indices, index, dim, symmetry_group.n_atoms)
    get_cartesian_indices!(cartesian_indices, index, dim, symmetry_group.n_atoms)

    matrix = zeros(T, [dim for i in 1:rank]..., n_symmetries_primitive)
    # Get all the transposed atomic indices
    @no_escape buffer begin
        irt_inv = @alloc(Int, n_atoms)

        first_generator = @alloc(T, [dim for i in 1:rank]...)
        first_generator .= 0.0
        first_generator[cartesian_indices...] = 1.0
        new_cartesian = @alloc(T, [dim for i in 1:rank]...)

        # Symmetrize the matrix if the atomic indices is the same
        for i in 1:rank
            for k in i+1:rank
                if atom_indices[i] == atom_indices[k]
                    # We need to symmetrize the matrix in the i-th and k-th dimension
                    # TODO: Symmetrize the matrix
                end
            end
        end

        for i in 1:n_symmetries

            # Invert the irt
            for k in 1:n_atoms
                irt_inv[symmetry_group.irt[i][k]] = k
            end

            for k in 1:rank
                transposed_indices[k, i] = irt_inv[atom_indices[k]]
            end
        end


        # Fill the matrix
        for i in 1:n_symmetries_primitive 
            apply_sym_fc!(view(matrix, (ntuple(_ -> Colon(), rank)..., i)),
        end

    end
end

@doc raw"""
    get_atomic_indices!(indices ::AbstractVector{Int}, index :: Int, dimension :: Int, nat :: Int)

Get the atomic indices of the given index in the given dimension.
This convert a linearized index to a multidimensional atomic index.
It only stores the atomic component of the index, not the cartesian one.
"""
function get_atomic_indices!(indices ::AbstractVector{Int}, index :: Int, dimension :: Int, nat :: Int)
    rank = length(indices)
    index -= 1
    for i in 1:rank
        indices[i] = index ÷ dimension + 1
        index = div(index, dimension * nat)
    end
end

@doc raw"""
    get_cartesian_indices!(indices ::AbstractVector{Int}, index :: Int, dimension :: Int, nat :: Int)

Get the cartesian indices
"""
function get_cartesian_indices!(indices ::AbstractVector{Int}, index :: Int, dimension :: Int, nat :: Int)
    rank = length(indices)
    index -= 1
    for i in 1:rank
        indices[i] = index % dimension + 1
        index = div(index, dimension * nat)
    end
end


function symmetrize_indices!(tensor :: Array{T, N}, i :: Int, j :: Int; buffer=default_buffer()) where {T, N}
    # Symmetrize the tensor in the i-th and j-th dimension
    @no_escape buffer begin

        for h in 1:N
            if h != i && h != j
                #TODO: SYmmetrize

            end
        end
    end
end
