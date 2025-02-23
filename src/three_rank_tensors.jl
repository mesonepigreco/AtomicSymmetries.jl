"""
Deal with 3-rank tensors
"""


@doc raw"""
    symmetrize_tensor!(fc :: AbstractArray{T, 3}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer())

Apply the symmetrization on a 3-rank tensor.
Look to the documentation of `symmetrize_fc!` with a 2-rank tensor for more details.
"""
function symmetrize_tensor!(fc :: AbstractArray{T, 3}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer()) where T

    dimensions = get_dimensions(symmetry_group)

    @no_escape buffer begin
        fc_cryst = @alloc(T, dimensions, dimensions, dimensions)
        transf_matrix = @alloc(T, dimensions, dimensions)
        inverse_transf_matrix = @alloc(T, dimensions, dimensions)

        # Prepare the crystal and cartesian transformation matrices
        get_cryst_cart_transform_matrix!(transf_matrix, cell; buffer=buffer, cart_to_cryst=true)
        get_cryst_cart_transform_matrix!(inverse_transf_matrix, cell; buffer=buffer, cart_to_cryst=false)


        # Perform the symmetrization
        tmp_fc = @alloc(T, size(fc)...)
        work = @alloc(T, dimensions, dimensions, dimensions)
        tmp_fc .= 0.0
        for i_sym in 1:length(sym.symmetries)
            smat = symmmetry_group.symmetries[i_sym]
            irt = sym.irt[i_sym]

            for i in 1:n_atoms, j in 1:n_atoms, k in 1:n_atoms
                i_s = irt[i]
                j_s = irt[j]
                k_s = irt[k]

                # Convert to crystalline
                @views apply_transformation!(fc_cryst, fc[dimensions*(i_s - 1) + 1: dimensions*i,
                                             dimensions*(j_s - 1) + 1 : dimensions*j,
                                             dimensions*(k_s - 1) + 1 : dimensions*k],
                                             transf_matrix)



                # Apply symmetries
                apply_transformation!(work, fc_cryst, smat)
                
                # go back in cartesian (store in fc_cryst)
                apply_transformation!(fc_cryst, work, inverse_transf_matrix)

                # Add to the temporary symmetrized tensor
                tmp_fc[dimensions*(i-1)+1 : dimensions*i,
                       dimensions*(j-1)+1 : dimensions*j,
                       dimensions*(k-1)+1 : dimensions*k] .+= fc_cryst
            end
        end
        tmp_fc ./= length(sym.symmetries)

        nothing
    end
end

function apply_transformation!(dest :: AbstractArray{T,3}, tensor :: AbstractArray{T, 3}, transform_matrix :: AbstractMatrix)
    dim = size(cell, 1)
    @assert size(tensor, 1) == dim

    dest .= 0.0 
    # Apply the transform matrix
    # TODO: Speedup
    for a in 1:dim, b in 1:dim, c in 1:dim,
        α in 1:dim, β in 1:dim, γ in 1:dim
        dest[a, b, c] += tensor[α, β, γ] * transform_matrix[α, a] * transform_matrix[β, b] * transform_matrix[γ, c]
    end
end


function cart_cryst_tensor_conversion!(dest :: AbstractArray{T,3}, tensor :: AbstractArray{T, 3}, cell :: AbstractMatrix{T}; cart_to_cryst = true, buffer=default_buffer(), transform_matrix=nothing) where T
    dim = size(cell, 1)
    @assert size(tensor, 1) == dim

    if transform_matrix == nothing
        transform_matrix = zeros(T, dim, dim)
        get_cryst_cart_transform_matrix!(transform_matrix, cell; buffer=buffer, cart_to_cryst=cart_to_cryst)
    end

    apply_transformation!(dest, tensor, transform_matrix)
    ## Apply the transform matrix
    #for a in 1:dim, b in 1:dim, c in 1:dim,
    #    α in 1:dim, β in 1:dim, γ in 1:dim
    #    dest[a, b, c] += tensor[α, β, γ] * transform_matrix[α, a] * transform_matrix[β, b] * transform_matrix[γ, c]
    #end
end


@doc raw"""
    get_tensor_generator!(generator :: AbstractArray{T, 3}, generator_index :: Int,
        symmetry_group :: Symmetries{U}, cell :: AbstractMatrix; 
        permute_sym :: Bool = true, normalize :: Bool = true)


Get the $i$-th generator of a 3-rank tensors. 
The generator is an element of the basis of the
vector subspace of 3-rank tensors that are invariant under the 
symmetry group.
Projecting any 3-rank tensor in the basis of the generators 
results in performing the symmetrization.

This function operates in place.

## Paramenters

- `generator` : The 3-rank tensor that will be filled with the generator.
- `generator_index` : The index of the generator to get.
- `symmetry_group` : The symmetry group.
- `cell` : The cell of the system.

### Optional Parameters

- `permute_sym` : If true, also the index permutation invariance is applied.
- `normalize` : If true, the generator is normalized.
- `buffer` : The Bumper.jl stack to use.

"""
function get_tensor_generator!(generator :: AbstractArray{T, 3}, generator_index :: Int,
        symmetry_group :: Symmetries{U}, cell :: AbstractMatrix; 
        permute_sym :: Bool = true, normalize :: Bool = true, buffer=defualt_buffer()) where T

    # Set the generator to zero
    generator .= 0

    n_modes = size(generator, 1)

    # Get the indices of the generator
    k_index = (generator_index - 1) % n_modes + 1
    residual = (generator_index - 1) ÷ n_modes + 1
    j_index = (residual - 1) % n_modes + 1
    i_index = (residual - 1) ÷ n_modes + 1

    if permute_sym
        generator[i_index, j_index, k_index] = 1
        generator[i_index, k_index, j_index] = 1
        generator[k_index, i_index, j_index] = 1
        generator[k_index, j_index, i_index] = 1
        generator[j_index, i_index, k_index] = 1
        generator[j_index, k_index, i_index] = 1
    else
        generator[generator_index] = 1
    end

    symmetrize_tensor!(generator, cell, symmetry_group; buffer=buffer)


    if normalize
        normvalue = norm(generator)
        generator ./= normvalue
    end
end


@doc raw"""
    get_tensor_from_generators!(tensor :: AbstractArray{T, 3}, generators :: AbstractVector{Int},
        coefficients :: AbstractVector{T}, symmetry_group :: Symmetries{U}, cell :: AbstractMatrix;
        buffer=default_buffer())

Get the real space 3-rank tensor from its projection on the basis of the generators.
The function stores the result in the `tensor` argumnet in place.

The tensor is assumed in Cartesian coordinates.

## Parameters

- `tensor` : The 3-rank tensor to fill.
- `generators` : The indices of the generators.
- `coefficients` : The coefficients of the generators.
- `symmetry_group` : The symmetry group.
- `cell` : The cell of the system.

### Optional Parameters

- `buffer` : The Bumper.jl stack to use.
- `cache_allocation` : The tensor to use as cache. It is used to store the result of the symmetrization of the generators.
"""
function get_tensor_from_generators!(tensor :: AbstractArray{T, 3}, generators :: AbstractVector{Int},
        coefficients :: AbstractVector{T}, symmetry_group :: Symmetries{U}, cell :: AbstractMatrix;
        buffer=default_buffer(); cache_allocation :: AbstractArray{T, 3} = zeros(T, size(tensor)...)) where {T, U}
    
    tensor .= 0

    for i in 1:length(generators)
        get_tensor_generator!(cache_allocation, generators[i], symmetry_group, cell; buffer=buffer)
        @simd for k in 1:length(tensor)
            tensor[k] += coefficients[i] * generator[k]
        end
    end
end


@doc raw"""
    get_coefficients_from_tensor!(coefficients :: AbstractVector{T}, tensor :: AbstractArray{T, 3},
        generators :: AbstractVector{Int}, symmetry_group :: Symmetries{U}, cell :: AbstractMatrix;
        buffer=default_buffer())

Get the coefficients of the generators from the real space 3-rank tensor.
The function stores the result in the `coefficients` argumnet in place.

The tensor is assumed in Cartesian coordinates.

## Parameters

- `coefficients` : The coefficients of the generators to fill.
- `tensor` : The 3-rank tensor.
- `generators` : The indices of the generators.
- `symmetry_group` : The symmetry group.
- `cell` : The cell of the system.

### Optional Parameters

- `buffer` : The Bumper.jl stack to use.
- `cache_allocation` : The tensor to use as cache. It is used to store the result of the symmetrization of the generators.
"""
function get_coefficients_from_tensor!(coefficients :: AbstractVector{T}, tensor :: AbstractArray{T, 3},
        generators :: AbstractVector{Int}, symmetry_group :: Symmetries{U}, cell :: AbstractMatrix;
        buffer=default_buffer(), cache_allocation :: AbstractArray{T, 3} = zeros(T, size(tensor)...)) where {T, U}

    for i in 1:length(generators)
        get_tensor_generator!(cache_allocation, generators[i], symmetry_group, cell; buffer=buffer)
        coefficients[i] = reshape(cache_allocation, :)' * reshape(tensor, :)
    end
end


@doc raw"""
    get_3rank_tensor_generators(symmetry_group :: Symmetries{U}, cell :: AbstractMatrix{T};
        buffer=default_buffer()) :: Vector{Int} where {T, U}

Get the generators of the 3-rank tensors that are invariant under the symmetry group.
The generators represent a basis of the vector space of 3-rank tensors that are invariant under the symmetry group.

The generators are stored as a vector of integers. Each integer represents a generator
whose 3-rank tensor representation can be retrived using the function
`get_tensor_generator!`.

## Parameters

- `symmetry_group` : The symmetry group.
- `cell` : The cell of the system.

### Optional Parameters

- `buffer` : The Bumper.jl stack to use.

## Returns

- `generators :: Vector{Int}` : The vector of the indices of the generators.
"""
function get_3rank_tensor_generators(symmetry_group :: Symmetries{U}, cell :: AbstractMatrix{T};
        buffer=default_buffer()) :: Vector{Int} where {T, U}

    n_dims = get_dimensions(symmetry_group)
    n_atoms = get_n_atoms(symmetry_group)
    n_modes = n_dims * n_atoms

    generators = Vector{Int}()

    @no_escape buffer begin
        # TODO: Use sparse representation
        tensor_cache = @alloc(T, n_modes, n_modes, n_modes)
        tensor_cache .= 0

        # The following array keeps track of which elements have been already generated.
        not_necessary_indices = @alloc(Int, n_modes * n_modes * n_modes)
        not_necessary_indices .= 0

        for i in 1:n_modes*n_modes*n_modes
            # Convert the index to the indices of the tensor
            c = (i - 1) ÷ n_modes + 1
            residual = (i - 1) ÷ n_modes + 1
            b = (residual - 1) % n_modes + 1
            a = (residual - 1) ÷ n_modes + 1

            if c < b || b < a
                continue
            end

            # Skip if i ∈ not_necessary_indices
            if not_necessary_indices[i] == 1
                continue
            end

            # Prepare the current generator
            get_tensor_generator!(tensor_cache, i, symmetry_group, cell;
                                  buffer=buffer, normalize=false)

            # Check if the tensor has a nonzero entry
            normvalue = sum(abs2, reshape(tensor_cache, :))
            if normvalue > 1e-8
                # Add the tensor
                push!(generators, i)
            end

            # Add all the entry that needs to be skipped
            for j in i:n_modes*n_modes*n_modes
                # Convert the index to the indices of the tensor
                c2 = (i - 1) ÷ n_modes + 1
                residual = (i - 1) ÷ n_modes + 1
                b2 = (residual - 1) % n_modes + 1
                a2 = (residual - 1) ÷ n_modes + 1

                if c2 < b2 || b2 < a2
                    continue
                end

                # Check wether the tensor is nonzero
                if abs(tensor_cache[j]) > 1e-8
                    # Check if this item is not already in the list
                    not_necessary_indices[j] = 1
                end
            end
        end
    end

    return generators
end


@doc raw"""
    get_asr_violation!(violation_tensor ::Array{T, 3}, coefficients :: AbstractVector{T}, generators :: AbstractVector{Int}, symmetry_group :: Symmetries, cell :: AbstractMatrix{U}; buffer=simple_buffer())  where {T,U}

Compute the tensor quantifying the violation of the acoustic sum rule.
For a 3-rank tensor, the ASR violation is quantified by a tensor of size

``
\Delta_{ab}^{\alpha\beta\gamma} = \sum_c T_{abc}^{\alpha\beta\gamma}
``

where the $a$ runs over atoms while greek letters over Cartesian coordinates.

## Parameters


- `violation_tensor::Array{T, 3}`, the tensor $\Delta_{ab}^{\alpha\beta\gamma}$. The dimension of the tensor is 3, 3nat, 3nat
- `coefficients :: AbstractVector{T}`, the components on the generators
- `generators :: AbstractVector{Int}`, the generators of the symmetry.
- `symmetry_group :: Symmetries` the symmetry group
- `cell :: AbstractMatrix{U}`, the primitive cell

### Optional arguments

- `buffer=default_buffer()` the Bumper.jl stack.

"""
function get_asr_violation!(violation_tensor :: Array{T, 3}, coefficients :: AbstractVector{T},
        generators :: AbstractVector{Int}, symmetry_group :: Symmetries, 
        cell :: AbstractMatrix
    end
