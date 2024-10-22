@doc raw"""
    get_vector_generators(symmetry_group::Symmetries{U};
                   func_apply_constraints! = nothing)
                    :: Vector{Int} where {T, U}

Get the generators of the symmetry group for the vectors.

# Arguments
- `symmetry_group::Symmetries{U}`: The symmetry group to be considered in the generator creation process.
- `func_apply_constraints!::Function` (optional): A function to apply constraints to the parameters. Defaults to `nothing`.

# Returns
- `Vector{Int}`: A vector of indices representing the set of independent generators.

# Description
This function generates a set of independent generators for a given symmetry group.
These generators proved a basis for the space of vectors that are invariant under the symmetry group.

The independence of a generator is determined by its norm and its linear independence from previously accepted generators. If a generator is found to be linearly dependent but not identical to a previous one, the function throws an error indicating that this scenario is not yet implemented.

The function ultimately returns a vector of indices representing the independent generators found during the process.

# Examples
```julia
symmetry_group = # Symmetries object
generators = get_generators(symmetry_group)
"""
function get_vector_generators( 
    symmetry_group::Symmetries{U};
    func_apply_constraints! = nothing,
    type = Float64) :: Vector{Int} where {U}
    n_dims = get_dimensions(symmetry_group)
    n_atoms = get_n_atoms(symmetry_group)
    n_modes = n_dims * n_atoms
    
    # Generate an (empty) vector of generators
    generators = Vector{Int}()
    generator = zeros(type, n_dims * n_atoms)

    old_vectors = Vector{Vector{type}}()

    # TODO: get the baseline generator
    for i in 1:n_modes
        get_vector_generator!(generator, i, symmetry_group;
            func_apply_constraints! = func_apply_constraints!,
            normalize=false)

        # Add the generator to the acceptable ones
        normvalue = norm(generator)
        if normvalue > 1e-8
            generator ./= normvalue

            # Check if the generator is independent from the others
            independent = true
            for j in 1: length(old_vectors)
                if norm(generator - old_vectors[j]) < 1e-8
                    independent = false
                    break
                end

                if sqrt(generator' * old_vectors[j]) > 1e-4
                    println("Current $i generator: $(generator)")
                    println("Original norm:", normvalue)
                    println("Old generators:")
                    println(old_vectors)
                    println("Linearly dependent with $j : $(sqrt(generator' * old_vectors[j]))")
                    throw("The generators are linearly dependent but not the same. This is not implemented yet.")
                end
            end

            if independent
                push!(generators, i)
                push!(old_vectors, copy(generator))
            end
        end
    end

    println("Found $(length(generators)) generators: $(generators)")

    return generators
end
function get_matrix_generators( 
    symmetry_group::Symmetries{U};
    func_apply_constraints! = nothing,
    type = Float64) :: Vector{Int} where {U}
    n_dims = get_dimensions(symmetry_group)
    n_atoms = get_n_atoms(symmetry_group)
    n_modes = n_dims * n_atoms
    
    # Generate an (empty) vector of generators
    generators = Vector{Int}()
    generator = zeros(type, n_dims * n_atoms, n_dims * n_atoms)
    old_vectors = Vector{Matrix{type}}()

    # TODO: get the baseline generator
    for i in 1:n_modes*n_modes
        # skip the symmetric part
        i_index = (i-1) % n_modes + 1
        j_index = (i-1) ÷ n_modes + 1
        # if i_index < j_index
        #     continue
        # end

        get_matrix_generator!(generator, i, symmetry_group;
            func_apply_constraints! = func_apply_constraints!,
            normalize=false)

        @debug "processing $i generator"

        # Add the generator to the acceptable ones
        normvalue = norm(generator)
        if normvalue > 1e-8
            generator ./= normvalue
            @debug "Found a generator in $i"

            # Check if the generator is independent from the others
            independent = true
            for j in 1: length(old_vectors)
                if norm(generator - old_vectors[j]) < 1e-8
                    independent = false
                    break
                end

                if sqrt(reshape(generator, :)' * reshape(old_vectors[j], :)) > 1e-4
                    println("Current $i generator: $(generator)")
                    println("Original norm:", normvalue)
                    println("Old generators:")
                    println(old_vectors)
                    println("Linearly dependent with $j : $(sqrt(generator' * old_vectors[j]))")
                    throw("The generators are linearly dependent but not the same. This is not implemented yet.")
                end
            end

            if independent
                push!(generators, i)
                push!(old_vectors, copy(generator))
            end
        end
    end

    println("Found $(length(generators)) generators: $(generators)")

    return generators
end



# @doc raw"""
#     get_generator!(generator::Vector{T}, generator_index::Int, symmetry_group :: Symmetries{T};
#         use_sqrt_representation=true,
#         optimize_struct=true,
#         optimize_nltransf=true) :: Vector{Int}
# 
# get the generator vector from the index.
# """
@doc raw"""
    get_generator!(generator::Vector{T}, generator_index::Int, n_modes::Int, n_layers::Int, symmetry_group::Symmetries{U};
                   use_sqrt_representation=true, optimize_struct=true, optimize_nltransf=true, 
                   func_apply_constraints!=nothing, baseline_generator=nothing, normalize=true) where {T, U}

Modify `generator` in-place to represent a specific generator of a transformation, subject to given constraints and symmetries.

The generator_index parameter ranges from 1 to N_max, where N_max represents the maximum dimension of the parameters, and it should be noted that some indices within this range may yield identical generators due to the underlying symmetries or constraints in the system.

# Arguments
- `generator::Vector{T}`: The generator vector to be modified.
- `generator_index::Int`: Index specifying which generator to construct.
- `n_modes::Int`: The number of modes in the system.
- `n_layers::Int`: The number of layers in the neural network model.
- `symmetry_group::Symmetries{U}`: The symmetry group to be imposed on the generator.
- `use_sqrt_representation::Bool` (optional): Flag to use the square root representation. Defaults to `true`.
- `optimize_struct::Bool` (optional): Flag to optimize the structure. Defaults to `true`.
- `optimize_nltransf::Bool` (optional): Flag to optimize nonlinear transformations. Defaults to `true`.
- `func_apply_constraints!::Function` (optional): A function to manually apply constraints to the parameters. Defaults to `nothing`.
- `baseline_generator::Vector{T}` (optional): A baseline generator for comparison. Defaults to `nothing`.
- `normalize::Bool` (optional): Flag to normalize the generator. Defaults to `true`.

# Description
This function constructs a generator vector that represents a transformation in a specified manner. It initializes the `generator` vector with zeros and sets the `generator_index` element to 1. The function builds a `scha` (Structured Component Histogram Analysis) and an `nltransf` (Non-Linear Transform) based on the specified parameters. The `scha` and `nltransf` are then used to set the parameters of the generator.

If `func_apply_constraints!` is provided, it is used to apply constraints to the `scha` and `nltransf`. If `baseline_generator` is provided, the function adjusts the generator relative to this baseline. Symmetry constraints from `symmetry_group` are imposed on the `nltransf` and `scha`.

Finally, the function updates the `generator` with the parameters obtained from `scha` and `nltransf`, optionally subtracting the baseline generator and normalizing the result.

# Examples
```julia
generator = zeros(Float64, 10)
get_generator!(generator, 2, 5, 3, my_symmetry_group)
```
"""
function get_vector_generator!(generator :: Vector{T}, generator_index:: Int, symmetry_group :: Symmetries{U};
    func_apply_constraints! =nothing,
    baseline_generator = nothing,
    normalize=true) where {T, U}

    # Prepare the starting vector
    generator .= 0
    generator[generator_index] = 1

    # Constrain manually the parameters indexed in constrained_parameters
    if func_apply_constraints! != nothing
        func_apply_constraints!(generator)
        if baseline_generator === nothing
            baseline_generator = zeros(T, length(generator))
            func_apply_constraints!(baseline_generator)
        end
    end

    # Apply the symmetry group
    symmetry_group.symmetrize_centroid!(generator)

    # Subtract the baseline generator
    if baseline_generator != nothing
        generator .-= baseline_generator
    end
    
    # Normalize the generator
    if normalize
        normvalue = norm(generator)
        @assert normvalue > 1e-8 "The norm of the $i-th generator is zero"
        generator ./= normvalue    
    end
end
function get_matrix_generator!(generator :: AbstractMatrix{T}, generator_index:: Int, symmetry_group :: Symmetries{U};
    func_apply_constraints! =nothing,
    baseline_generator = nothing,
    normalize=true) where {T, U}

    # Prepare the starting vector
    generator .= 0
    generator[generator_index] = 1

    # Constrain manually the parameters indexed in constrained_parameters
    if func_apply_constraints! != nothing
        func_apply_constraints!(generator)
        if baseline_generator === nothing
            baseline_generator = similar(generator)
            baseline_generator .= 0
            func_apply_constraints!(baseline_generator)
        end
    end

    # Apply the symmetry group
    symmetry_group.symmetrize_fc!(generator)

    # Subtract the baseline generator
    if baseline_generator != nothing
        generator .-= baseline_generator
    end
    
    # Normalize the generator
    if normalize
        normvalue = norm(generator)
        @assert normvalue > 1e-8 "The norm of the $i-th generator is zero"
        generator ./= normvalue    
    end
end




@doc raw"""
    get_params_from_generators!(params :: Vector{T}, generators::Vector{Int}, coefficients :: Vector{T})

Return the parameters from the generators and the coefficients.
"""
function get_centroids_from_generators!(centroids:: AbstractVector{T}, generators::AbstractVector{Int}, coefficients :: AbstractVector{T},
    symmetries :: Symmetries{U},
    n_modes::Int;
    func_apply_constraints! =nothing,
    ) where {T, U}

    centroids .= 0
    generator = similar(centroids)

    for i in 1:length(generators)
        get_vector_generator!(generator, generators[i], n_modes, symmetries;
                              func_apply_constraints! = func_apply_constraints!)
        
        centroids .+= coefficients[i] * generator
    end

    if func_apply_constraints! != nothing
        func_apply_constraints!(params)
    end
end
function get_fc_from_generators!(fc :: AbstractMatrix{T}, generators::AbstractVector{Int}, coefficients :: AbstractVector{T},
    symmetries :: Symmetries{U};
    func_apply_constraints! =nothing,
    ) where {T, U}

    fc .= 0
    generator = similar(fc)

    for i in 1:length(generators)
        get_matrix_generator!(generator, generators[i], symmetries;
                              func_apply_constraints! = func_apply_constraints!)
        
        fc .+= coefficients[i] * generator
    end

    if func_apply_constraints! != nothing
        func_apply_constraints!(fc)
    end
end

@doc raw"""
    get_coefficients_from_params!(coefficients :: Vector{T}, params :: Vector{T}, generators :: Vector{Int},
        n_modes::Int, n_layers :: Int;
        use_sqrt_representation=true,
        optimize_struct=true,
        optimize_nltransf=true) where {T}
"""
function get_coefficients_from_vector!(coefficients :: AbstractVector{T}, vector:: AbstractVector{T}, 
    generators :: AbstractVector{Int},
    symmetries :: Symmetries{U},
    func_apply_constraints! =nothing) where {T, U}

    generator = similar(vector)
    for i in 1: length(generators)
        get_vector_generator!(generator, generators[i], 
            symmetries;
            func_apply_constraints! = func_apply_constraints!)
        
        coefficients[i] = generator' * vector
    end
end
function get_coefficients_from_fc!(coefficients :: AbstractVector{T}, fc :: AbstractMatrix{T}, 
    generators :: AbstractVector{Int},
    symmetries :: Symmetries{U},
    func_apply_constraints! =nothing) where {T, U}

    generator = similar(fc)
    for i in 1: length(generators)
        get_matrix_generator!(generator, generators[i], 
            symmetries;
            func_apply_constraints! = func_apply_constraints!)
        
        coefficients[i] = reshape(generator, :)' * reshape(fc, :) 
    end
end


