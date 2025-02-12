@doc raw"""
    mutable struct Symmetries{T}

The structure containing the symmetries of the system.

Once the symmetries have been initialized, 
the symmetrize_fc! and symmetrize_centroid! functions can be used to symmetrize the force constant matrix and the centroid.

The exchange_symmetry is a vector of length n_particles, 
where each element identify the id of the particle.
If two ids are equal, the particles are indistinguishable.

irt[i][j] is the index of the atom that is equivalent to the j-th atom before the symmetry is applied.
The expression is
$$
v_{\text{irt[i]}} = S v_i
$$

The name irt stands for "index of the representative of the transformation".
and it is in line with the notation used in the Quantum Espresso and the CellConstructor codes.

"""
mutable struct Symmetries{T}
    symmetries :: Vector{Matrix{T}}
    dimension :: Int
    n_particles :: Int
    exchange_symmetry :: Union{Vector{Int}, Nothing}
    irt :: Union{Vector{Vector{Int}}, Nothing}
    symmetrize_fc! :: Union{Function, Nothing}
    symmetrize_centroid! :: Union{Function, Nothing}
    enforce_noninteracting :: Vector{Int}
    translations :: Vector{Vector{T}}
end
get_nsymmetries(sym :: Symmetries) = length(sym.symmetries)
get_dimensions(sym :: Symmetries) = sym.dimension
get_n_atoms(sym :: Symmetries) = sym.n_particles
Base.length(sym :: Symmetries) = get_nsymmetries(sym)


# Override the Base.isempty function to check if the Symmetries object is empty
function Base.isempty(sym :: Symmetries{T}) where {T} 
    return isempty(sym.symmetries)
end


@doc raw"""
    update_symmetry_functions!(sym :: Symmetries{T}) :: where {T}

Initialize the symmetries of the system by setting the
symmetrize_fc! and symmetrize_centroid! functions.

Note that the symmetries must be set before calling this function.
"""
function update_symmetry_functions!(sym :: Symmetries{T}) where {T}
    function sym_fc!(fc :: AbstractMatrix{U}; buffer=default_buffer()) where {U}
        @no_escape buffer begin
            my_fc = @alloc(U, size(fc)...)
            # my_fc = similar(fc)
            my_fc .= 0.0

            # Apply the symmetries
            for i in 1:length(sym.symmetries)
                symmat = sym.symmetries[i]
                irt =sym.irt[i]
                apply_sym_fc!(my_fc, fc, symmat, sym.dimension, irt; buffer=buffer)
            end
            fc .= my_fc 
            fc ./= length(sym.symmetries)

            nothing #Avoid returning
        end

        # Apply the exchange symmetry
        if !isnothing(sym.exchange_symmetry)
            apply_exchange_symmetry!(fc, sym.exchange_symmetry, sym.dimension)
        end

        # Check if to force noninteracting objects 
        if length(sym.enforce_noninteracting) > 0
            enforce_noninteracting!(fc, sym.enforce_noninteracting, sym.dimension)
        end
    end

    function sym_centroid!(centroid :: AbstractVector{U}; apply_translations=false) where {U}
        my_centroid = similar(centroid)
        my_centroid .= 0.0

        # Check if translations are correctly initialized
        if apply_translations 
            if length(sym.translations) != length(sym.symmetries)
                error("The number of translations must be equal to the number of symmetries")
            end
        end

        for i in 1:get_nsymmetries(sym)
            symmat = sym.symmetries[i]
            irt = sym.irt[i]
            apply_sym_centroid!(my_centroid, centroid, symmat, sym.dimension, irt)
        end
        centroid .= my_centroid ./ length(sym.symmetries)

        # Apply the exchange symmetry
        if !isnothing(sym.exchange_symmetry)
            apply_exchange_symmetry!(centroid, sym.exchange_symmetry, sym.dimension)
        end
    end

    # Update the symmetrization functions
    sym.symmetrize_fc! = sym_fc!
    sym.symmetrize_centroid! = sym_centroid!

    # Update irt (index of the representative of the transformation)
    # if sym.n_particles > 0
    #    update_irt!(sym)
    # end
end    


@doc raw"""
    symmetrize_fc!(fc :: AbstractMatrix{T}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer()) where {T}

Symmetrize the force constant matrix `fc` using the symmetry group `symmetry_group`.
This function assumes that the force constant matrix is in cartesian coordinates,
opposite to `symmetry_group.symmetrize_fc!`, which assumes that the force constant matrix is in crystal coordinates.
Note that to symmetrize cartesian coordinates, also the primitive cell is required (`cell`).

The function operates in place, meaning that
the final result overwrites the input force constant matrix `fc`.

This function exploits Bumper.jl stack allocation
to avoid memory allocation.
The stack can be manually specified as an optional keyword argument `buffer`.
"""
function symmetrize_fc!(fc :: AbstractMatrix{T}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer()) where {T}
    if isnothing(symmetry_group.symmetrize_fc!)
        error("Symmetry group not initialized")
    end

    @no_escape buffer begin
        fc_cryst = @alloc(T, size(fc)...)

        cart_cryst_matrix_conversion!(fc_cryst, fc, cell; cart_to_cryst = true, buffer=buffer)
        symmetry_group.symmetrize_fc!(fc_cryst; buffer=buffer)
        cart_cryst_matrix_conversion!(fc, fc_cryst, cell; cart_to_cryst = false, buffer=buffer)

        nothing
    end
end

@doc raw"""
    symmetrize_vector!(vector :: AbstractVector{T}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer()) where {T}

Symmetrize the vector `vector` using the symmetry group `symmetry_group`.
The vector has a length of `dim * nat`, where `nat` is the number of atoms in the system.
It is assumed to represent a quantity that is invariant under translations (e.g. a force or a displacement).
If you want to symmetrize a quantity that is not invariant under translations (e.g. atomic positions),
use `symmetrize_positions!` instead.

This function assumes that the `vector` is provided in cartesian coordinates,
opposite to `symmetry_group.symmetrize_centroid!`, which assumes that the `vector` is in crystal coordinates.
Note that to symmetrize cartesian coordinates, also the primitive cell is required (`cell`).

The function operates in place, meaning that
the final result overwrites the input `vector`.

This function exploits Bumper.jl stack allocation
to avoid memory allocation.
The stack can be manually specified as an optional keyword argument `buffer`.
"""
function symmetrize_vector!(vector :: AbstractVector{T}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer()) where T
    if isnothing(symmetry_group.symmetrize_centroid!)
        error("Symmetry group not initialized")
    end

    ndims = symmetry_group.dimension

    @no_escape buffer begin
        vector_cryst = @alloc(T, length(vector))

        get_crystal_coords!(reshape(vector_cryst, ndims, :), reshape(vector, ndims, :), cell; buffer=buffer)
        symmetry_group.symmetrize_centroid!(vector_cryst)
        get_cartesian_coords!(reshape(vector, ndims, :), reshape(vector_cryst, ndims, :), cell)

        nothing
    end
end



@doc raw"""
    symmetrize_positions!(positions :: AbstractMatrix{T}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer()) where {T}

symmetrize an atomic coordinates in real space. 
This subroutie symmetrizes a system with Cartesian coordinats (positions)
using the specified symmetry group (that must include translations).

If you want to symmetrize a quantity that is invariant under translations (e.g. a force or a displacement),
use `symmetrize_vector!` instead.

The function operates in place, meaning that
the final result overwrites the input positions.

This function exploits Bumper.jl stack allocation
to avoid memory allocation, you can manually specify the stack buffer as an optional keyword argument `buffer`.
"""
function symmetrize_positions!(positions :: AbstractMatrix{T}, cell :: AbstractMatrix, symmetry_group :: Symmetries; buffer=default_buffer()) where {T}
    # Check the presence of translations
    if length(symmetry_group.translations) != length(symmetry_group.symmetries)
        error("The number of translations must be equal to the number of symmetries to symmetrize the positions")
    end

    # Convert to crystal coordinates
    @no_escape buffer begin 
        crystal_coords = @alloc(T, size(positions)...)
        get_crystal_coords!(crystal_coords, positions, cell; buffer=buffer)

        tmp_centroids = @alloc(T, length(crystal_coords))
        tmp_results = @alloc(T, length(crystal_coords))
        tmp_results .= 0.0
        
        # Loop over the symmetries
        for i in 1:length(symmetry_group.symmetries)
            tmp_centroids .= 0.0
            apply_sym_centroid!(tmp_centroids, reshape(crystal_coords, :),
                                symmetry_group.symmetries[i],
                                symmetry_group.dimension,
                                symmetry_group.irt[i];
                                translation = symmetry_group.translations[i])

            # Compute the δ from the centroid
            for j in 1:length(tmp_results)
                delta = tmp_centroids[j] - crystal_coords[j]
                delta -= round(delta)
                if abs(delta)>0.2
                    println()
                    println("Found a value of δ = $delta; Sym $i, coord $j")
                    println("translation = ", symmetry_group.translations[i])
                    iat = (j-1) ÷ 3 + 1
                    i_transf = findfirst(isequal(iat), symmetry_group.irt[i])
                    println("original atom = ", iat)
                    println("corresponding atom = ", i_transf)
                    println("original = ", crystal_coords[3*(iat-1)+1:3*iat])
                    println("corresponding original = ", crystal_coords[3*(i_transf-1)+1:3*i_transf])

                    println("symmetry = ", symmetry_group.symmetries[i])
                    println("transformed = ", tmp_centroids[3*(iat-1)+1:3*iat])
                end
                tmp_results[j] += delta
            end
        end
        tmp_results ./= length(symmetry_group.symmetries)

        # Add the δ from the original position
        for i in 1:length(crystal_coords)
            crystal_coords[i] += tmp_results[i]
        end

        get_cartesian_coords!(positions, crystal_coords, cell)
        nothing # <-- avoid returning
    end
end

@doc raw"""
    update_irt!(sym :: Symmetries{T}; coordinates :: Union{Nothing, Vector{T}} = nothing) where {T}

Update the irt (index of the representative of the transformation) attribute of the Symmetries object.
IRT is a vector of vectors, where irt[i][j] is the index of the atom that is equivalent to the j-th atom before the symmetry is applied.

If the coordinates are not provided, all symmetries are assumed not to change the position of the atoms.
"""
function update_irt!(sym :: Symmetries{U}; coordinates :: Union{Nothing, Vector{T}} = nothing) where {U, T}
    # Check if the number of particles is set
    @assert sym.n_particles > 0 "The number of particles is not set"

    if isnothing(coordinates)
        sym.irt = [collect(1:sym.n_particles) for i in 1:length(sym.symmetries)]
    else
        # Raise not implemented error
        error("Using a symmetry that change the atomic position is not implemented yet")
    end
end


@doc raw"""
    add_symmetry!(sym :: Symmetries{T}, symm :: Matrix{T}; update :: Bool = true, check_existing :: Bool = false) :: where {T}

Add a symmetry to the system.

If update is true, the symmetrize_fc! and symmetrize_centroid! functions are updated.
If check_existing is true, the symmetry is only added if it is not already in the list of symmetries.
"""
function add_symmetry!(sym :: Symmetries{T}, symm :: AbstractMatrix{U}; 
        update :: Bool = true, check_existing :: Bool = false, irt = nothing) where {T, U}
    if check_existing
        for s in sym.symmetries
            if isapprox(s, symm)
                return
            end
        end
    end

    if !isnothing(irt)
        push!(sym.irt, irt)
    end

    push!(sym.symmetries, symm)
    if update
        update_symmetry_functions!(sym)
    end
end


@doc raw"""
    apply_sym_fc!(result :: Matrix{T}, fc :: Matrix{T}, sym :: Matrix{T}, tions :: Int, irt :: Vector{Int}) where {T}

Apply the symmetry ``sym`` to the force constant matrix ``fc`` and add the result to ``result``.

irt indicates how the symmetry maps particle into each other.
"""
function apply_sym_fc!(result :: AbstractMatrix{T}, fc :: AbstractMatrix{T}, sym :: AbstractMatrix{U}, dimensions :: Int, irt :: AbstractVector{Int}; buffer=default_buffer()) where {T, U}
    # Use mul! to avoid allocating memory
    n_atoms = size(fc, 1) ÷ dimensions

    @no_escape buffer begin
        work = @alloc(T, dimensions, dimensions)
        #work = zeros(T, (dimensions, dimensions))
        work .= 0
        for i ∈ 1:n_atoms
            i_s = irt[i]
            for j in 1:n_atoms 
                j_s = irt[j]
                @views mul!(work, fc[dimensions*(i_s-1) + 1: dimensions*i_s, dimensions*(j_s-1)+1 : dimensions*j_s], sym, 1.0, 0.0)
                @views mul!(result[dimensions*(i-1) + 1: dimensions*i, dimensions*(j - 1) + 1: dimensions*j], 
                    sym', work, 1.0, 1.0)
            end
        end
        nothing
    end
end



@doc raw"""
    apply_sym_centroid!(result :: Vector{T}, centroid :: Vector{T}, sym :: Matrix{T}, dimensions :: Int, irt :: Vector{Int};
    translation = nothing) where {T}

Apply the symmetry ``sym`` to the centroid vector ``centroid`` and add the result to ``result``.

irt indicates how the symmetry maps particle into each other.
translation, if present, is a vector of which all atoms are translated after the symmetry operation.
"""
function apply_sym_centroid!(result :: AbstractVector{T}, centroid :: AbstractVector{T}, sym :: AbstractMatrix{U}, dimensions :: Int, irt :: AbstractVector{Int};
    translation = nothing) where {T,U}
    # Use mul! to avoid allocating memory
    n_atoms = length(centroid) ÷ dimensions
    if translation == nothing
        for i ∈ 1:n_atoms 
            j = irt[i]
            @views mul!(result[dimensions*(j-1) + 1: dimensions*j], sym, centroid[dimensions*(i-1) + 1: dimensions*i], 1.0, 1.0)    
        end
    else
        partial_result = zeros(T, dimensions)
        for i ∈ 1:n_atoms 
            j = irt[i]
            @views mul!(partial_result, sym, centroid[dimensions*(i-1) + 1: dimensions*i], 1.0, 0.0)    
            
            # Convert the partial result into the primitive cell
            # @views to_primitive_cell_cryst!(partial_result, [0.0, 0.0, 0.0])

            result[dimensions*(j-1) + 1: dimensions*j] .+= partial_result .+ translation
        end
    end
end


@doc raw"""
    apply_exchange_symmetry!(fc :: Matrix{T}, exchange_symmetry :: Vector{Int}, dims :: Int =3) where {T}
    apply_exchange_symmetry!(centroid :: Vector{T}, exchange_symmetry :: Vector{Int}, dims :: Int =3) where {T}

Apply the exchange symmetry to the force constant matrix or centroid vector.

TODO: use the set_parameters! function to avoid allocating memory for exchange_symmetry
"""
function apply_exchange_symmetry!(fc :: AbstractMatrix{T}, exchange_symmetry :: AbstractVector{Int}, dims :: Int) where {T}
    fc_dim = size(fc, 1)
    work1 = zeros(T, (dims, dims))
    work2 = zeros(T, (dims, dims))

    nat = length(exchange_symmetry)
    @assert fc_dim == dims * nat "The force constant matrix and exchange_symmetry do not match ($fc_dim != $dims * $nat; $nat atoms, $dims dimensions)"

    done_ids = Int[]
    count = 0
    for i in 1:length(exchange_symmetry)
        if exchange_symmetry[i] in done_ids
            continue
        end

        indistinguishable_ids = findall(x -> exchange_symmetry[x] == exchange_symmetry[i], 
            1:nat)

        # Create the diagonal block
        work1 .= 0
        for k in indistinguishable_ids
            work1 .+= view(fc, dims*(k-1) + 1: dims*k, dims*(k-1) + 1: dims*k)
        end
        work1 ./= length(indistinguishable_ids)

        # Create the off-diagonal block
        work2 .= 0
        count = 0
        for k in indistinguishable_ids
            for l in indistinguishable_ids
                if k > l
                    work2 .+= view(fc, dims*(k-1) + 1: dims*k, dims*(l-1) + 1: dims*l)
                    count += 1
                end
            end
        end
        work2 ./= count

        # Fill the matrix
        for k in indistinguishable_ids
            fc[dims*(k-1) + 1: dims*k, dims*(k-1) + 1: dims*k] .= work1
            for l in indistinguishable_ids
                if k > l
                    fc[dims*(k-1) + 1: dims*k, dims*(l-1) + 1: dims*l] .= work2
                    fc[dims*(l-1) + 1: dims*l, dims*(k-1) + 1: dims*k] .= work2'
                end
            end
        end 

        push!(done_ids, exchange_symmetry[i])
    end
end
function apply_exchange_symmetry!(centroid :: AbstractVector{T}, exchange_symmetry :: AbstractVector{Int}, dims :: Int) where {T}
    centroid_dim = length(centroid)
    work1 = zeros(T, dims)

    nat = length(exchange_symmetry)
    @assert centroid_dim == dims * nat "The centroid vector and exchange_symmetry do not match"

    done_ids = Int[]
    for i in 1:length(exchange_symmetry)
        if exchange_symmetry[i] in done_ids
            continue
        end

        indistinguishable_ids = findall(x -> exchange_symmetry[x] == exchange_symmetry[i], 
            1:nat)

        # Create the diagonal block
        work1 .= 0
        for k in indistinguishable_ids
            work1 .+= view(centroid, dims*(k-1) + 1: dims*k)
        end
        work1 ./= length(indistinguishable_ids)

        # Fill the matrix
        for k in indistinguishable_ids
            centroid[dims*(k-1) + 1: dims*k] .= work1
        end 

        push!(done_ids, exchange_symmetry[i])
    end
end


@doc raw""" 
    enforce_noninteracting!(fc :: Matrix{T}, list_of_particles :: Int, dimension :: Int)

Enforce that the force constant matrix is zero for the off-diagonal elements among particles in ``list_of_particles``.
This creates a noncorrelated wavefunction between the particles in the list.
""" 
function enforce_noninteracting!(fc :: AbstractMatrix{T}, list_of_particles :: AbstractVector{Int}, dimension :: Int) where {T}
    for i in list_of_particles
        for j in list_of_particles
            if i != j
                fc[dimension*(i-1) + 1: dimension*i, dimension*(j-1) + 1: dimension*j] .= 0
            end
        end
    end
end



@doc raw"""
    is_initialized(sym :: Symmetries{T}) :: Bool where {T}

Check if the symmetries have been initialized.
"""
function is_initialized(sym :: Symmetries{T}) :: Bool where {T}
    if isnothing(sym.symmetrize_fc!) 
        return false 
    end

    if isnothing(sym.symmetrize_centroid!) 
        return false
    end

    if isnothing(sym.irt)
        return false
    end

    if sym.dimension == 0
        return false
    end

    return true
end


function get_empty_symmetry_group(T :: Type) :: Symmetries{T}
    return Symmetries{T}([], 0, 0, nothing, [], nothing, nothing, [], [])
end


@doc raw"""
    complete_symmetry_group!(symmetries :: Symmetries{T})

Complete the symmetry group by adding the inverse symmetry operations
and the compositions of the symmetry operations, until the group is closed.

Since symmetry are unitary transformations, 
the inverse symmetry operation is the transpose of the symmetry operation.
"""
function complete_symmetry_group!(symmetries :: Symmetries{T}) where {T}
    # Add all possible inverse symmetry operations
    for i ∈ 1:length(symmetries.symmetries)
        symmat = symmetries.symmetries[i]

        # Check if the inverse symmetry operation is already present
        if !any(x -> isapprox(x, symmat'), symmetries.symmetries)
            # Add the inverse symmetry operation
            tmpsym = similar(symmat)
            tmpsym .= symmat'
            add_symmetry!(symmetries, tmpsym; update = false)
        end
    end

    # Look for compositions of the symmetry operations
    index_i = 1
    index_j = 1
    composed_sym = zeros(T, size(symmetries.symmetries[1]))
    # The while loops are used to continue the search even if new symmetries are added
    while index_i <= length(symmetries.symmetries)
        index_j = index_i
        while index_j <= length(symmetries.symmetries)
            @views sym_i = symmetries.symmetries[index_i]
            @views sym_j = symmetries.symmetries[index_j]
            composed_sym .= sym_i * sym_j

            #println("sym_$index_i = $sym_i , sym_$index_j = $sym_j ; total_sym = $composed_sym")

            # Check if the composition of the symmetry operations is already present
            if !any(x -> isapprox(x, composed_sym), symmetries.symmetries)
                #println("Adding")
                # Add the composition of the symmetry operations (new allocation is necessary)
                add_symmetry!(symmetries, copy(composed_sym); update = false)
            end

            index_j += 1
        end
        index_i += 1
    end

    # Update the symmetry functions
    update_symmetry_functions!(symmetries)
end


@doc raw"""
    get_full_inversion_symmetry_group(T :: Type, dims :: Int =2) :: Symmetries{T}


Return the group of inversion symmetries along the x, y, and z axis.
This suppose that atoms are not exchanged by the symmetry.
"""
function get_full_inversion_symmetry_group(T :: Type, dims :: Int =2; n_atoms=1) :: Symmetries{T}
    syms = get_empty_symmetry_group(T)
    syms.n_particles = n_atoms
    syms.dimension = dims

    # Add the identity matrix
    add_symmetry!(syms, Matrix{T}(I, dims, dims); update = false)

    # Add the inversion operations
    for i in 1:dims
        id_mat = Matrix{T}(I, dims, dims)
        id_mat[i,i] = -1
        add_symmetry!(syms, id_mat)
    end

    complete_symmetry_group!(syms)

    return syms
end


@doc raw"""
    get_spherical_symmetry_group(T :: Type, dims :: Int = 3) :: Symmetries{T}

Return the symmetry group of all rotation and inversion operations in
the dimension of the system.
"""
function get_spherical_symmetry_group(T :: Type, dims :: Int = 3) :: Symmetries{T}
    symmetries = get_empty_symmetry_group(T)

    # Set the dimension of the system
    symmetries.dimension = dims

    # Add the identity matrix
    add_symmetry!(symmetries, Matrix{T}(I, dims, dims); update = false)

    # Add the inversion operation
    add_symmetry!(symmetries, -Matrix{T}(I, dims, dims); update = false)

    # Add the rotation operations
    for i ∈ 1:dims
        for j ∈ 1:dims
            if i != j
                # Create the rotation matrix
                rotmat = Matrix{T}(I, dims, dims)
                rotmat[i, i] = 0.0
                rotmat[j, j] = 0.0
                rotmat[i, j] = 1.0
                rotmat[j, i] = -1.0

                # Add the rotation matrix
                add_symmetry!(symmetries, rotmat; update = false, check_existing = true)
            end
        end
    end

    # Complete the symmetry group
    complete_symmetry_group!(symmetries)

    return symmetries
end


@doc raw"""
    get_cylindrical_symmetry_group(T :: Type, axis :: Int) :: Symmetries{T}

The cylindrical symmetry group is the group of all rotation around one axis, plus all inversions.
"""
function get_cylindrical_symmetry_group(T :: Type, axis :: Int)
    symm_group = get_full_inversion_symmetry_group(T, 3)
    
    # Add a rotation around the chosen axis
    rotmat = Matrix{T}(I, 3, 3)
    i = mod(axis, 3) + 1
    j = mod(axis + 1, 3) + 1
    rotmat[i, i] = 0.0
    rotmat[j, j] = 0.0
    rotmat[i, j] = 1.0
    rotmat[j, i] = -1.0
    add_symmetry!(symm_group, rotmat; update = false, check_existing = true)

    # Complete the symmetry group
    complete_symmetry_group!(symm_group)
    return symm_group
end

@doc raw"""
    get_identity_symmetry_group(T :: Type) :: Symmetries{T}


This function returns the symmetry group containing only the identity matrix.
"""
function get_identity_symmetry_group(T :: Type; 
        dims :: Int = 2,
        n_atoms :: Int = 1,
        translations :: Bool = false) :: Symmetries{T}
    syms = get_empty_symmetry_group(T)
    syms.n_particles = n_atoms
    syms.dimension = dims

    # Add the identity matrix
    add_symmetry!(syms, Matrix{T}(I, dims, dims); 
                  update = false, 
                  irt = collect(1:n_atoms))


    complete_symmetry_group!(syms)

    if translations 
        push!(syms.translations, zeros(T, dims))
    end

    return syms
end



@doc raw"""
    get_irt!(irt, coords, matrix, translation)

Get the irt index.

# Arguments

- `irt::Array{Int}`: The irt index (inplace modified)
- `coords::AbstractMatrix`: The atomic positions in the cell (crystallographic coordinates), with shape `(3, N)`.
- `matrix::AbstractMatrix`: The rotation matrix.
- `translation::AbstractVector`: The translation vector.
"""
function get_irt!(irt, coords, matrix, translation)
    new_coords = matrix * coords

    nat = size(coords, 2)
    ndims = size(coords, 1)

    # debugvalue = sum(translation.^2) > 1e-8
    debugvalue = false
    dist = 0

    for i in 1:nat 
        min_dist = 1000.0
        min_j = 1
        for j in 1:nat
            nval = 0
            for k in 1:ndims
                dist = new_coords[k, i] - coords[k, j] + translation[k]
                nval += (dist - round(dist))^2
            end

            if debugvalue
                println("i: $i, j: $j, distance: $nval, notransl: $(norm(new_coords[:, i] - coords[:, j])); min_dist: $min_dist")

            end

            if nval < min_dist
                min_j = j
                min_dist = nval
            end

            if min_dist < 1e-5
                break
            end
        end
        if min_dist > 0.1
            println("The distance between the atoms is too large: $min_dist")
            println("Atom $i: ", coords[:, i])
            println("Atom $min_j: ", coords[:, min_j])
            error("Error while initializing the symmetry group.")
        # else
        #     println("IRT[$min_j] = $i")
        #     println("min_dist = $min_dist")
        #     println("original coords = ", coords[:, min_j])
        #     println("transformed coords = ", new_coords[:, i])
        #     println("matrix = ", matrix)
        #     println()
        end

        #irt[min_j] = i
        irt[i] = min_j
    end
end
