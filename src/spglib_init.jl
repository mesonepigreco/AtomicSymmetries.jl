# This files init the symmetries from the spglib library
# Initialize a symmetry group starting from spglib

@doc raw"""
    get_symmetry_group_from_spglib(positions :: AbstractMatrix, cell :: AbstractMatrix, types :: Array{Int}; 
        symprec :: Float64 = 1e-6,
        type :: Type = Float64, spglib_py_module = nothing) :: Symmetries

Build a symmetry group from the spglib library.
Optionally, this function can be called with a spglib python module.
In this way, the python module will be used to get the symmetry operations (since the julia spglib module is buggy).

# Arguments

- `positions::AbstractMatrix`: The atomic positions in the cell (crystallographic coordinates), with shape `(3, N)`.
- `cell::AbstractMatrix`: The cell matrix with shape `(3, 3)`.
- `types::Array{Int}`: The atomic types.

Optional arguments:
- `symprec::Float64`: The symmetry precision.
- `type::Type`: The numerical precision type for the symmetry operations.
- `spglib_py_module`: The spglib python module. If not provided, the default julia spglib module is used.

Alternatively, you can pass an ASE Atoms object.
"""
function get_symmetry_group_from_spglib(positions::AbstractMatrix{<: Real}, cell::AbstractMatrix{<:Real}, types::Vector{<:Int};  symprec::Float64 = 1e-6, type::Type = Float64, spglib_py_module = nothing) :: Symmetries
    nat = size(positions, 2)
    ndim = size(positions, 1)

    newpos = collect(eachcol(positions))
    newcell = collect(eachcol(cell))

    @assert ndim == 3 "Only 3D systems are supported for SPGLIB symmetry groups."

    # Build the SPGLIB cell object
    if spglib_py_module != nothing
        @debug "Using the python spglib module"
        spglib_cell = (cell', positions', types)
        symmetry_py_dict = spglib_py_module.get_symmetry(spglib_cell, symprec=symprec)
        n_sym = size(symmetry_py_dict["rotations"], 1)
        Rmat = zeros(type, 3, 3, n_sym)
        Tmat = zeros(type, 3, n_sym)
        R = []
        T = []
        for i in 1:n_sym
            Rmat[:, :, i] .= symmetry_py_dict["rotations"][i, :, :]'
            Tmat[:, i] .= symmetry_py_dict["translations"][i, :]
            push!(R, Rmat[:, :, i])
            push!(T, Tmat[:, i])
        end
    else
        @debug "Building the SPGLIB cell object"
        spglib_cell = Cell(newcell, newpos, types)
        
        # Get the symmetry operations
        @debug "Getting the symmetry operations"
        R, T = Spglib.get_symmetry(spglib_cell, symprec)
    end

    # Create the new symmetry group
    @debug "Creating the symmetry group"
    sym_group = get_identity_symmetry_group(type; dims=ndim, n_atoms=nat, 
                                            translations = true)

    # Add the symmetry operations excluding the first identity
    @debug "Adding the symmetry operations"
    for i in 2:length(R)
        irt = zeros(Int, nat)
        get_irt!(irt, positions, R[i], T[i])

        # Add the symmetry operation
        add_symmetry!(sym_group, R[i], update=false; irt=irt)
        push!(sym_group.translations, T[i])
    end

    round_translations!(sym_group)

    @debug "Updating the symmetry functions"
    update_symmetry_functions!(sym_group)


    return sym_group
end


@doc raw"""
    round_translation!(translation :: AbstractVector; max_value :: Int = 8)

round the translations to the exact values.
Fractional translation can be 1/2 1/3 2/3 1/4 ... 
The maximum division is given by max_value parameter.
"""
function round_translation!(translation :: AbstractVector{T}; max_value :: Int=10, threshold = 1e-8) where T
    # Round the translation
    for i in 1:length(translation)
        thr = 1000.0
        good_j = 1
        for j in 2:max_value
            if abs(translation[i] * j - round(translation[i] * j)) < thr
                thr = abs(translation[i] * j - round(translation[i] * j))
                good_j = j
            end

            if thr < threshold
                break
            end
        end
        translation[i] = round(translation[i] * good_j) / good_j
    end
end

@doc raw"""
    round_translations!(symmetry_group:: Symmetries; max_value :: Int = 10)

This subroutine corrects the translations to assume exact fractional values.

Fractional translation can be 1/2 1/3 2/3 1/4 ... 
The maximum division allowed is given by max_value parameter.
"""
function round_translations!(symmetry_group :: Symmetries; max_value :: Int = 10, kwargs...)
    for i in 1:length(symmetry_group.translations)
        round_translation!(symmetry_group.translations[i]; max_value = max_value, kwargs...)
    end
end

