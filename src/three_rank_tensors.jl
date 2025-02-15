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
        
    # Apply the transform matrix
    for a in 1:dim, b in 1:dim, c in 1:dim,
        α in 1:dim, β in 1:dim, γ in 1:dim
        dest[a, b, c] += tensor[α, β, γ] * transform_matrix[α, a] * transform_matrix[β, b] * transform_matrix[γ, c]
    end
end

