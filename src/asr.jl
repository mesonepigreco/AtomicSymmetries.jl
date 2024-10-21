@doc raw"""
ASRconstraint!

Apply the ASR constraint to a rank-`rank` tensor of dimension `dimension`.
The application works as 

```julia
my_asr! = ASRconstraint!(dimension, 2)
my_tensor = rand(nat * dimension, nat * dimension)
my_asr!(my_tensor)
```

This will apply the ASR constraint to the tensor `my_tensor` in place.
"""
struct ASRVectorConstraint!
    dimension::Int
end
function (asr::ASRVectorConstraint!)(tensor::AbstractVector)
    # Get the dimension of the tensor
    dimension = asr.dimension

    # Get the number of atoms
    nat = length(tensor) ÷ dimension

    mytensor = reshape(tensor, dimension, nat)

    # Apply the ASR constraint
    for i in 1:nat
        @views mytensor[:, i] .-= sum(mytensor[:, i]) / dimension
    end
end

# Apply the ASR constraint to a rank-2 tensor
struct ASRMatrixConstraint!
    dimension::Int
end
function (asr::ASRMatrixConstraint!)(matrix :: AbstractMatrix)
    nat = size(matrix, 1) ÷ asr.dimension
    proj = zeros(eltype(matrix), asr.dimension * nat)
    proj_view = reshape(proj, asr.dimension, nat)
    copy_mat = copy(matrix)
    for t in 1:asr.dimension
        # Get the projector
        proj .= 0
        @views proj_view[t, :] .= 1.0 / √(asr.dimension * nat)

        # subtract the projector
        matrix_projected = proj' * matrix * proj
        matrix .-= matrix_projected * (proj * proj')
    end
end



