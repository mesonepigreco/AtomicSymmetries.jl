# AtomicSymmetries.jl Documentation

A package to enforce symmetry constraints on atomic positions in a crystal structure.
It also works on dynamical matrices and force constants.

The package also exploits a memory efficient representation of the symmetry independent components of a tensor.
It can be used in pairs with `Spglib.jl` and the python version of spglib via PyCall for symmetry recognition.


## The symmetry group

The key object in this package is the `Symmetries` type, which represents a symmetry group.
The ``Symmetries`` is a container for the symmetry operations under which the specific crystal
is invariant.

```@docs
Symmetries
```

To build the symmetry group, one can create the following constructor

```julia 
symmetry_group = get_identity_symmetry_group(Float64)
```

This will create a symmetry group with the identity operation only.
The symmetry group can be extended by adding new symmetry operations with the function `add_symmetry!`.

```julia
add_symmetry!(symmetry_group, symmetry_operation; irt = nothing)
```

The ``symmetry_operation`` is a 3x3 matrix representing the rotation part of the symmetry operation.
The optional argument ``irt`` is a vector representing how the symmetry operation exchanges the atoms.
In particular the simmetry maps atom ``i`` to atom ``irt[i]``. If it is not provided, 
the symmetry operation is assumed not to exchange the atoms.

Once the fundamental symmetry operations have been added, the symmetry group can closed
by completing the irreducible representation. This is done by the function ``complete_symmetry_group!``.

```julia
complete_symmetry_group!(symmetry_group)
```

Once completed, the symmetry group can be exploited to enforce symmetry constraints on 
displacement vectors and 2-rank tensors as

```julia
symmetry_group.symmetrize_centroid!(vector)
symmetry_group.symmetrize_fc!(matrix)
```

### Build your own symmetry group (API)

To build a custom symmetry group, you can exploit the following subroutines

```@docs
complete_symmetry_group!
AtomicSymmetries.add_symmetry!
AtomicSymmetries.get_identity_symmetry_group
AtomicSymmetries.get_cylindrical_symmetry_group
AtomicSymmetries.get_spherical_symmetry_group
AtomicSymmetries.get_full_inversion_symmetry_group
```

### Structure symmetrization

Sometimes it is useful to symmetrize the atomic positions of a crystal structure.
To follow the correct Wyckoff positions. 
This is achieved with the function `symmetrize_positions!`.

```julia
symmetrize_positions!(cartesian_coords, cell, symmetry_group)
```

This function takes the atomic positions in cartesian coordinates, the cell matrix and the symmetry group.
The cell must be column-wise, i.e., each column is a primitive vector of the cell.

It is also possible just to apply the translational symmetries to either a vector or a force constants matrix.
This is achieved calling the subroutine

```julia
apply_translations!(vector, translations)
apply_translations!(matrix, translations)
```
where the `translations` is a vector of vector of Ints, each element represent a translations which maps each respective atomic index in the corresponding one. The `translations` object can be obtained from a symmetry group using the subroutine `get_translations`.

Here the APIs for these calls

```@docs
symmetrize_positions!
get_translations
apply_translations!
```

## The symmetry generators

From the symmetry group, we can obtain a vectorial subspace that is invariant under the symmetry operations.
A basis of this subspace is given by the symmetry generators.
For vectorial quantities, the generators can be obtained as

```julia
generators = get_vector_generators(symmetry_group)
```

For 2-rank tensors, the generators can be obtained as

```julia
generators = get_matrix_generators(symmetry_group)
```

The generators are vectors of indexes that can be used to build the symmetry independent components of a tensor.
This allows to store each generator as a 64-bit integer, which is more memory efficient than storing the full tensor.
The full vector/2-rank tensor can be retriven with the function ``get_vector_generator!``/``get_matrix_generator!``.

```julia
# Retrive the first element from the generators
i = 1
vector = zeros(3)
get_vector_generator!(vector, generators[i], symmetry_group)
```
And analogously for 2-rank tensors.


The generators can be used to project any vector or 2-rank tensor in the symmetry invariant subspace.

```julia
coefficients = zeros(length(generators))
my_vector = ...
get_coefficients_from_vector!(coefficients, my_vector, generators, symmetry_group)
```

The previous function projects the vector `my_vector` in the symmetry invariant subspace and stores the coefficients in the vector `coefficient`.
The coefficients can be used to reconstruct the original vector (symmetrized) as

```julia
final_vector = similar(my_vector)
get_centroids_from_generators!(final_vector, generators, coefficients, symmetry_group)
```

The same works for 2-rank tensors.
```julia
# Get the coefficients of the matrix projected in the symmetric subspace
coefficients = zeros(length(generators))
my_matrix = ...
get_coefficients_from_fc!(coefficients, my_matrix, generators, symmetry_group)

# And reconstruct bach the matrix from the coefficients
final_matrix = similar(my_matrix)
get_fc_from_generators!(final_matrix, generators, coefficients, symmetry_group)
```

### Generators (API)

```@docs
AtomicSymmetries.get_vector_generators
AtomicSymmetries.get_matrix_generators
AtomicSymmetries.get_vector_generator!
AtomicSymmetries.get_coefficients_from_vector!
AtomicSymmetries.get_centroids_from_generators!
AtomicSymmetries.get_coefficients_from_fc!
AtomicSymmetries.get_fc_from_generators!
```


## Spglib integration

The symmetry group can be directly constructed exploiting Spglib to recognize the symmetry operations
of a crystal structure.

```julia
get_symmetry_group_from_spglib(positions::AbstractMatrix{<: Real}, cell::AbstractMatrix{<:Real}, types::Vector{<:Int};  symprec::Float64 = 1e-6, type::Type = Float64, spglib_py_module = nothing) :: Symmetries
```

Here the arguments are the atomic positions (in crystal coordinates), the cell matrix and the atomic types.
Optionally, the symprec parameter can be used to set the tolerance for the symmetry recognition (passed to Spglib).

Since the `Spglib.jl` implementation is much less mature than the python version, 
if needed, it is possible to pass the module of the python version of spglib to the function to replace the Spglib.jl implementation
with the Python API from the official spglib package. This requires `PyCall.jl` to be installed.

### Spglib API

```@docs
get_symmetry_group_from_spglib
```

### How to get crystal coordinates

The code also allows for a quick conversion between cartesian and crystal coordinates.
Assuming `cartesian_positions` is a 3xNat matrix of atomic positions in cartesian coordinates and `cell` is the 3x3 cell matrix (each column is a primitive vector), the crystal coordinates can be obtained as

```julia
crystal_positions = similar(cartesian_positions)
get_crystal_coords!(crystal_positions, cartesian_positions, cell)
```
Optionally, a `Bumper.jl` buffer can be passed to the function to avoid memory allocations, otherwise, the `default_buffer()` is retrived.

The cartesian coordinates can be obtained as

```julia
get_cartesian_coords!(cartesian_positions, crystal_positions, cell)
```
This function does not require any memory allocation.

Notably, a much faster, nonallocating implementation can be used if the reciprocal vectors are available,
which can perform transformation of both real space and q space. Reciprocal vectors
can be computed using the `get_reciprocal_lattice!` subroutine, while the crystal-cartesian conversion
is obtained through `cryst_cart_conv!`.

```@docs
cryst_cart_conv!
get_reciprocal_lattice!
```


## Filter symmetries

It is possible to filter symmetries incompatible with a given external perturbation.
At this stage, only linear perturbations are supported.
For example, to filter the symmetries that are not compatible with a perturbation along the x direction, one can use

```julia
filter_invariant_symmetries!(symmetry_group, [1.0, 0.0, 0.0])
```

All symmetry operations not leaving the perturbation vector invariant are removed from the symmetry group.
Since version 0.2, it is possible to parse a vector of size n_dimension * n_atoms, 
with a different displacement vector acting on each atom.

```@docs
filter_invariant_symmetries!
```

## Symmetry sparsification

A new feature available since version 0.7 is the possibility to get the symmetry matrices as CSC sparse matrices.
To sparsify a symmetry group, you just need to use the `sparse` method from the SparseArrays library (novel dependency of 0.7)

```julia
Using SparseArrays
sparse_symmetry_group = sparse(symmetry_group)

# Apply the symmetry 5th operation on a vector v (previosly defined)
v_new = apply_sparse_symmetry(sparse_symmetry_group.symmetries[5], v)
```

Notably, this is differentiable via Zygote and Enzyme, so it allows to implement symmetrization in a differentiable way.
v could also be a series of vector as a Matrix where each vector is stored as a column of v.

```@docs
apply_sparse_symmetry
```
