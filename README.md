# AtomicSymmetries.jl

A package to enforce symmetry constraints on atomic positions in a crystal structure.
It also works on dynamical matrices and force constants.

The package also exploits a memory efficient representation of the symmetry independent components of a tensor.
It can be used in pairs with Spglib.jl and the python version of spglib via PyCall for symmetry recognition.


## The symmetry group

The key object in this package is the `Symmetries` type, which represents a symmetry group.
The ``Symmetries`` is a container for the symmetry operations under which the specific crystal
is invariant.

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

The previous function projects the vector ``my_vector`` in the symmetry invariant subspace and stores the coefficients in the vector ``coefficient``.
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


## Spglib integration

The symmetry group can be directly constructed exploiting Spglib to recognize the symmetry operations
of a crystal structure.

```julia
get_symmetry_group_from_spglib(positions::AbstractMatrix{<: Real}, cell::AbstractMatrix{<:Real}, types::Vector{<:Int};  symprec::Float64 = 1e-6, type::Type = Float64, spglib_py_module = nothing) :: Symmetries
```

Here the arguments are the atomic positions (in crystal coordinates), the cell matrix and the atomic types.
Optionally, the symprec parameter can be used to set the tolerance for the symmetry recognition (passed to Spglib).

Since the ``Spglib.jl`` implementation is much less mature than the python version, 
if needed, it is possible to pass the module of the python version of spglib to the function to replace the Spglib.jl implementation
with the Python API from the official spglib package. This requires ``PyCall`` to be installed.

### How to get crystal coordinates

The code also allows for a quick conversion between cartesian and crystal coordinates.
Assuming ``cartesian_positions`` is a 3xNat matrix of atomic positions in cartesian coordinates and ``cell`` is the 3x3 cell matrix (each column is a primitive vector), the crystal coordinates can be obtained as

```julia
crystal_positions = similar(cartesian_positions)
get_crystal_coords!(crystal_positions, cartesian_positions, cell)
```
Optionally, a ``Bumper.jl`` buffer can be passed to the function to avoid memory allocations, otherwise, the ``default_buffer()`` is retrived.

The cartesian coordinates can be obtained as

```julia
get_cartesian_coords!(cartesian_positions, crystal_positions, cell)
```
This function does not require any memory allocation.


## Filter symmetries

It is possible to filter symmetries that are not compatible with a given external perturbation.
At this stage, only linear perturbations are supported.
For example, to filter the symmetries that are not compatible with a perturbation along the x direction, one can use

```julia
filter_invariant_symmetries!(symmetry_group, [1.0, 0.0, 0.0])
```

All symmetry operations that do not leave invariant the perturbation vector are removed from the symmetry group.
