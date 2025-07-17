# Fourier symmetrization

The symmetrization can be performed also in Fourier space.
This is implemented now for force-constant dynamical matrices and vectors (displacements, forces, ...).

A vector is transformed from real to q-space with the following convention:

$$
v_k(\vec q) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{-i 2\pi \vec R\cdot \vec q} v_k(\vec R)
$$

Note the sign of the Fourier and the normalization prefactor.

## Fourier transform

The API to perform the fourier transform occur mainly with `vector_r2q!`, `vector_q2r!` which, respectively, trasform a vector from real to q space and vice-versa. Transformation of matrices occur with `matrix_r2q!`, `matrix_q2r!`. All these operations are inplace. The matrices are assumed in crystal coordinates, but in this case it should not matter.

The detailed API calls are

```@docs
vector_r2q!
vector_q2r!
matrix_r2q!
matrix_q2r!
```


## Applications of symmetries in Fourier space

The application of symmetries in Fourier space must also account how points in q space are mapped by the symmetry operations.

For this, the important information about how q points are related by symmetries
needs to be computed and stored.
This identification is performed by the helper function `get_irt_q!`, which identifies, for a given symmetry operation, the i->j mapping between q points. Q points mapped into themselves by the same set of symmetry operations form the socalled small-group of ``q``, while the set of ``q`` points mapped by all the symmetries of a crystal is called the star of the ``q`` point.
Due to time-inversion symmetry, the dynamical matrix must also satisfy the condition

$$
D(q) = D^\dagger(-q + G)
$$

therefore it is necessary also to keep track, for each q point, which one is the corresponding ``-q + G`` in the mesh. This mapping is computed by the helper function `get_minus_q!`. All these information needs to be stored when applying symmetries. Therefore we defined a new Symmetries struct that ihnerits from the `GenericSymmetries` called `SymmetriesQSpace`

```@docs
SymmetriesQSpace
AtomicSymmetries.get_irt_q!
AtomicSymmetries.get_minus_q!
```

### Application of symmetries

Applying a symmetry means transforming a vector or a matrix (already in q-space) into a new vector (matrix). If the vector (matrix) is invariant under that transformation, then that transformation belong to the symmetry group.

To apply the symmetry to a matrix we use the `apply_symmetry_vectorq!`. For the matrix, we use `apply_symmetry_matrixq!`. Both these function modify in-place the first argument, storing the result of the transformation there.
Note that, since symmetries are stored in crystalline components, both the vector and the matrix must be in crystalline components. 

```@docs
apply_symmetry_vectorq!
apply_symmetry_matrixq!
```

### Enforcing symmetries

One of the most useful operation to do is enforce a specific matrix or vector in q-space to satisfy a given symmetry group.

This can be implemented by applying the complete irreducible representation of the symmetry group.



 
