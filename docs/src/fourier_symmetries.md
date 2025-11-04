# Symmetries in Fourier space

From version 0.8, `AtomicSymmetries.jl` provided the possibility to apply symmetries
to vector and matrices directly in Fourier space.
This is implemented now for force-constant dynamical matrices and vectors (displacements, forces, ...).

A vector is transformed from real to q-space with the following convention:

``
\displaystyle
\tilde v_k(\vec q) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{-i 2\pi \vec R\cdot \vec q} v_k(\vec R)
``

``
\displaystyle
v_k(\vec R) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{i 2\pi \vec R\cdot \vec q} \tilde v_k(\vec q)
``


Note the sign of the Fourier and the normalization prefactor. 
This convention allows for correctly transforming the matrices, however, it introduces a size inconsistency on the vectors.
If we have a periodic vector in the cell, its ``q`` fourier transformed counterpart will be ``\sqrt {N_q}`` times
the value in the primitive cell. So be carefull when extracting ``\Gamma`` point data from periodic vectors.

With this convention, we recover the standard rule for the matrices.

``
\displaystyle
\tilde \Phi_{ab}(\vec q) = \sum_{\vec R} e^{2\pi i \vec q\cdot \vec R}\Phi_{a;b + \vec R}
``

``
\displaystyle
\Phi_{ab} = \frac{1}{N_q} \sum_{\vec q}
\tilde\Phi_{ab}(\vec  q) e^{2i\pi \vec q\cdot[\vec R(a) - \vec R(b)]}
``

Note that these transformation of matrices and vector are consistent so that matrices and vector written as outer product can be consistently transformed

``
\displaystyle
\Phi(\vec R) = \sum_i\sum_{\vec R} \vec v_i(\vec R_1) \otimes \vec v_i(\vec R_1 + \vec R)
``

``
\displaystyle
\tilde \Phi(\vec q) = \sum_i \vec {\tilde v}_i(\vec q) \otimes \vec {\tilde v_i}(-\vec q)
``


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

``
D(q) = D^\dagger(-q + G)
``

therefore it is necessary also to keep track, for each q point, which one is the corresponding ``-q + G`` in the mesh. This mapping is computed by the helper function `get_minus_q!`. All these information needs to be stored when applying symmetries. Therefore we defined a new Symmetries struct that ihnerits from the `GenericSymmetries` called `SymmetriesQSpace`. Note that, to initialize the symmetries in q-space, we **must** use the symmetries object (`Symmetries`) evaluated in the primitive cell. The correct initialization of symmetries could be checked with the subroutine `check_symmetries`, which will spot if a different cell has been employed when initializing the symmetries.

```@docs
SymmetriesQSpace
AtomicSymmetries.get_irt_q!
AtomicSymmetries.get_minus_q!
AtomicSymmetries.check_symmetries
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

This can be implemented by applying the complete irreducible representation of the symmetry group. Symmetrization of an ent `\phi` is obtained as

``
\Phi = \frac{1}{N}\sum_{i=1}^N S_i(\phi)
``

where ``S_i`` is the symmetry operation. The two functions performing the symmetrization are `symmetrize_matrix_q!` and `symmetrize_vector_q!`. Also in this case, the dynamical matrix must be provided in crystalline coordinates.

To symmetrize vector and matrices already provided in cartesian coordinates,
we must use the appropriate subroutines `symmetrize_vector_cartesian_q!` and
`symmetrize_matrix_cartesian_q!`. 
These two subroutines correctly convert the vector/matrix in crystal coordinates 
before applying the symmetries, and then convert the symmetrized result back in cartesian space.
They are the most used subroutines to perform symmetrization in q-space,
the equivalent of `symmetrize_vector!` and `symmetrize_fc!` for real space.


Here the complete API

```@docs
symmetrize_vector_q!
symmetrize_matrix_q!
symmetrize_vector_cartesian_q!
symmetrize_matrix_cartesian_q!
```


## Manipulating q points

The fourier transform depends on the knowledge of few vectors:
`q_points`, `itau`, and `R_lat` (evenutally `translations`, for inverse 
fourier transform into a matrix).

All these properties can be evaluated from the core source.
For example, to obtain the supercell to which the q points are commensurate, 
we can use the `get_supercell` method as

```@docs
get_supercell!
```

Analogously, we can get the translations `R_lat` as

```@docs
get_R_lat!
```
