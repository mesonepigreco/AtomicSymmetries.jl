# Symmetries in Fourier space

From version 0.8, `AtomicSymmetries.jl` provided the possibility to apply symmetries
to vector and matrices directly in Fourier space.
This is implemented now for force-constant dynamical matrices and vectors (displacements, forces, ...).

A vector is transformed from real to q-space with the following convention:

```math
\tilde v_k(\vec q) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{-i 2\pi \vec R\cdot \vec q} v_k(\vec R)
```

```math
v_k(\vec R) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{i 2\pi \vec R\cdot \vec q} \tilde v_k(\vec q)
```


Note the sign of the Fourier and the normalization prefactor. 
This convention allows for correctly transforming the matrices, however, it introduces a size inconsistency on the vectors.
If we have a periodic vector in the cell, its ``q`` fourier transformed counterpart will be ``\sqrt {N_q}`` times
the value in the primitive cell. So be carefull when extracting ``\Gamma`` point data from periodic vectors.

With this convention, we recover the standard rule for the matrices.

```math
\tilde \Phi_{ab}(\vec q) = \sum_{\vec R} e^{2\pi i \vec q\cdot \vec R}\Phi_{a;b + \vec R}
```

```math
\Phi_{ab} = \frac{1}{N_q} \sum_{\vec q}
\tilde\Phi_{ab}(\vec  q) e^{2i\pi \vec q\cdot[\vec R(a) - \vec R(b)]}
```

Note that these transformation of matrices and vector are consistent so that matrices and vector written as outer product can be consistently transformed

```math
\Phi(\vec R) = \sum_i\sum_{\vec R} \vec v_i(\vec R_1) \otimes \vec v_i(\vec R_1 + \vec R)
```

```math
\tilde \Phi(\vec q) = \sum_i \vec {\tilde v}_i(\vec q) \otimes \vec {\tilde v_i}(-\vec q)
```

Notably, this convention introduces two main properties that must be handled with care.
The ``\Gamma`` value of the fourier transform is not the average over the supercell of the same
quantity. If you want to obtain the average, you must divide by ``\sqrt {N_q}`` (the number of q-points).
If the `R_lat` is not centered around zero, and the coordinates passed as `v_sc` are absolute values of positions,
then the ``\Gamma`` value of the fourier transform will be shifted by a total translation which is the average of the translations of the supercell lattice
vectors.
This can be avoided by either removing the corner of the supercell from the positions before performing the fourier transform, by centering R_lat around 0,
or by removing this translational average *a posteriori* using the method `shift_position_origin!`.



## Fourier transform

The API to perform the fourier transform occur mainly with `vector_r2q!`, `vector_q2r!` which, respectively, trasform a vector from real to q space and vice-versa. Transformation of matrices occur with `matrix_r2q!`, `matrix_q2r!`. All these operations are inplace. The matrices are assumed in crystal coordinates, but in this case it should not matter.

To shift the origin for the fourier transformed absolute positions, use the method `shift_position_origin!` as

The detailed API calls are

```@docs
vector_r2q!
vector_q2r!
matrix_r2q!
matrix_q2r!
shift_position_origin!
```


## Symmetries in Q space

The application of symmetries in Fourier space must also account how points in q space are mapped by the symmetry operations.

For this, the important information about how q points are related by symmetries
needs to be computed and stored.
This identification is performed by the helper function `get_irt_q!`, which identifies, for a given symmetry operation, the i->j mapping between q points. Q points mapped into themselves by the same set of symmetry operations form the socalled small-group of ``q``, while the set of ``q`` points mapped by all the symmetries of a crystal is called the star of the ``q`` point.
Due to time-inversion symmetry, the dynamical matrix must also satisfy the condition

``
D(q) = D^\dagger(-q + G)
``

therefore it is necessary also to keep track, for each q point, which one is the corresponding ``-q + G`` in the mesh. This mapping is computed by the helper function `get_minus_q!`. All these information needs to be stored when applying symmetries. Therefore we defined a new Symmetries struct that ihnerits from the `GenericSymmetries` called `SymmetriesQSpace`. Note that, to initialize the symmetries in q-space, we **must** use the symmetries object (`Symmetries`) evaluated in the primitive cell. The correct initialization of symmetries could be checked with the subroutine `check_symmetries`, which will spot if a different cell has been employed when initializing the symmetries.

Since the q points must be passed in crystal coordinates, it may be useful to get the reciprocal lattice, which can be done with ``get_reciprocal_lattice!`` (see section on crystal coordinates for the API)

### Application of symmetries

Applying a symmetry means transforming a vector or a matrix (already in q-space) into a new vector (matrix). If the vector (matrix) is invariant under that transformation, then that transformation belong to the symmetry group.

Notably, the symmetries in the supercell are always the symmetries in the primitive cell times all possible translations operated by the lattice vectors compatible with the chosen supercell. On the contrary, the symmetries in q space are always only equal to the symmetries in the primitive cell.
The reason is that translations are automatically incorporated in the q space representation by the block theorem:

```math
D(q, q') = D(q)\delta(q - q')
```

This means that applying each symmetry operation in ``q`` space is equivalent to averaging the result of the same symmetry operation in the supercell averaging among all possible translations.

The application of a symmetry in q space can be performed by considering how the force-constant matrix transform in real space under a symmetry operation ``S``.

```math
S[\tilde\Phi_{ab}(\bm q)] = \sum_{\bm R} e^{2\pi i \bm q\cdot \bm R}S^\dagger\Phi_{S\bm a, S(\bm b + \bm R)}S
```
The transformation also changes which atoms are considered. However, we must be careful with the convention adopted for the Fourier transform. In fact, we have
that $\bm a$ and $\bm b$ are the positions on the atom in the primitive cell considered. The vectors $S\bm a$ and $S(\bm b + \bm R)$ may not correspond to atoms in the primitive cell, but rather folded in the supercell. 
To solve this issue, we need to define, for each symmetry operation, which atom in the primitive cell is mapped into which other atom in the primitive cell. This is indicated with ``s(a)`` and ``s(b)``. Also, we need to consider
what is the translation vector ``\bm t_{s,a}`` that brings the vector ``S\bm a`` inside the primitive cell. With this information, we can rewrite the transformation as
```math
\bm t_{s,a} = S\bm a - \bm{s(a)}
```

```math
S[\tilde\Phi_{ab}(\bm q)] = \sum_{\bm R} e^{2\pi i \bm q\cdot \bm R}S^\dagger\Phi_{s(a) + \bm t_{s,a}, s(b) + \bm t_{s, b} + S\bm R}S
```

Exploiting the translational invariance, we can remove the $\bm t_{s,a}$ vector from the first index of the supercell force constant matrix, and rewrite the expression as


```math
S[\tilde\Phi_{ab}(\bm q)] = \sum_{\bm R} e^{2\pi i \bm q\cdot \bm R}S^\dagger\Phi_{s(a), s(b) + \bm t_{s, b} - \bm t_{s, a} + S\bm R}S
```
By defining ``\bm R' = \bm t_{s, b} - \bm t_{s, a} + S\bm R``, we can rewrite the summation in ``\bm R'`` as


```math
S[\tilde\Phi_{ab}(\bm q)] = \sum_{\bm R'} e^{2\pi i \bm q\cdot S^{-1} (\bm R' + \bm t_{s,a} - \bm t_{s,b})}S^\dagger\Phi_{s(a), s(b) + \bm R'}S
```
Since we work in crystal coordinates and reciprocal vectors, ``S^{-1}\neq S^\dagger``. Therefore, we have

```math
S[\tilde\Phi_{ab}(\bm q)] = \sum_{\bm R'} e^{2\pi i [(\bm S^{-1})^\dagger\bm q]\cdot(\bm R' + \bm t_{s,a} - \bm t_{s,b})}S^\dagger\Phi_{s(a), s(b) + \bm R'}S
```

Which is equivalent to the Fourier transform of the dynamical matrix at the transformed q-point ``(\bm S^{-1})^\dagger\bm q``, times a phase factor.
This is how symmetries operates in q space:
```math
\bm S_\text{recip} = \left(\bm S^{-1}\right)^\dagger
```

```math
S[\tilde\Phi_{ab}(\bm q)] = S^\dagger \tilde\Phi_{s(a)s(b)}(S_\text{recip}\bm q) S e^{2\pi i (S_\text{recip}\bm q)\cdot ( \bm t_{s,a} - \bm t_{s,b})}
```

Note that the ``S_\text{recip}q`` vector in the phase factor and in the dynamical matrix can be always folded back into the first Brilluin zone. In fact the dynamical matrix is periodic in the reciprocal vector, while the phase factor is multiplied by a direct lattice vector. Thus, by adding a reciprocal lattice vector ``\bm G`` to ``S_\text{recip}\bm q``, the phase factor is multiplied by ``e^{2\pi i \bm G\cdot ( \bm t_{s,a} - \bm t_{s,b})}``, which is always equal to 1.

This transformation for each q point is operated by the subroutine `apply_symmetry_matrixq!`. Both these function modify in-place the first argument, storing the result of the transformation there. 
Note that, since symmetries are stored in crystalline components, both the vector and the matrix must be in crystalline components. This makes it also important that the ``\bm q`` points are provided in crystalline coordinates, to correctly compute the phase factor and the transformed ``S\bm q``.

```@docs
SymmetriesQSpace
apply_symmetry_vectorq!
apply_symmetry_matrixq!
AtomicSymmetries.get_irt_q!
AtomicSymmetries.get_minus_q!
AtomicSymmetries.check_symmetries
```

## Enforcing symmetries

One of the most useful operation to do is enforce a specific matrix or vector in q-space to satisfy a given symmetry group.

This can be implemented by applying the complete irreducible representation of the symmetry group. Symmetrization of an ent `\Phi` is obtained as

``
\Phi = \frac{1}{N}\sum_{i=1}^N S_i(\Phi)
``

where ``S_i`` is the symmetry operation. The two functions performing the symmetrization are `symmetrize_matrix_q!` and `symmetrize_vector_q!`. Also in this case, the dynamical matrix must be provided in crystalline coordinates.

To symmetrize vector and matrices already provided in cartesian coordinates,
we must use the appropriate subroutines `symmetrize_vector_cartesian_q!` and
`symmetrize_matrix_cartesian_q!`. 
These two subroutines correctly convert the vector/matrix in crystal coordinates 
before applying the symmetries, and then convert the symmetrized result back in cartesian space.
They are the most used subroutines to perform symmetrization in q-space,
the equivalent of `symmetrize_vector!` and `symmetrize_fc!` for real space.

The Hermitianity is not automatically imposed by the symmetrization procedure.
This allows to symmetrize matrices that are not necessarily hermitian, for example, the cross correlation matrices between different quantities.
Hermitianity and time-reversal symmetry can be imposed with the subroutine `impose_hermitianity_q!`, which enforces the condition.
The time-reversal symmetry corresponds to the condition that the original matrix in real space is real-valued.


Here the complete API

```@docs
symmetrize_vector_q!
symmetrize_matrix_q!
symmetrize_vector_cartesian_q!
symmetrize_matrix_cartesian_q!
impose_hermitianity_q!
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


