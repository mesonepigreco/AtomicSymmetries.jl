# Symmetries in Fourier space

From version 0.8, `AtomicSymmetries.jl` provided the possibility to apply symmetries
to vector and matrices directly in Fourier space.
This is implemented now for force-constant dynamical matrices and vectors (displacements, forces, ...).

A vector is transformed from real to q-space with the following convention:

```math
\tilde v_a(\vec q) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{-i 2\pi \vec q\cdot (\vec R + \vec\tau_a)} v_a(\vec R)
```

```math
v_a(\vec R) = \frac{1}{\sqrt{N_q}} \sum_{\vec q} e^{i 2\pi \vec q\cdot(\vec R + \vec\tau_a)} \tilde v_a(\vec q)
```

where ``\vec\tau_a`` is the position of the atom ``a`` inside the primitive cell,
so that ``\vec R + \vec\tau_a`` is the equilibrium position of the atom in the supercell.

The phase factor employs the atomic positions ``\vec R + \vec\tau_a`` and not only
the origin ``\vec R`` of the cell each atom belongs to. This is a gauge choice:
the two conventions are related by the atom-dependent rephasing
``e^{-2i\pi \vec q\cdot\vec\tau_a}``. The atomic-position gauge is adopted
(since version 0.12) because it makes the Fourier transformed quantities smooth
functions of ``\vec q``, which is much better suited for interpolating
quantities in q space. The price to pay is that quantities in this gauge are
**not periodic in the Brillouin zone**:

```math
\tilde v_a(\vec q + \vec G) = e^{-2i\pi \vec G\cdot \vec\tau_a}\, \tilde v_a(\vec q)
```

This must be carefully accounted for whenever a q point is folded back into the
grid by a reciprocal lattice vector ``\vec G`` (this occurs when applying symmetries,
see below).

Note the sign of the Fourier and the normalization prefactor. 
This convention allows for correctly transforming the matrices, however, it introduces a size inconsistency on the vectors.
If we have a periodic vector in the cell, its ``q`` fourier transformed counterpart will be ``\sqrt {N_q}`` times
the value in the primitive cell. So be carefull when extracting ``\Gamma`` point data from periodic vectors.

With this convention, we recover the standard rule for the matrices.

```math
\tilde \Phi_{ab}(\vec q) = \sum_{\vec R} e^{2\pi i \vec q\cdot (\vec R + \vec\tau_a - \vec\tau_b)}\Phi_{a + \vec R;b}
```

```math
\Phi_{ab} = \frac{1}{N_q} \sum_{\vec q}
\tilde\Phi_{ab}(\vec  q) e^{2i\pi \vec q\cdot[(\vec R(a) + \vec\tau_a) - (\vec R(b) + \vec\tau_b)]}
```

As for the vectors, the matrix in this gauge is not periodic in the Brillouin zone:

```math
\tilde\Phi_{ab}(\vec q + \vec G) = e^{2i\pi \vec G\cdot(\vec\tau_a - \vec\tau_b)}\, \tilde\Phi_{ab}(\vec q)
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
If the coordinates passed as `v_sc` are absolute values of positions,
then the ``\Gamma`` value of the fourier transform will be shifted by a total translation which is the average of the equilibrium positions.
This can be avoided by passing `absolute_positions = true` (which subtracts the
equilibrium positions ``\vec R + \vec\tau_a`` before transforming),
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

(where, in the atomic-position phase gauge, the folding by ``G`` introduces the block
phase factor ``e^{2i\pi \vec G\cdot(\vec\tau_a - \vec\tau_b)}``, see `impose_hermitianity_q!`),
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

The application of a symmetry in q space can be performed by considering how the force-constant matrix transforms in real space under a symmetry operation ``\{S | \vec v\}`` (rotation ``S`` plus fractional translation ``\vec v``, in crystal coordinates).

In the atomic-position phase gauge the derivation is particularly simple.
The symmetry maps the equilibrium position of the atom ``a`` in the cell ``\vec R``
into the equilibrium position of the atom ``s(a)`` in another cell:

```math
S(\vec R + \vec\tau_a) + \vec v = \vec R\,' + \vec \tau_{s(a)}
```

Since the phase factors of the Fourier transform are computed exactly from these
equilibrium positions, the phases follow the atoms through the symmetry
operation, and the transformation of the dynamical matrix takes the form

```math
S[\tilde\Phi_{ab}(\bm q)] = S^\dagger\, \tilde\Phi_{s(a)s(b)}(S_\text{recip}\bm q)\, S,
\qquad
\bm S_\text{recip} = \left(\bm S^{-1}\right)^\dagger
```

with **no phase factor** associated with the fractional translations: the phases
``e^{-2i\pi (S_\text{recip}\bm q)\cdot \vec v}`` picked up by the two displacement
vectors cancel between the two atomic indices (this is one of the advantages of
this gauge with respect to the lattice one, where a phase factor
``e^{2\pi i \bm q\cdot(\bm t_{s,a} - \bm t_{s,b})}`` involving the unit-cell
translations ``\bm t_{s,a} = S\vec\tau_a + \vec v - \vec\tau_{s(a)}`` appears).

However, in this gauge the dynamical matrix is **not periodic** in the reciprocal
lattice. The vector ``S_\text{recip}\bm q`` may fall outside the q grid, and it is
folded back into the grid point ``\bm q_{\text{grid}}`` by a reciprocal lattice
vector ``\bm G``:

```math
S_\text{recip}\bm q = \bm q_{\text{grid}} + \bm G
```

The folding introduces the phase factor

```math
\tilde\Phi_{s(a)s(b)}(\bm q_\text{grid} + \bm G) = e^{2\pi i \bm G\cdot(\vec\tau_{s(a)} - \vec\tau_{s(b)})}\, \tilde\Phi_{s(a)s(b)}(\bm q_\text{grid})
```

which is computed and applied automatically by `apply_symmetry_matrixq!`.
For this reason, the atomic positions inside the primitive cell (in crystal
coordinates) must be provided when initializing `SymmetriesQSpace`, and they
must be the same positions employed in the Fourier transform (choosing a
different periodic image of an atom changes the gauge).

For vectors, the fractional translation phase does not cancel, and the
transformation reads

```math
S[\tilde v_{a}(\bm q)] = e^{-2\pi i (S_\text{recip}\bm q)\cdot \vec v}\,
e^{2\pi i \bm G\cdot \vec\tau_{s(a)}}\, S\, \tilde v_{a}(\bm q)
```

These phases are applied by `apply_symmetry_vectorq!` when the positions, the
q points and the fractional translation are provided (they are automatically
provided when using the general interface `rotate_vector!`).

The application of symmetries is handled by the general function `rotate_vector!` and `rotate_dynamical_matrix!` or `rotate_matrix!` that works exactly like for real space symmetries, with the same general interface. However, we also provide specific q-space only functions. Note that, while the `rotate_*` functions works in cartesian space, the following one expects symmetries in real space.

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


