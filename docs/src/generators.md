# Symmetry Generators

## Introduction

A central problem in computational physics is to represent tensors —
force constants, higher-order interatomic force constants (IFCs), Born effective charges, etc. —
that respect the symmetry of a crystal.
Rather than storing every component of a rank-``k`` tensor of size ``(d \cdot N_\text{at})^k``
and enforcing symmetry *a posteriori*, one can work in the symmetry-invariant subspace directly.
A **generator** is a basis element of this subspace: any tensor that is invariant under the
crystal symmetry group can be written as a linear combination of its generators.

This approach has two main advantages:

1. **Dimensionality reduction.** The number of independent generators ``n_\text{gen}`` is
   typically much smaller than ``(d \cdot N_\text{at})^k``, so describing a tensor by
   its ``n_\text{gen}`` coefficients is far more compact.
2. **Symmetry by construction.** Any linear combination of generators automatically
   satisfies all symmetry constraints, so no symmetrization step is ever needed after
   reconstruction.


## Mathematical background

### Symmetry operations on tensors

Consider a crystal with ``N_\text{at}`` atoms in ``d`` spatial dimensions.
A rank-``k`` tensor ``T`` has components indexed by ``k`` composite indices
``\mu_j = (\text{atom}_j, \alpha_j)`` where ``\alpha_j \in \{1, \ldots, d\}``
is the Cartesian direction.

A symmetry operation ``S`` of the crystal group ``G`` acts on a rank-``k`` tensor as

```math
\bigl[S(T)\bigr]_{\mu_1\ldots\mu_k}
= \sum_{\nu_1\ldots\nu_k}
  \prod_{j=1}^{k} R_{\nu_j \mu_j}\;
  T_{\hat S^{-1}(\nu_1)\ldots\hat S^{-1}(\nu_k)},
```

where ``R`` is the ``d \times d`` rotation (or improper rotation) matrix associated
with ``S``, and ``\hat S^{-1}`` denotes the atom permutation induced by the inverse of
the symmetry operation: atom ``a`` is mapped to atom ``\hat S^{-1}(a)``.
For rank 2 this reduces to the familiar ``S' \Phi S`` rotation of a force-constant
matrix with atom relabelling.

### The symmetrization projector

The **symmetrization projector**

```math
\mathcal{P} = \frac{1}{|G|} \sum_{S \in G} S
```

is a linear idempotent operator (``\mathcal{P}^2 = \mathcal{P}``) whose image is the
subspace of ``G``-invariant tensors.
Its rank, ``n_\text{gen} = \text{tr}(\mathcal{P})``, equals the dimension of the invariant subspace.

### Generators as a basis

A set of ``n_\text{gen}`` tensors ``\{g_i\}_{i=1}^{n_\text{gen}}``
that form a basis of ``\text{Im}(\mathcal{P})`` are called **generators**.
Any ``G``-invariant tensor ``T`` admits a unique decomposition

```math
T = \sum_{i=1}^{n_\text{gen}} \alpha_i \, g_i,
```

where the coefficients ``\alpha_i`` are obtained by projecting ``T`` onto the basis:
``\alpha_i = \langle g_i, T \rangle`` (since the generators are orthonormal under the
Frobenius inner product).

### Constructing generators

Each generator is obtained by **symmetrizing a seed tensor**.
A seed tensor ``e^{(\mu)}`` has a single nonzero entry equal to 1 at position
``\mu = (\mu_1,\ldots,\mu_k)`` (plus its permutation-symmetric copies if
full permutation symmetry is assumed).
Its symmetrization

```math
g^{(\mu)} = \mathcal{P}\bigl(e^{(\mu)}\bigr)
```

is either zero (the seed probes a symmetry-forbidden component) or a valid generator.
By iterating over all seeds and discarding zero and linearly dependent results, one
builds a complete orthonormal basis.

A key property simplifies the independence check: in crystal coordinates, the
inner product between any two normalized generators is exactly ``0`` or ``\pm 1``.
No Gram–Schmidt orthogonalization is needed.

### Permutation symmetry

When the tensor has **full permutation symmetry** among its ``k`` indices —
as is the case for the ``k``-th order energy derivatives (mixed partial derivatives commute) —
only sorted atom ``k``-tuples ``a_1 \leq a_2 \leq \ldots \leq a_k``
need to be considered as seeds.
This eliminates redundancy from the index permutation group.

### Centrosymmetric crystals and odd-rank tensors

Under the spatial inversion ``\mathcal{I}`` (with ``R = -\mathbb{I}_d``),
a rank-``k`` tensor transforms with a factor ``(-1)^k``.
If the crystal group contains inversion **and** every atom sits on an inversion
center (so the atom permutation is the identity), then the symmetrized
odd-rank tensor satisfies ``\mathcal{P}(T) = -\mathcal{P}(T) = 0``.
In this case, the invariant subspace for odd-rank fully-symmetric tensors
is trivially empty (zero generators).


## Compact representation

A full rank-``k`` tensor over ``N_\text{at}`` atoms in ``d`` dimensions
has ``(d \cdot N_\text{at})^k`` components.
The [`Generator`](@ref) struct avoids materializing this full tensor:
it stores only the nonzero ``d^k`` Cartesian blocks, one per distinct
target atom ``k``-tuple.
The memory cost is ``\mathcal O(n_\text{targets} \cdot d^k)``
instead of ``\mathcal O((d \cdot N_\text{at})^k)``,
which becomes essential for higher ranks and large supercells.

All generator blocks are stored in **Cartesian coordinates**.


## Rank-1 and rank-2 generators (index-based API)

For rank-1 (vectors) and rank-2 (matrices) tensors, a lightweight API is available
that represents each generator by a single integer index.

### Vectors

```julia
generators = get_vector_generators(symmetry_group, cell)
```

The generators are a `Vector{Int}` of seed indices. The full generator vector
can be recovered with

```julia
vector = zeros(Float64, n_modes)
get_vector_generator!(vector, generators[i], symmetry_group)
```

Projection and reconstruction:

```julia
coefficients = zeros(length(generators))
get_coefficients_from_vector!(coefficients, my_vector, generators, symmetry_group)

final_vector = similar(my_vector)
get_centroids_from_generators!(final_vector, generators, coefficients, symmetry_group)
```

### Matrices

```julia
generators = get_matrix_generators(symmetry_group, cell)
```

Projection and reconstruction:

```julia
coefficients = zeros(length(generators))
get_coefficients_from_fc!(coefficients, my_matrix, generators, symmetry_group, cell)

final_matrix = similar(my_matrix)
get_fc_from_generators!(final_matrix, generators, coefficients, symmetry_group, cell)
```

### Index-based generators API

```@docs
AtomicSymmetries.get_vector_generators
AtomicSymmetries.get_matrix_generators
AtomicSymmetries.get_vector_generator!
AtomicSymmetries.get_coefficients_from_vector!
AtomicSymmetries.get_centroids_from_generators!
AtomicSymmetries.get_coefficients_from_fc!
AtomicSymmetries.get_fc_from_generators!
```


## Arbitrary-rank generators (compact API)

The [`Generator`](@ref) struct generalizes the generator concept to tensors of
any rank ``k``.
Unlike the index-based API above, each generator directly stores its Cartesian
blocks, enabling fast contraction and reconstruction without ever materializing
the full ``(d \cdot N_\text{at})^k`` tensor.

### Finding generators

```julia
generators = get_tensor_generators(symmetry_group, cell; rank=2)
```

This returns a `Vector{Generator}`. For rank 2, the result is consistent
with `get_matrix_generators` (same number of independent generators,
same invariant subspace). For rank 3 and higher, this is the only available
interface.

### Symmetrizing a tensor

The function `symmetrize_tensor!` applies the symmetrization projector to a
rank-``k`` tensor in Cartesian coordinates, modifying it in place.
For rank 1 and 2, it delegates to the optimized `symmetrize_vector!` and
`symmetrize_fc!` implementations.

```julia
symmetrize_tensor!(tensor, cell, symmetry_group)
```

### Projection and reconstruction

Given a set of generators and a Cartesian tensor, the coefficients are
obtained by projection:

```julia
coeffs = zeros(length(generators))
get_coefficients_from_tensor!(coeffs, tensor_cart, generators, cell)
```

The tensor can be reconstructed from the coefficients:

```julia
tensor = zeros(Float64, ntuple(_ -> n_modes, rank))
reconstruct_tensor!(tensor, generators, coeffs, cell)
```

### Contracting generators

Generators can be contracted with vectors and matrices without building the
full tensor. This is useful, for example, to compute ``\Phi^{(3)}_{a\alpha} \cdot u``
(a rank-3 IFC contracted with a displacement vector, yielding an effective force-constant matrix).

Contract a rank-``k`` generator with a Cartesian vector along one index,
producing a rank-``(k-1)`` result:

```julia
contract_generator_vector!(result_blocks, result_target_atoms,
                           gen, vector_cart, contract_index)
```

Contract with a Cartesian matrix along two indices,
producing a rank-``(k-2)`` result:

```julia
contract_generator_matrix!(result_blocks, result_target_atoms,
                           gen, matrix_cart, idx1, idx2)
```

### Applying a single symmetry operation

To apply one symmetry operation to a rank-``k`` tensor in crystal coordinates:

```julia
apply_symmetry_tensor!(result, tensor_cryst, S, dim, irt)
```

This generalizes `apply_sym_fc!` (rank 2) and `apply_sym_centroid!` (rank 1).
The result is *added* to `result`, not overwritten.

### Arbitrary-rank generators API

```@docs
Generator
get_tensor_generators
symmetrize_tensor!
apply_symmetry_tensor!
reconstruct_tensor!
get_coefficients_from_tensor!
contract_generator_vector!
contract_generator_matrix!
AtomicSymmetries.contract_generators
```
