# Acoustic sum rule

One fundamental aspect of symmetries of phonons is the so called acoustic sum rule.
This properties stems from the invariance of physical properties to global translations, and it is reflected to the fact that acoustic phonons, at the center of the Brilluin zone, have a vaninshing energy.

Therefore, the force constant matrix in real space, has always a kernel of dimension equal to the number of acoustic phonons: i.e. the number of operations under which the system is completely invariant.

The application of the acoustic sum rule depends on the dimensionality of the system, as well as the number of dimension. `AtomicSymmetries.jl` implements acoustic sum rules via defining a set of rules to compute it.

The standard acoustic sum rule is imposed by defining a structure `ASRConstraint!` which contains the dimension of the system. Then, the structure can act as a method to vectors and force constant matrices imposing the ASR condition on them.


```@docs
ASRConstraint!
```

It is possible to automatically identify translational modes with `translation_mask!` subroutine, which returns false to modes that are purely translational.
This is very important to selectively optimize only a subset of modes.

```@docs
translation_mask!
```
