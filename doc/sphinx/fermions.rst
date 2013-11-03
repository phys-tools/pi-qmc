Fermions and the Fixed Node Approximation
=========================================

Exact Fermions
--------------

Example: Two fermions in a one-dimensional simple harmonic oscillator
`````````````````````````````````````````````````````````````````````
For two distinguishable particles in a haromonic confinining potential
in three dimensions, the imaginary time propagator is

.. math::
   :label: shok3d

   K(\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_1', \mathbf{r}_2';\tau)
   = \left(\frac{m\omega}{2\pi\hbar\sinh\omega\tau}\right)^3
   \exp\left(-\frac{m\omega((r_1^2 + r_2^2 + r_1'^2 + r_2'^2)
   \cosh\omega\tau
   -2\mathbf{r}_1\cdot\mathbf{r}_1'
   -2\mathbf{r}_2\cdot\mathbf{r}_2')}
   {2\hbar\sinh\omega\tau}\right)

The partition function is the trace of the propagator for 
:math:`\tau = \beta\hbar`,

.. math::
   :label: shoz3d

   Z &= \operatorname{tr} K \\
   &= \int d\mathbf{r}_1^3 \int d\mathbf{r}_1^3 \,
   K(\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_1, \mathbf{r}_2; \beta\hbar) \\
   &= \left(\frac{m\omega}{2\pi\hbar\sinh\beta\hbar\omega}\right)^3
   \int d\mathbf{r}_1^3 \int d\mathbf{r}_1^3 
   \exp\left(-\frac{2m\omega(\cosh\beta\hbar\omega -1)(r_1^2 +r_2^2)}
   {2\hbar\sinh\omega\tau}\right)\\
   &= \left(2(\cosh\beta\hbar\omega-1)\right)^{-3} \\
   &= \left[2\sinh\left(\frac{\hbar\omega}{2k_BT}\right)\right]^{-6}

For identical particles, we need to symmetrize the states for
fermions, or antisymmetrize the states for bosons.
The trace of the permuted propagator is,

.. math::

   Z_P &= \operatorname{tr} PK \\
   &= \int d\mathbf{r}_1^3 \int d\mathbf{r}_1^3 \,
   K(\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_2, \mathbf{r}_1; \beta\hbar) \\
   &= \left(\frac{m\omega}{2\pi\hbar\sinh\beta\hbar\omega}\right)^3
   \int d\mathbf{r}_1^3 \int d\mathbf{r}_1^3 
   \exp\left(-\frac{2m\omega(\cosh\beta\hbar\omega(r_1^2 +r_2^2)
   -\mathbf{r}_1\cdot\mathbf{r}_2)}
   {2\hbar\sinh\omega\tau}\right)

Next we change coordinates to 
:math:`\mathbf{R}=\frac{1}{2}(\mathbf{r}_1+\mathbf{r}_2)`
and :math:`\mathbf{r}=\mathbf{r}_1-\mathbf{r}_2`.
Then :math:`r_1^2+r_2^2 = R^2 + r^2/2`
and :math:`\mathbf{r}_1\cdot\mathbf{r}_2 = R^2 - r^2/2`,
and we find,

.. math::
   :label: shozp3d

   Z_P = \left[2\cosh\left(\frac{\hbar\omega}{2k_BT}\right)
   \sinh\left(\frac{\hbar\omega}{2k_BT}\right)\right]^{-3}
