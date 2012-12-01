Dynamic Response Functions
==========================

We define a dynamic correlation function as

.. math::
   :label: chiAB

   \chi_{AB}(\tau) = -(1/\hbar)\langle A(\tau)B(0)\rangle,

where :math:`A` and :math:`B` are operators.

Density-density response
------------------------

Since pi-qmc uses a position basis, we often collect density fluctuations in
real space. However, most textbook descriptions of density fluctioans are in
k-space, and results for homogeneous systems are often best represented in
k-space. Here we give a brief summary of common definitions for pedagogical
purposes. For simplicity we write all formulas for spinless particles.

The dimensionsless Fourier transform of the density operator is (Eqs. 1.11 and
1.66 of Giuliani and Vignale)

.. math::
   :label: nq

   n_{\mathbf q} 
   &= \sum_j e^{-i\mathbf{k}\cdot\mathbf{r}_j} \\
   &= \sum_{\mathbf{k}} a^\dagger_{\mathbf{k}-\mathbf{q}}a_{\mathbf{k}}.

Note that :math:`n_0=N`,
the total number of particles. To get back to real-space density use

.. math::
   :label: nr

   n(\mathbf{r}) = \frac{1}{V} 
   \sum_{\mathbf{q}} n_{\mathbf{q}} e^{i\mathbf{q}\cdot\mathbf{r}}.

For each of these, we define frequency-dependent density operators,

.. math::
   :label: nriwn

   n(\mathbf{r},i\omega_n) 
   = \int_0^{\beta\hbar} n(\mathbf{r},\tau) e^{i\omega_n\tau} d\tau,

and

.. math::
   :label: nqiwn

   n_{\mathbf{q}}(i\omega_n) = \int_0^{\beta\hbar} n_{\mathbf{q}}
   e^{i\omega_n\tau} d\tau,

where :math:`i\omega_n = 2\pi ink_BT/\hbar`
are the Matsubara frequencies. Within the pi-qmc code, these
frequency-dependent densities are easily calculated with fast Fourier
transforms, which are most efficient when the number of slices is a power of
two.

Real-space response
```````````````````

The imaginary-frequency response of the density to an external perturbation is
given by (Ch 3.3 of Guiliani and Vignale),

.. math::
   :label: deltan

   \delta n(\mathbf{r},i\omega_n) 
   = \int d\mathbf{r} \chi_{nn}(\mathbf{r},\mathbf{r}'', i\omega_n)
   V_{\text{ext}}(\mathbf{r}', i\omega_n).

In k-space this takes the convienent form,

.. math::
   :label: deltanq

   \delta n(\mathbf{q}, i\omega_n) 
   = \sum_{\mathbf{q}'} \chi_{nn}(\mathbf{q}, \mathbf{q}', i\omega_n)
   V_{\text{ext}}(\mathbf{q}',i\omega_n).

where the external potential in k-space satisfies

.. math::
   :label: Vr

   V_{\text{ext}}(\mathbf{r}') 
   = \frac{1}{V} \sum_{\mathbf{q}'} V_{\text{ext}}(\mathbf{q}')
   e^{i\mathbf{q}'\cdot\mathbf{r}'},

and

Vext(**q**') = ∫ d**q**' e-i**q**'·**r**' Vext(**r**').

These response functions are related to imaginary-frequency dynamic correlation
functions,

χnn(**r**,**r**',iωn) = -(1/βħ²)〈n(**r**,iωn) n(**r**',-iωn)〉,

and

χnn(**q**,**q**',iωn) = -(1/βħ²V)〈n**q**(iωn) n-**q**'(-iωn)〉.

For a homogeneous system,

χnn(**q**,**q**',iωn) = -(1/βħ²V)〈n**q**(iωn) n-**q**(-iωn)〉 δ**qq**'.

Structure factor
````````````````

The dynamic structure factor S(**k**,iωn) measures the density response of the
system,

S(**k**,iωn) = -(V/ħN) χnn(**k**,**k**,iωn)

The static structure factor is defined for equal time, not for ωn → 0,

S(**k**) = (1/N) 〈n**k**(τ=0) n-**k**(τ=0)〉.

In terms of χnn(**q**,**q**',iω), the static structure factor is given by
(*prefactor is wrong*)

S(**k**) = -(V/ħN) Σωn χnn(**k**,**k**,iωn) e-iωnτ.

Polarizability
``````````````
