*************************************
Example: A Simple Harmonic Oscillator
*************************************

Introduction: Quantum Simple Harmonic Oscillator
================================================

Pi-qmc is a computer program in the language of C++ that predicts the behavior of quantum, or atomically small, particles. These particles do not always behave in the same manner as larger bodies and must be studied using unique equations. Quantum particles behave more and more like larger bodies as their energy increases.

Pi-qmc is able to easily predict the motion of these particles by summing the potential motions found using a random walk where the motions are predicted using statistical equations. Normally, the behavior of quantum systems would have to be found using a far more difficult process where the energies of each eigenstate are summed and divided by the number of possible states in order to normalize the function and find the average energy. A simple example of a quantum system where this process could be applied is the Quantum Harmonic Oscillator (QHC). The mathematics behind this process start with calculating the total energy of the system.

Background
----------

The total energy, or Hamiltonian, of a Quantum Harmonic Oscillator is represented by the function

.. math:: 
   :label: sho_hamilt

   \hat{H}=\frac{\hat{\rho} ^{2}}{2m}+(1/2)m\omega^{2}\hat{x}^{2}

The kinetic energy of the quantum harmonic oscillator is represented by the first term while the stored energy is represented by the second term. This may seem complex at first glance, but upon further inspection, it becomes apparent that this Hamiltonian is nearly identical to Hamiltonian of a regular harmonic oscillator which is represented by the function

.. math:: 
   :label: sho_rearrangedhamilt

   E=\frac{1}{2}mv^{2}+\frac{1}{2}kx^{2}

In the Hamiltonian of a Quantum Harmonic Oscillator, the mass multiplied by the velocity can be substituted in place of the momentum and the first term can be rearranged to resemble the first term of the regular harmonic oscillator’s Hamiltonian. Similarly, the square root of the spring constant divided by the mass can be substituted in for the wavelength and the second term can be rearranged to resemble the second term of the regular harmonic oscillator. After both terms are rearranged, it becomes clear that the two Hamiltonians are identical.


Eigenstates of the Simple Harmonic Oscillator
---------------------------------------------

The Eigenstate of Quantum Harmonic Oscillator is the energy level it occupies. These energy levels are represented by whole numbers (n=0, n=1, n=2, …) and are separated by consistent amounts energy that increase with each consecutive energy level.  As the eigenstate increases, the period of the quantum harmonic oscillator decreases and the wave behaves differently. The energy of the eigenstate, or eigenvalue can be found using the function

.. math:: 
   :label: sho_eigenstatesol

   E_{n}=\hbar\omega(n+\frac{1}{2})

According to Heisenberg’s Uncertainty Principle, the momentum and the position of a quantum particle cannot be known at the same time. Therefore, we can surmise that the energy of the eigenvalue will never be equal to zero since the no energy means no movement which, in turn, means that both the momentum and the position of the quantum particle will be known. When we substitute in a 0 for the eigenvalue, we find that the function simplifies to

.. math:: 
   :label: sho_zeroenergy

   E_{0}=\frac{1}{2}\hbar\omega

which represents the lowest possible energy of the Eigenstate and agrees with Heisenberg’s Uncertainty principle.
When the eigenstate is multiplied by the Hamiltonian, the resulting value is the eigenvalue as shown in the function

.. math:: 
   :label: sho_eigenvalue

   \hat{H}\psi(x)=E\psi(x)


The Classical Partition Function
--------------------------------

.. index:: 
   pair: simple harmonic oscillator; partition function

The Quantum Mechanical Partition Function
-----------------------------------------

In one dimension, the partition function of the simple harmonic oscillator is

.. math::
   :label: shoz

   Z = \sum_{n=0}^\infty e^{-\beta\hbar\omega (n+\frac{1}{2})}
   = \left[2\sinh\left(\frac{\hbar\omega}{2k_BT}\right)\right]^{-1}

For N oscillators in D dimensions, the partition function is

.. math::
   :label: shoznd

   Z = \left[2\sinh\left(\frac{\hbar\omega}{2k_BT}\right)\right]^{-ND}

The Helmholtz free energy is

.. math::
   :label: sho_free_energy

   F = E - TS = -k_BT \ln(Z) 
   = ND k_BT \ln\left[2 \sinh\left(\frac{\hbar\omega}{2k_BT}\right)\right].

The total energy is

.. math:: 
   :label: sho_energy

   E = -\frac{d}{d\beta} \ln(Z) = ND \frac{\hbar\omega}{2} 
   \coth\left(\frac{\hbar\omega}{2k_BT}\right).

Normalizing the Function
------------------------

The function is now normalized by dividing by the partition function

.. math:: 
   :label: sho_norm
   z=\sum_{n=0}^{\infty}e^{-\frac{E_{n}}{kt}}

This function can be simplified to the form

.. math:: 
   :label: sho_simpnorm
   z=\frac{1}{2sin(h)\frac{\beta\hbar\omega}{2}} \beta=\frac{1}{kt}

Our final equation is now  

.. math:: 
   :label: sho_endnorm
   P_{n}=\frac{e^{\frac{-E_{n}}{kt}}}{z}

The Density Matrix
------------------

For one particle in a haromonic confinining potential
in three dimensions, the imaginary time propagator is

.. math::
   :label: shok3d

   K(\mathbf{r}, \mathbf{r}';\tau)
   = \left(\frac{m\omega}{2\pi\hbar\sinh\omega\tau}\right)^3
   \exp\left(-\frac{m\omega((r^2 + r'^2)
   \cosh\omega\tau - 2\mathbf{r}\cdot\mathbf{r}')}
   {2\hbar\sinh\omega\tau}\right)

The diagonal of the density matrix is the probability density,

.. math::
   :label: shorho3d

   \rho(\mathbf{r}) & = K(\mathbf{r}, \mathbf{r}';\tau) \\
   & = \left(\frac{m\omega}{2\pi\hbar\sinh\omega\tau}\right)^3
   \exp\left(-\frac{m\omega((r^2 + r'^2)
   \cosh\omega\tau - 2\mathbf{r}\cdot\mathbf{r}')}
   {2\hbar\sinh\omega\tau}\right)

Simulating a Simple Harmonic Oscillator
---------------------------------------

Setting up the Input file
-------------------------

You can simulate a simple harmonic oscillator with the following
input file, which we copy from the ``test/system/sho/`` directory.

.. literalinclude:: ../../test/system/sho/pimc.xml
   :language: xml
   :linenos:

Here we have set the temperature :math:`k_BT = 1.0` (line 5),
and we are simulating a harmonic oscillator with :math:`\hbar\omega = 0.5`
Hartrees (line 8).
From Eq. :eq:`sho_energy`, the energy of this oscillator in three
dimensions should be :math:`\frac{3}{2}\coth\frac{1}{4}\approx 3.06224`.

Calculating the Energy
----------------------

Calculating the Density
-----------------------

Calculating the Polarizability
------------------------------

