Example: A Simple Harmonic Oscillator
=====================================

Background
----------

Eigenstates of the Simple Harmonic Oscillator
`````````````````````````````````````````````

The Classical Partition Function
````````````````````````````````

.. index:: 
   pair: simple harmonic oscillator; partition function

The Quantum Mechanical Partition Function
`````````````````````````````````````````

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

The Density Matrix
``````````````````

Simulating a Simple Harmonic Oscillator
---------------------------------------

Setting up the Input file
`````````````````````````

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
``````````````````````

Calculating the Density
```````````````````````

Calculating the Polarizability
``````````````````````````````

