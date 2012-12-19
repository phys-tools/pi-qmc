********
Overview
********

This is a quantum simulation program from the Shumway Research Group, which
focuses on applications to nanoscience and technology. Path integral Monte
Carlo (PIMC) simulates particles (often electrons and ions) by directly
sampling the canonical partition function. In the path integral formulation of
quantum statistical mechanics developed by Richard Feynman, particles get
represented by closed imaginary-time trajectories of length 
:math:`\hbar/k_BT`. PIMC
simulations are able to compute total energies, correlation functions, charge
distribution, and linear response functions for thermal equilibrium. As in many
quantum Monte Carlo methods, PIMC has efficient scaling with system size, often
order N\ :sup:`2` or N\ :sup:`3`.

Our application, pi-qmc, 
is well suited for modeling conduction electrons and holes
in quantum dots, quantum wires, and quantum wells. For quantum dots and wires,
we often generate realistic confining potentials using `qdot-tools`. We are
also testing and developing pi for ab initio calculations, but at this point
only hydrogen and helium systems work well.
