#Introduction

An exciton in a three-dimensional harmonic potential with no Coulomb 
interaction has the wavefunction

Ψ(r<sub>e</sub>, r<sub>h</sub>) 
= (2πσ<sub>e</sub><sup>2</sup>)<sup>3/4</sup>
  (2πσ<sub>h</sub><sup>2</sup>)<sup>3/4</sup>
  exp(
  -r<sub>e</sub><sup>2</sup> / 4 σ<sub>e</sub><sup>2</sup>
  -r<sub>e</sub><sup>2</sup> / 4 σ<sub>e</sub><sup>2</sup> )

Here σ is the width of the probability density,

σ = sqrt(ℏ / 4 m ω).

The exciton recombination matrix element is

I = ∫ Ψ(r,r) d<sup>3</sup>r

For the non-interacting exciton wavefunction, we find

|I|<sup>2</sup> = (2 σ<sub>e</sub> σ<sub>h</sub> 
  / (σ<sub>e</sub><sup>2</sup> + σ<sub>h</sub><sup>2</sup>))<sup>3</sup>

#Sample feature test

The input file in this directory simulates an exciton with
m<sub>e</sub> = 1.0, m<sub>h</sub> = 0.5,
ℏω<sub>e</sub> = 2.0, and ℏω<sub>e</sub> = 0.5.

For these values, we find
σ<sub>e</sub> = 0.5, σ<sub>h</sub> = 1.4142,
and |I|<sup>2</sup> = 0.2483.

The pi-qmc code outputs the average value of
1 / (|I|<sup>2</sup> + 1)
which should be 0.80108.

