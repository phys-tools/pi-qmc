// $Id$
/*  Copyright (C) 2008-2009 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifndef __EwaldSum
#define __EwaldSum
class MultiLevelSampler;
class SectionChooser;
class Paths;
class SuperCell;
#include <vector>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <complex>

/** Class for creating and evaluating Ewald sums.
The Ewald sum technique enables fast evaluation of a periodic
sum of long-range potentials. At the heart of the method
is a function @f$f(r)@f$ that behaves as @f$1/r@f$ for large r,
but is smooth for small r.
The traditional choice for this function is
@f$ f(r) = \operatorname{erf}(\kappa r)/r@f$.
The periodic sum of this smooth function may be evaluated
efficiently in k-space.
@f[
\begin{split}
\frac{1}{2}\sum_{\mathbf{L}}\sum_{ij}' 
q_iq_j f(|\mathbf{r}_i-\mathbf{r_j}+\mathbf{L}|)
= &\frac{1}{2V}\sum_{\mathbf{k}\ne 0} f(|\mathbf{k}|)
\left|\sum_j q_j e^{i\mathbf k\cdot \mathbf r_j}\right|^2\\
-&\frac{1}{2}\sum_j q_j^2 f(r\rightarrow 0)
+\frac{1}{2V}f(k\rightarrow 0)\left|\sum_j q_j \right|^2.
\end{split},
@f]
where the prime on the sum indicates that the i=j term is
omitted for L=0.
The left hand side of this equation may be subtracted from any
sum over pair actions or potentials to remove 1/r tails and
periodic images, then the quantity is easily added back in using
the k-space sum on the right hand side.
Two choices for f(r) are implemented as sub-classes: 
TradEwaldSum makes the traditional choice that uses an error function,
and OptEwaldSum uses a polynomial fit to connect smoothly to the
1/r tail and optimized to converge quickly in k-space.
@version $Revision$
@author John Shumway. */
class EwaldSum {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<Vec,1> VArray;
  typedef std::complex<double> Complex;
  typedef blitz::Array<Complex,2> CArray2;
  /// Constructor calcuates the k-vectors for a given rcut and kcut.
  EwaldSum(const SuperCell&, const int npart, 
           const double rcut, const double kcut);
  /// Virtual destructor.
  virtual ~EwaldSum();
  /// Returns @f$ f(r) @f$, used to cancel tails on actions or potentials.
  virtual double evalFR(const double r) const=0;
  /// Returns @f$ f(r\rightarrow 0) @f$, used in evalSelfEnergy.
  virtual double evalFR0() const=0;
  /// Returns @f$ f(k) @f$, used to set up vk array.
  virtual double evalFK(const double k) const=0;
  /// Returns @f$ f(k\rightarrow 0) @f$, used in evalSelfEnergy for interaction
  /// with neutralizing background if system has a net charge.
  virtual double evalFK0() const=0;
  /// Evaluate the long range sum.
  double evalLongRange(const VArray& r) const;
  /// Evaluate the self energy using evalFR0 and evalFK0 virtual methods.
  /// You must call this function again if you change the charge array.
  void evalSelfEnergy();
  /// Get a reference to the charge array.
  Array& getQArray() {return q;}
  /// Set the long range table using the evalFK virtual method.
  void setLongRangeArray();
protected:
  /// The SuperCell.
  const SuperCell &cell;
  /// The number of particles.
  const int npart;
  /// Real space cutoff.
  const double rcut;
  /// K-space cutoff.
  const double kcut;
  /// Self energy;
  double selfEnergy;
  /// The particle charges.
  mutable Array q;
  /// The particle positions.
  mutable VArray pos;
  /// Constants.
  static const double PI;
  /// Integer limits of k-space sums.
  IVec ikmax;
  /// Reciprical space lattice spacing.
  Vec deltak;
  /// Number of k-vectors.
  int totk;
  /// Stored values of k-space potential.
  Array vk;
  /// Arrays to tabulate @f$e^{ik_xx},e^{ik_yy},e^{ik_zz}@f$
  mutable CArray2 eikx, eiky, eikz;
  /// Calculate the long range part.
  double calcLongRangeUtau(VArray& r) const;
  /// The prefactor on the k-space sum, 1/2V.
  double oneOver2V;
};
#endif
