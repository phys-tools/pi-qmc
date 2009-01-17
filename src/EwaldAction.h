// $Id: EwaldAction.h,v 1.9 2006/10/18 17:08:18 jshumwa Exp $
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

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
#ifndef __EwaldAction_h_
#define __EwaldAction_h_
class MultiLevelSampler;
class SectionChooser;
class Paths;
class SuperCell;
#include "PairAction.h"
#include "Action.h"
#include <vector>
#include <blitz/tinyvec.h>
#include <complex>

/** Class for creating and evaluating Ewald evaluation of long-range action.
The most common Ewald summation technique splits a sum of @f$1/r@f$ 
potentials into short-range real space and long-range k-space terms,
@f[
\newcommand{\bfr}{\mathbf{r}}
\newcommand{\bfk}{\mathbf{k}}
\sum_{i<j} \frac{q_iq_j}{|\bfr_{ij}|}
=\sum_{i<j}\frac{q_iq_j\operatorname{erfc}(\kappa|\bfr_{ij}|)}{|\bfr_{ij}|}
+\frac{2\pi}{V}\sum_{\bfk\ne 0}\frac{e^{-\frac{|\bfk|^2}{4\kappa^2}}}{|\bfk|^2}
  \left|\sum_j q_j e^{i\bfk\cdot\bfr_j}\right|^2
-\frac{\kappa}{\sqrt\pi}\sum_j q_j^2
-\frac{2\pi}{3V}\left|\sum_j q_j\bfr_j\right|^2
@f]
Thus the following quantity is zero,
@f[
\newcommand{\bfr}{\mathbf{r}}
\newcommand{\bfk}{\mathbf{k}}
0=-\sum_{i<j}\frac{q_iq_j\operatorname{erf}(\kappa|\bfr_{ij}|)}{|\bfr_{ij}|}
+\frac{2\pi}{V}\sum_{\bfk\ne 0}\frac{e^{-\frac{|\bfk|^2}{4\kappa^2}}}{|\bfk|^2}
  \left|\sum_j q_j e^{i\bfk\cdot\bfr_j}\right|^2
-\frac{\kappa}{\sqrt\pi}\sum_j q_j^2
-\frac{2\pi}{3V}\left|\sum_j q_j\bfr_j\right|^2,
@f]
and we add it to the input actions.
@todo Make a new implementation using optimized Ewald breakup, as described
in <a href="http://dx.doi.org/10.1006/jcph.1995.1054">
Natoli and Ceperley, J. Comp. Phys. <b>117</b>, 171-178 (1995).</a>
@bug Only works for NDIM=3.
@version $Revision: 1.9 $
@author John Shumway. */
class EwaldAction : public Action {
  typedef std::vector<PairAction*> SRActionContainer;
  typedef SRActionContainer::iterator SRActIter;
  typedef SRActionContainer::const_iterator ConstSRActIter;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array1;
  typedef blitz::Array<Vec,1> VArray1;
  typedef std::complex<double> Complex;
  typedef blitz::Array<Complex,2> CArray2;
public:
  /// Constructor
  EwaldAction(const SimulationInfo&, const double rcut, const double kcut, 
              const int nreserve=1);
  /// Virtual destructor deletes all Action objects.
  virtual ~EwaldAction();
  /// Setup the long range and short range action.
  void setup();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
       double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  /// Add an action object.
  void addAction(PairAction* a) {actions.push_back(a);}
  /// Initialize for a sampling section.
  virtual void initialize(const SectionChooser&);
  /// Accept last move.
  virtual void acceptLastMove();
private:
  /// The SuperCell.
  const SuperCell &cell;
  /// Pointers to the short-range Action objects.
  SRActionContainer actions;
  /// Real space cutoff.
  const double rcut;
  /// K-space cutoff.
  const double kcut;
  /// Screening parameter.
  const double kappa;
  /// The timestep.
  const double tau;
  /// The number of particles.
  const int npart;
  /// The particle charges.
  Array1 q;
  /// The particle positions.
  mutable VArray1 pos;
  /// The self-energy term, @f$\frac{\kappa}{\sqrt\pi}\sum_j q_j^2@f$.
  double selfEnergy;
  /// Constants.
  static const double PI;
  /// Integer limits of k-space sums.
  IVec ikmax;
  /// Reciprical space lattice spacing.
  Vec deltak;
  /// Number of k-vectors.
  int totk;
  /// Stored values of @f$e^{-\frac{k^2}{4\kappa^2}}/k^2@f$.
  Array1 expoverk2;
  /// Arrays to tabulate @f$e^{ik_xx},e^{ik_yy},e^{ik_zz}@f$
  mutable CArray2 eikx, eiky, eikz;
  /// Calculate the long range action.
  double calcLongRangeUtau(VArray1& r) const;
  /// The the @f$2\pi/V@f$ prefactor on the k-space sum.
  double twoPiOverV;
};
#endif
