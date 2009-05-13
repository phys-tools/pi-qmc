//$Id$
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
#ifndef __SpinFixedPhaseAction_h_
#define __SpinFixedPhaseAction_h_

#include "DoubleAction.h"
#include <blitz/array.h>
#include <blitz/tinymat.h>
#include <vector>
class SimulationInfo;
class Paths;
class Species;
class SpinPhaseModel;

/**  Class for fixed-phase action.
In the fixed-phase approximation, the phase of a trial function
@f$\rho_T@f$ (defined by the SpinPhaseModel) contributes to the hamiltonian,
@f[ \sum_{j=0}^N \frac{\hbar^2}{2m_j}\left|
\nabla_j\arg\rho_T(R,R')-\frac{e}{c}A_j(r_j)\right|^2. @f]
This hamiltonian acts on the same slice as @f$R@f$, but references 
coordinate positions on slice @f$R'@f$.

This action needs careful treatment, since the phase of the trial function 
is discontinuous for real-valued functions. Rather that relying solely on the 
value of the phase gradients a the ends of the time interval, we fit 
a cubic function to the phase gradients and values, and integrate
along a straight trajectory. Let @f$\partial_{\parallel}\phi@f$ 
be the partial derivative of @f$\phi@f$ along the direction from @f$\mathbf{r}'@f$ 
to @f$\mathbf{r}'@f$. We fit the phase and gradient of the phase to the value 
and gradient of the phase at points @f$\mathbf{r}'@f$ and @f$\mathbf{r}@f$,
using a cubic polynomial for the gradient along the link
and linear interpolation for gradient components perpendicular to the link
directions. We find the following action,
@f[
\newcommand{\br}{\mathbf{r}}
\newcommand{\bA}{\mathbf{A}}
\begin{split}
S_{\text{FP}} =  \frac{\tau\hbar^3}{2m}
\Big\{&
\frac{6}{5d^2}\Big[\phi(\br)-\phi(\br')\Big]^2
+\frac{1}{15}\Big[
2 [\partial_\parallel \phi(\br)]^2+
2 [\partial_\parallel \phi(\br')]^2
-\partial_\parallel\phi(\br)\partial_\parallel\phi(\br')+
\frac{3}{d}\big(-\phi(\br)+\phi(\br')\big)
\big(\partial_\parallel\phi(\br)+\partial_\parallel\phi(\br')\big)\Big]\\
&+\frac{1}{3}\Big[\Big|\nabla_\perp\phi(\br)-\frac{e}{\hbar c}\bA_\perp(\br)
+\nabla_\perp\phi(\br')-\frac{e}{\hbar c}\bA_\perp(\br')\Big|^2
-\Big(\nabla_\perp\phi(\br)-\frac{e}{\hbar c}\bA_\perp(\br)\Big)
\cdot\Big(\nabla_\perp\phi(\br')-\frac{e}{\hbar c}\bA_\perp(\br')\Big)\Big]\\
&+\frac{e}{\hbar cd}\Big(A_\parallel(\br)-A_\parallel(\br')\Big)
\Big(\phi(\br)-\phi(\br')\Big)
+\frac{e^2}{6\hbar^2c^2}\Big(A_\parallel(\br)-A_\parallel(\br')\Big)
\Big(\partial_\parallel\phi(\br)-\partial_\parallel\phi(\br')\Big)\\
&+\frac{e^2}{3\hbar^2c^2}\Big(
\big(A_\parallel\phi(\br)+A_\parallel(\br')\big)^2
-A_\parallel(\br)A_\parallel(\br')\Big)\Big\}\\
\end{split}
@f]
For real valued functions, this action is @f$3\tau\hbar^3\pi^2/5md^2@f$ when
the propagator crosses a node, with an average value @f$(3/5)\hbar\pi^2@f$.

@version $Revision$
@author Daejin Shin and John Shumway */
class SpinFixedPhaseAction : public DoubleAction {
public:
  /// Constants and typedefs.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<double,4> SVec;
  typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Vec,1> VArray;
  typedef blitz::Array<SVec,1> SArray;
  typedef blitz::Array<double,2> Matrix;
  typedef blitz::Array<Vec,2> VMatrix;
  typedef blitz::Array<SVec,2> SMatrix;
  typedef blitz::Array<Mat,2> MMatrix;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Constructor.
  SpinFixedPhaseAction(const SimulationInfo&, const Species&,
                  SpinPhaseModel*, const int maxlevel=10);
  /// Destructor.
  ~SpinFixedPhaseAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const DoubleMLSampler&,int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const;
  /// Initialize for a sampling section.
  virtual void initialize(const DoubleSectionChooser&);
  /// Accept last move.
  virtual void acceptLastMove();
private:
  /// The time step.
  const double tau;
  /// Total number of particles.
  const int npart;
  /// Number of particles of this type of fermion.
  const int nSpeciesPart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  /// Mass.
  const double mass;
  /// Number of slices.
  int nslice;
  /// Storage for the current slice.
  mutable VArray r1,r2;
  /// Storage for spins in the current slice.
  mutable SArray s1,s2;
  /// Flag for checking if this species is being moved.
  bool notMySpecies;
  /// Stored phase values.
  mutable Array phi, newPhi;
  /// Stored gradients of phase values.
  mutable std::vector<VArray> gradPhi1, gradPhi2, newGradPhi1, newGradPhi2;
  /// Stored spin gradients of phase values.
  mutable std::vector<SArray> sgradPhi1, sgradPhi2, newSGradPhi1, newSGradPhi2;
  /// Stored vector potential values.
  mutable std::vector<VArray> vecPot1, vecPot2, newVecPot1, newVecPot2;
  /// Determinants.
  mutable Array dmValue, newDMValue;
  /// Distances to the nodes.
  mutable Array dist, newDist;
  /// Array to cache forces for getBeadAction.
  mutable VArray force;
  /// The SpinPhaseModel.
  SpinPhaseModel* phaseModel;
  /// Function for evaluating the cubic interpolated action.
  double action(const double phi0, const double phi1, 
                const VArray& r1, const VArray& r2,
                const SArray& s1, const SArray& s2,
                const VArray& gradPhi01, const VArray& gradPhi11,
                const VArray& gradPhi02, const VArray& gradPhi12,
                const SArray& sgradPhi01, const SArray& sgradPhi11,
                const SArray& sgradPhi02, const SArray& sgradPhi12,
                const VArray& a01, const VArray& a11,
                const VArray& a02, const VArray& a12,
                const double tau, const bool reject=false) const;
  /// Pi.
  static const double PI;
  /// The speed of light.
  static const double C;
};
#endif
