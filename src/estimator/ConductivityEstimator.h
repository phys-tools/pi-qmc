// $Id$
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
#ifndef __ConductivityEstimator_h_
#define __ConductivityEstimator_h_
#include "stats/BlitzArrayBlkdEst.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <fftw3.h>
class Paths;
class SimulationInfo;
class MPIManager;
/// Conductivity estimator for an inhomogeneous system.
///
/// @f[
/// \chi_{jj}(x,x';\tau) = -\langle j(x,\tau) j(x',0)\rangle_\beta
/// @f]
///
/// @f[
/// j(x_0)=\frac{q\hbar}{2im}[\partial_x\delta(x-x_0)+\delta(x-x_0)\partial_x]
/// @f]
/// We find it simplest to apply the estimator between beads, in the
/// middle of the free particle propagator, @f$G_0@f$, over the 
/// imaginary-time interval @f$\Delta\tau = \tau-\tau'@f$ :
/// @f[
/// T_\tau j(x,(\tau+\tau')/2)G_0(x,\tau;x,\tau') = 
/// iq\frac{x-x'}{\Delta\tau}
/// \frac{e^{-\frac{m(x+x'-2x_0)^2}{\Delta\tau}}}{\sqrt{\pi\Delta\tau/2m}}
/// G_0(x,\tau;x,\tau')
/// @f]
/// Even with interactions, we use the same free-particle estimator, as
/// this estimator satisfies the continuity equation and has the correct
/// coarse-grained response.
///
/// @version $Revision$
/// @author John Shumway 
class ConductivityEstimator : public BlitzArrayBlkdEst<3>, public LinkSummable {
public:
  typedef blitz::Array<std::complex<double>,3> CArray3;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  /// Constructor.
  ConductivityEstimator(const SimulationInfo& simInfo,
    const int nfreq, const int nbin, const int ndbin, const int nstride, 
    MPIManager *mpi);
  /// Virtual destructor.
  virtual ~ConductivityEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const blitz::TinyVector<double,NDIM>& start,
                          const blitz::TinyVector<double,NDIM>& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  // Get value of kinetic energy estimator.
  //virtual double calcValue() {return value=etot/enorm;}
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  ///
  const int npart, nslice, nfreq, nbin, ndbin, nstride;
  const double tau, tauinv, massinv, dx, dxinv;
  Array q;
  CArray3 temp;
  fftw_plan fwd, rev;
  MPIManager *mpi;
};

#endif
