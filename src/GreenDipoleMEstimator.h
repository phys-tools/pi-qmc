//$Id: GreenDipoleMEstimator.h,v 1.2 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __GreenDipoleMEstimator_h_
#define __GreenDipoleMEstimator_h_
#include "stats/BlitzArrayBlkdEst.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <blitz/array.h>
#include <string.h>
#include <iostream>
#include <fftw3.h>
class Paths;
class Action;
class DoubleAction;
class SimulationInfo;
class Species;
class MPIManager;
/** dipole moment estimator using  temperature Green function method.
 *  @author Daejin Shin  */
class GreenDipoleMEstimator : public BlitzArrayBlkdEst<1>, public LinkSummable {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<std::complex<double>,1> CArray;
   typedef blitz::Array<std::string,1> SArray;

  /// Constructor.
  GreenDipoleMEstimator(const SimulationInfo& simInfo, const int index,
		        MPIManager *mpi);
  /// Virtual destructor.
  virtual ~GreenDipoleMEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int nfirst);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of coulomb energy estimator.
  /// Clear value of dipole energy estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// Parameters;
  const int npart, nslice;
  /// Tau
  const double tau;
  /// Dipole Moment.
  double dipole;
  /// The charges of the particles.
  Array q;
  /// The component of the dipole moment to measure (0->x, 1->y, 2->z)
  int index;
  SArray names;
  CArray temp;
  /// FFT plan
  fftw_plan fwd, rev;
  MPIManager *mpi;
};

#endif
