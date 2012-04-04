// $Id$
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#ifndef __DensDensEstimator_h_
#define __DensDensEstimator_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "stats/ArrayBlockedEstimator.h"
#include "stats/BlitzArrayBlkdEst.h"
#include "stats/MPIManager.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>
#include "SimulationInfo.h"
#include "Species.h"
#include "util/SuperCell.h"
#include "action/Action.h"
#include "Paths.h"
#include <blitz/tinyvec-et.h>
#include <vector>
#include <fftw3.h>
class Distance;
/** Calculates single particle densities for many different geometries.
 *  Implements several options for studying fluctations.
 *  @version $Revision$
 *  @author John Shumway  */
class DensDensEstimator : public LinkSummable, 
                          public BlitzArrayBlkdEst<2*NDIM+1> {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef std::vector<Distance*> DistArray;
  typedef blitz::Array<int,NDIM> IArrayNDIM;
  typedef std::complex<double> Complex;

  /// Constructor.
  DensDensEstimator(const SimulationInfo& simInfo, const std::string& name,
      const Species *s1, const Species *s2,
      const Vec &min, const Vec &max, const IVec &nbin,
      const IVecN &nbinN, const DistArray &dist, int nstride,
      MPIManager *mpi); 

  /// Virtual destructor.
  virtual ~DensDensEstimator();
  
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);

  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
      const int ipart, const int islice, const Paths &paths);
  
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  
  /// Clear value of estimator.
  virtual void reset();

  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths);

private:
  const int nslice,nfreq,nstride,ntot;
  const Vec min;
  const Vec deltaInv;
  const IVec nbin;
  const IVecN nbinN;
  const DistArray dist;
  SuperCell cell;
  const double tau;
  blitz::Array<Complex,NDIM+1> temp1;
  blitz::Array<Complex,NDIM+1> temp2;
  int ifirst1, npart1, ifirst2, npart2;
  fftw_plan fwd1, fwd2;
  MPIManager *mpi;
  blitz::Array<Complex,2>* temp1_;
  blitz::Array<Complex,2>* temp2_;
  blitz::Array<float,3>* value_;
};
#endif
