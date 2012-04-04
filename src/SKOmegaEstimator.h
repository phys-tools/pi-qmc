// $Id: SKOmegaEstimator.h 393 2011-08-04 14:39:41Z john.shumwayjr $
/*  Copyright (C) 2011 John B. Shumway, Jr.

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
#ifndef __SKOmegaEstimator_h_
#define __SKOmegaEstimator_h_
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
#include "SuperCell.h"
#include "action/Action.h"
#include "Paths.h"
#include <blitz/tinyvec-et.h>
#include <vector>
#include <fftw3.h>

/** Calculates dynamic structure factor for a translationally invariant system.
 *  @version $Revision: 393 $
 *  @author John Shumway  */
class SKOmegaEstimator : public LinkSummable, 
                         public BlitzArrayBlkdEst<NDIM+3> {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef blitz::TinyVector<int, NDIM+3> IVecN;
  typedef blitz::Array<int,NDIM> IArrayNDIM;
  typedef blitz::Array<int,1> IArray;
  typedef std::complex<double> Complex;

  /// Constructor.
  SKOmegaEstimator(const SimulationInfo& simInfo, const std::string& name,
      const IVec &nbin, const IVecN &nbinN, int nstride, MPIManager *mpi); 

  /// Virtual destructor.
  virtual ~SKOmegaEstimator();
  
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
  const int nslice,nfreq,nstride,npart,nspec,ntot;
  const Vec min;
  const Vec deltaInv;
  const IVec nbin;
  SuperCell cell;
  const double tau;
  blitz::Array<Complex,NDIM+2> temp;
  IArray speciesIndex;
  fftw_plan fwd;
  MPIManager *mpi;
  blitz::Array<Complex,3>* temp_;
  blitz::Array<float,4>* value_;
};
#endif
