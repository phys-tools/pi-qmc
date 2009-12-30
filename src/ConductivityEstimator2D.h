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
#ifndef __ConductivityEstimator2D_h__
#define __ConductivityEstimator2D_h__
#include "stats/BlitzArrayBlkdEst.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <fftw3.h>
class Paths;
class SimulationInfo;
class MPIManager;

class ConductivityEstimator2D : public BlitzArrayBlkdEst<7>, public LinkSummable {
public:
  typedef blitz::Array<std::complex<double>,7> CArray7;
  typedef blitz::Array<std::complex<double>,4> CArray4;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  /// Constructor.
  ConductivityEstimator2D(const SimulationInfo& simInfo, const double xmin, const double xmax, 
    const double ymin, const double ymax,  const int nfreq, const int nxbin, const int nybin,
    const int nxdbin, const int nydbin, const int nstride, MPIManager *mpi);
  /// Virtual destructor.
  virtual ~ConductivityEstimator2D();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const blitz::TinyVector<double,NDIM>& start,
                          const blitz::TinyVector<double,NDIM>& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  const int npart, nslice, nfreq, nstride, nxbin, nybin, nxdbin, nydbin;
  const double tau, tauinv, massinv, dx, dy, dxinv, dyinv, xmin, xmax, ymin, ymax;
//  double dx[2], dxinv[2], min[2], max[2];
//  int nbin[2], ndbin[2];
  Array q;
  CArray7 temp;
  CArray4 jx,jy;
  fftw_plan fwdx,fwdy;
  MPIManager *mpi;
};

#endif
