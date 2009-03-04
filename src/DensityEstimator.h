// $Id: DensityEstimator.h 12 2009-02-07 23:32:51Z john.shumwayjr $
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
#ifndef __DensityEstimator_h_
#define __DensityEstimator_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "stats/BlitzArrayBlkdEst.h"
#include "stats/MPIManager.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <blitz/array.h>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Action.h"
#include "Paths.h"
#include <blitz/tinyvec-et.h>
#include <vector>
/** Calculates single particle densities for many different geometries.
 *  Implements several options for studying fluctations.
 *  @version $Revision: 12 $
 *  @author John Shumway  */
class DensityEstimator : public LinkSummable {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  /// Base class for distance functions.
  class Dist {public: 
    virtual double operator()(const Vec &r, const SuperCell &cell)=0;
  };
  /// Distance taken from radial separation.
  class Radial : public Dist {
  public:
    Radial(int idir=-1, Vec center=Vec(0.)) 
      : mask(1.0), center(center) {
      if (idir!=-1) mask(idir)=0;
    }
    virtual double operator()(const Vec &r, const SuperCell &cell) {
      Vec delta=r-center;
      double radius2=0;
      for (int i=0; i<NDIM; ++i) radius2 += delta(i)*delta(i)*mask(i);
      return sqrt(radius2);
    }
    Vec mask, center;
  };
  /// Distance taken from cartesian position of particle.
  class Cart1 : public Dist { public:
    Cart1(int idim, double min=0.) : idim(idim), min(min) {};
    int idim;
    virtual double operator()(const Vec &r,
                              const SuperCell &cell) {
      return r[idim]-min;
    };
    double min;
  };
  typedef std::vector<Dist*> DistArray;
  typedef std::Blitz<double,1> Array;

  /// Constructor.
  DensityEstimator(const SimulationInfo& simInfo, const std::string& name,
      const Species &s1, int N, const Vec &min, const Vec &max,
      const IVec &nbin, const DistArray &dist, MPIManager *mpi); 

  /// Virtual destructor.
  virtual ~DensityEstimator();
  
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
  const int N;
  const Vec min;
  const Vec deltaInv;
  const IVec nbin;
  const DistArray dist;
  SuperCell cell;
  //Array temp;
  int ifirst, npart;
  MPIManager *mpi;
};
#endif
