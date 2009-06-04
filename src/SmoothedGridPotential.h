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
#ifndef __SmoothedGridPotential_h
#define __SmoothedGridPotential_h
#include "MultiLevelSampler.h"
#include "DisplaceMoveSampler.h"
#include "SimulationInfo.h"
#include "Paths.h"
#include "SuperCell.h"
#include "Action.h"
#include "PeriodicGaussian.h"
#include <hdf5.h>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <fftw3.h>

class MultiLevelSampler;
class Paths;
template <int TDIM> class Beads;
class SimulationInfo;
class PeriodicGuassian;

/** Class for getting action and potential from a grid with quantum smoothing.
  * As described in section 3.5 of Feynman's <em>Statistical Mechanics</em>,
  * we can smooth the potential with a Gaussian spread @f$(12mk_BT)^{-1/2}@f$,
  * @f[
  * U(r;\tau) = \left({\frac{6m}{\pi\tau}}\right)^{3/2}
  *        \int V(r-r')\exp\left[-\frac{6m(r-r')^2}{\tau}\right]\;d^3r'
  * @f]
  * @todo Set up a potential table.
  * @todo Have simInfo contain the maximum number of levels.
  * @todo use advanced FFTW interface to only plan once (one FFT, one iFFT)
  * @todo use real DFT
  * @todo read in mass
  * @todo make scalar and vector mass routines
  * @bug Hard coded for SiGe ve and vh input.
  * @bug Hard coded for NDIM=3.
  * @version $Revision$
  * @author Matthew Harowitz. */
  
  /// Constants
  const double PI=acos(-1.);
  const double TWOPI=2.*PI;
  const double PI2=PI*PI;
  const double HATOEV=27.21138;
  const double EVTOHA=1./HATOEV;

class SmoothedGridPotential : public Action {
public:

  /// Typedefs.
  typedef blitz::TinyVector<double, NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef blitz::TinyVector<PeriodicGaussian*, NDIM> PGVec;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,NDIM> Array3;
  typedef blitz::Array<std::complex<double>, NDIM> CArray3;
  typedef blitz::Array<Array3,1> Array3Array;

  /// Constructor by providing an HDF5 file name and SimulationInfo.
  SmoothedGridPotential(const SimulationInfo& simInfo, const int maxLevel,
                        const std::string& filename);
  /// Destructor.
  ~SmoothedGridPotential();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int ilevel);
  virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
 /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int ilevel) const;
  /// Calculate the and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
          double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
  /// The number of particles.
  const int npart;
  /// The grid dimensions.
  IVec n;
  /// Inverse grid spacing.
  double b;
  /// The number of levels stored on grids.
  int nlevel;
  /// temporary electron grids.
  Array3 vetemp;
  ///temporary hole grids.
  Array3 vhtemp;
  /// The electron action grids.
  Array3Array vegrid;
  /// The hole action grids.
  Array3Array vhgrid;
  /// Evaluate the potential at a point.
  double v(Vec, const int ipart, int ilevel) const;
  /// The grid index that associates the particles and temperatures to the grids.
  IArray vindex;
  ///Peridic Gaussians for hole smoothing
  PGVec pgh;
  ///Peridic Gaussians for electron smoothing
  PGVec pge;
};
#endif
