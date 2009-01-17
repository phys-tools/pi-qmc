// $Id: DensDensEstimator.h,v 1.2 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __DensDensEstimator_h_
#define __DensDensEstimator_h_
#include "stats/BlitzArrayBlkdEst.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <fftw3.h>
class Paths;
class Action;
class DoubleAction;
class SimulationInfo;
class MPIManager;
/** DensDens estimator for a homogeneous system.
 *  @version $Revision: 1.2 $
 *  @author John Shumway  */
class DensDensEstimator : public BlitzArrayBlkdEst<3>, public LinkSummable {
public:
  typedef blitz::Array<std::complex<double>,3> CArray3;
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  DensDensEstimator(const SimulationInfo& simInfo,
    const Action*, const DoubleAction*,
    const int nbin, const int ndbin, MPIManager *mpi);
  /// Virtual destructor.
  virtual ~DensDensEstimator();
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
  const Action *action;
  const DoubleAction *doubleAction;
  ///
  const int npart, nslice, nbin, ndbin;
  const double tauinv, massinv, dx, dxinv;
  CArray3 temp;
  IArray ninbin,ninbinbuff;
  fftw_plan fwd, rev;
  MPIManager *mpi;
};

#endif
