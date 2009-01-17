// $Id: AngularMomentumEstimator.h,v 1.2 2007/01/16 20:29:17 jshumwa Exp $
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
#ifndef __AngularMomentumEstimator_h_
#define __AngularMomentumEstimator_h_
#include "stats/ScalarEstimator.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <string.h>
#include <math.h>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
class Paths;
class PhaseModel;
class SimulationInfo;
/** Angular Momentum estimator.
 * Hard code for 2-dimensional systems.
 * @f[ L_z = \hbar (x \partial_y\phi_T - y \partial_x\phi_T) @f]
 *  @version $Revision: 1.2 $
 *  @author Daejin Shin  */
class AngularMomentumEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.
  AngularMomentumEstimator(const SimulationInfo& simInfo, PhaseModel*);
  /// Virtual destructor.
  virtual ~AngularMomentumEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of coulomb energy estimator.
  virtual double calcValue() {return value=angMtot/angMnorm;}
  /// Clear value of coulomb energy estimator.
  virtual void reset() {angMtot=angMnorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// AngularMomentum energy.
  double angM;
  /// The accumulated energy.
  double angMtot;
  /// The normalization of the accumualted energy.
  double angMnorm;
  /// The phase model.
  PhaseModel* phaseModel;
  /// The number of particles.
  const int npart;
  /// The number of slices.
  const int nslice;
  /// Arrays for coordinates and phase gradient.
  VArray r1, r2, grad;
};


#endif
