// $Id$
/*  Copyright (C) 2010 John B. Shumway, Jr.

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
#ifndef __EMARateAction_h_
#define __EMARateAction_h_
class MultiLevelSampler;
class SectionSamplerInterface;
class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
class Species;
#include "action/Action.h"
#include <cstdlib>
#include <blitz/array.h>
#include <vector>

/** Class to modify the action to allow for electron-hole recombination.
 * To sample recombination rates, we work with a larger ensemble that
 * contains recombining paths. This larger ensemble can be described by
 * an additional action term,
 * @f[ S_C = -\hbar\ln\left(1+Ce^{-\frac{S_E'-S_E}{\hbar}}\right), @f]
 * where
 * @f[ S_E'-S_E = \frac{m_h|r_{e,0}-r_{h,N_{T}-1}|^2}{2\Delta\tau}
 *              + \frac{m_h|r_{e,1}-r_{h,0}|^2}{2\Delta\tau}
 *              - \frac{m_h|r_{h,0}-r_{h,N_{T}-1}|^2}{2\Delta\tau}
 *              - \frac{m_h|r_{e,1}-r_{e,0}|^2}{2\Delta\tau}. @f]
 * @version $Revision$
 * @author John Shumway */
class EMARateAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<bool,1> BArray;

  EMARateAction(const SimulationInfo&, const Species&, const Species&,
    double C);
  virtual ~EMARateAction();

  virtual double getActionDifference(const SectionSamplerInterface&,
                                     const int level);
  /// Calculate the difference in action.
  virtual double getActionDifference(const Paths&, const VArray &displacement,
    int nmoving, const IArray &movingIndex, int iFirstSlice, int iLastSlice);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead (defaults to no
  /// contribution).
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;

    double testableGetActionDifference(const SectionSamplerInterface&,
                                       const int level);
private:
  /// The inverse of the timestep.
  const double invTau;
  /// The electron species.
  const Species& species1;
  /// The hole species.
  const Species& species2;
  /// The index of the recombining electron.
  const int index1;
  /// The index of the recombining hole.
  const int index2;
  /// The weight parameter "C" used to optimize sampling.
  const double C;
  /// The mass of the electron.
  Vec mass1;
  /// The mass of the hole.
  Vec mass2;
  /// The total number of slices in the path.
  const int nPathSlice;


};
#endif
