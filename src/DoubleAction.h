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
#ifndef __DoubleAction_h_
#define __DoubleAction_h_
class DoubleMLSampler;class DoubleDisplaceMoveSampler;
class DoubleSectionChooser;
class Paths;
#include <cstdlib>
#include <blitz/array.h>

/// Virtual base class for calculating the fermi action that depends
/// on two slices
/// Action objects can compute differences in action from
/// a MultiLevelDoubleSampler, the total action, and the action and
/// its derivatives at a bead.
/// @version $Revision$
/// @author John Shumway.
class DoubleAction {
public:
  /// Typdefs and constants.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Vec,1> VArray;
  /// Virtual destructor.
  virtual ~DoubleAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const DoubleMLSampler&,
                                     const int level) {return 0;}
  /// Calculate the difference in action.
  virtual double getActionDifference(const Paths&, const VArray &displacement,
    int nmoving, const IArray &movingIndex, int iFirstSlice, int iLastSlice) {
    return 0;
  }
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const=0;
  /// Calculate action and derivatives at a bead (defaults to no
  /// contribution).
  virtual void getBeadAction(const Paths&, int ipart, int islice,
          double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const=0;
  /// Initialize for a sampling section.
  virtual void initialize(const DoubleSectionChooser&) {};

  /// Accept last move.
  virtual void acceptLastMove() {};
};
#endif
