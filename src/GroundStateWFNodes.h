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
#ifndef __GroundStateWFNodes_h_
#define __GroundStateWFNodes_h_

#include "Action.h"
class Wavefunction;
#include <blitz/array.h>
#include <vector>
class SimulationInfo;
class Paths;
class SuperCell;
class Species;

/** Ground state nodes.
  * @version $Revision$
  * @author John Shumway */
class GroundStateWFNodes : public Action {
public:
  /// Typdefs and constants.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  GroundStateWFNodes(const Species& species, const Species& refSpecies);
  /// Virtual destructor.
  virtual ~GroundStateWFNodes() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&, int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const {
    return 0;}
  /// Calculate action and derivatives at a bead (no contribution).
  virtual void getBeadAction(const Paths&, int ipart, int islice,
          double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
    u=utau=ulambda=0; fm=0.; fp=0.;}
  /// Initialize for a sampling section.
  virtual void initialize(const SectionChooser&) {};
  /// Accept last move.
  virtual void acceptLastMove() {};
private:
  const int npart;
  const int ifirst;
  const int iref;
  bool isMoving0,isMoving1,isRefMoving;
  int indref,ind0,ind1;
};
#endif
