// $Id: EwaldImagePairAction.h 214 2009-12-10 22:39:46Z saadAK $
/*  Copyright (C) 2004-2006, 2009 John B. Shumway, Jr.

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
#ifndef __EwaldImagePairAction_h_
#define __EwaldImagePairAction_h_
class MultiLevelSampler;
class DisplaceMoveSampler;
class Paths;
class Species;
class SimulationInfo;
#include "PairAction.h"
#include <blitz/array.h>
#include <vector>

class EwaldImagePairAction :public PairAction {
public:
  /// Construct by providing the species and EmpiricalPairAction.
    EwaldImagePairAction(const Species&, const Species&, const EmpiricalPairAction&,
			 const SimulationInfo&, const int norder, 
			 const double rmin, const double rmax, const int ngpts,
			 const int nImages);
    
  /// Virtual destructor.
  virtual ~EwaldImagePairAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the difference in action.
  virtual double getActionDifference(const Paths&, const VArray &displacement,
    int nmoving, const IArray &movingIndex, int iFirstSlice, int nslice);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  void findBoxImageVectors(std :: vector<double>); 

  friend class EwaldAction;
protected:
  double sphereR;
  const int nImages;
  std :: vector<std :: vector<double> > boxImageVecs;
};
#endif
