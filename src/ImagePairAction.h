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
#ifndef __ImagePairAction_h_
#define __ImagePairAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
class Paths;
class Species;
class SimulationInfo;
#include "PairAction.h"
#include <blitz/array.h>

/** Class for calculating pair action between particles.
* We tabulate and store the action as a radial function, with angular
* dependence stored as a power series. We take advantage of the extra
* symmetry of the coulomb action. The variables are
* @f[ q_{ij} = \frac{r_{ij} + r'_{ij}}{2} @f]
* and
* @f[ s^2_{ij} = \frac{|r_{ij}- r'_{ij}|^2}{q^2}. @f]
* The action and its tau derivative are stored as radial grids,
* @f[ u(q,s)=\sum_{k=0}^{N_{order}} u_k(q)s^{2k} @f]
* and
* @f[ \dot{u}(q,s)=\sum_{k=0}^{N_{order}} \dot{u}_k(q)s^{2k}. @f]
* The derivatives of the action (for the VirialEnergyEstimator) are given by
* @f[ \;-\nabla_i u(r_{ij},r'_{ij}) = -
* \frac{\partial u}{\partial q} \nabla q -
* \frac{\partial u}{\partial (s^2)} \nabla s^2, @f]
* where
* @f[ \nabla_i q_{ij} = \frac{\mathbf{r}_{ij}}{2r_{ij}} @f]
* and
* @f[ \nabla_i s^2_{ij} = \frac{\mathbf{r}-\mathbf{r}'_{ij}}{q_{ij}^2}-
*                         \frac{s^2\mathbf{r}_{ij}}{q_{ij}r_{ij}}. @f]
* The partial derivatives of the action are straightforward to calculate:
* @f[ \frac{\partial u}{\partial q}  
*     =\sum_{k=0}^{N_{order}} u'_k(q)s^{2k} @f]
* and
* @f[ \frac{\partial u}{\partial (s^2)}  
*     =\sum_{k=1}^{N_{order}} k u_k(q)s^{2(k-1)}. @f]
* @version $Revision$
* @todo Make a version for non-coulomb actions.
* @author John Shumway. */
class ImagePairAction : public PairAction {
public:
  /// Typedefs.
  typedef blitz::TinyVector<int,NDIM> IVec;
  /// Construct by providing the species and dmu filename.
  ImagePairAction(const Species&, const Species&, const std::string& filename,
    const SimulationInfo&, const int norder, const IVec nimage,
    const bool isDMD, int exLevel);
  /// Construct by providing the species and dmu filename.
  ImagePairAction(const Species&, const Species&, const EmpiricalPairAction&,
    const SimulationInfo&, const int norder, const IVec nimage,
    const double rmin, const double rmax, const int ngpts, int exLevel);
  /// Virtual destructor.
  virtual ~ImagePairAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
 /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  friend class EwaldAction;
private:
  /// The number of images.
  const IVec nimage;
};
#endif
