// $Id$
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
#ifndef __PairAction_h_
#define __PairAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
class Paths;
class Species;
class SimulationInfo;
class PairIntegrator;
#include "Action.h"
#include <blitz/array.h>

/** Class for calculating pair action between particles.
* We tabulate and store the action as a radial function, with angular
* dependence stored as a power series. We take advantage of the extra
* symmetry of the coulomb action. The variables are
* @f[ q_{ij} = \frac{r_{ij} + r'_{ij}}{2} @f]
* and
* @f[ s^2_{ij} = \frac{|r_{ij}- r'_{ij}|^2}{q^2}. @f]
* For general potentials (not 1/r), there is an additional coordinate
* @f[ z_{ij} = \frac{r_{ij}-r'_{ij}}{q}. @f]
* The action and its tau derivative are stored as radial grids,
* @f[ u(q,s)=\sum_{k=0}^{N_{order}}
*            \sum_{l=0}^k u_{kl}(q)s^{2k}z^{2l}. @f]
* and
* @f[ \dot{u}(q,s)=\sum_{k=0}^{N_{order}}
*                  \sum_{l=0}^k \dot{u}_{kl}(q)s^{2k}z^{2l}. @f]
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
* @author John Shumway. */
class PairAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<double,2> Array2;
  typedef blitz::Array<double,3> Array3;
  /// Helper class for constructing from emprical action.
  class EmpiricalPairAction{public: 
    virtual double u(double r, int iorder) const=0;
    virtual double utau(double r, int iorder) const=0;
  };
  /// Construct by providing the species and squarer pidmu filename.
  PairAction(const Species&, const Species&, const std::string& filename,
             const SimulationInfo&, const int norder, 
             const bool hasZ);
  /// Construct by providing the species and dmu filename.
  PairAction(const Species&, const Species&, const std::string& filename,
             const SimulationInfo&, const int norder, 
             const bool hasZ, const bool isDMD);
  /// Construct by providing the species and EmpiricalPairAction.
  PairAction(const Species&, const Species&, const EmpiricalPairAction&,
             const SimulationInfo&, const int norder, 
             const double rmin, const double rmax, const int ngpts,
             const bool hasZ);
  /// Construct by providing the species and PairIntegrator.
  PairAction(const Species&, const Species&, PairIntegrator&,
             const SimulationInfo&, const int norder, 
             const double rmin, const double rmax, const int ngpts);
  /// Virtual destructor.
  virtual ~PairAction() {}
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
  /// Write data tables to disk (defaults to speciesNames.dm[eu]).
    void write(const std::string &filename, const bool hasZ) const;
  friend class EwaldAction;
protected:
  /// The timestep.
  const double tau;
  /// The number of gridpoints.
  int ngpts;
  /// The first grid radial value.
  double rgridinv;
  /// The log of the ratio of consecutive grid points.
  double logrratioinv;
  /// The action grid (radial coord is logrithmic).
  Array3 ugrid;
  /// Evalute the diagonal action from the grid.
  double u00(double r) const;
  /// Evalute the off-diagonal action from the grid.
  double uk0(double q, double s2) const;
  /// Evalute the off-diagonal action from the grid.
  double uk0(double q, double s2, double z2) const;
  /// Evaluate the derivatives of the action from the grid.
  void uk0CalcDerivatives(double q, double s2, double &u,
            double &utau, double &uq, double &us2) const;
  /// Evaluate the derivatives of the action from the grid.
  void uk0CalcDerivatives(double q, double s2, double z2, double &u,
            double &utau, double &uq, double &us2, double &uz2) const;
  /// The species.
  const Species &species1, &species2;
  /// Indicies of first particles of each species.
  int ifirst1,ifirst2;
  /// Number of each particle for species in the interaction.
  int npart1,npart2;
  /// Order of the off-diagonal expansion.
  int norder;
  /// Flag for using the format from DMD downloaded from Ceperley site. 
  bool isDMD;
  /// Flag for including sum over z.
  bool hasZ;
  
};
#endif
