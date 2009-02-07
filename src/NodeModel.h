//$Id$
/*  Copyright (C) 2004-2008 John B. Shumway, Jr.

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
#ifndef __NodeModel_h_
#define __NodeModel_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <blitz/array.h>
#include <fstream>
class DoubleMLSampler;

/** Base class for node models used in the fixed node approximation. 
The NodeModel is a function @f$\rho_0(R(t),R(t+\beta/2)@f$ that is
used by FixedNodeAction.
 * @author John Shumway */
class NodeModel {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<Vec,1> VArray;
  typedef blitz::Array<Vec,2> VMatrix;
  /** Base class for matrix updates.
   The trial density matrix is
   @f[
   \rho_{ij} = \rho_T(r_j, r_j).
   @f] 
   The inverse matrix satisfies @f$1=\rho^{-1}\rho@f$, so a change in
   the density matrix changes the determinant by a factor,
   @f[
   \frac{|\rho'|}{|\rho|} = |\rho^{-1}\rho'|.
   @f]
   This matrix is the identity matrix, except for rows with moving
   particles. Using the notation @f$j_i@f$ to indicate the index j
   of the ith moving particle, we construct the nmoving-by-nmoving
   matrix
   @f[
   a_{ij} = \sum_k \rho^{-1}_{i_i k}\rho'_{k j_j}
   @f]
   so @f$|\rho'|/|\rho|=|a|@f$.
   We calculate and store the vectors @f$\phi_{ij}=\rho'_{i_ij}@f$ for
   each of the moving particles j.

   To calculate the distance to the node, we need to temporarily calculate
   a new updated matrix, The update formula for a changed row @f$j_j@f$ is,
   @f[
   (\rho^{-1})'_{ik} =
   \begin{cases}
   \frac{\rho_{ik}^{-1}}{b_{j_j,j}};& i=j_j\\
   \rho_{ik}^{-1} - \frac{\rho_{jk}^{-1} b_{ij}}{b_{j_j,j}};& i\ne j_j
   \end{cases}
   @f]
   where @f$ b_{ij} = \sum_k \rho^{-1}_{ik} \phi_{kj}@f$.
   */
  class MatrixUpdate {
    public:
    virtual double evaluateChange(const DoubleMLSampler&, int islice)=0;
    virtual void evaluateNewInverse(const int islice)=0;
    virtual void evaluateNewDistance(const VArray &r1, const VArray &r2,
            const int islice, Array &d1, Array &d2)=0;
    virtual void acceptLastMove(int nslice)=0;
  };
  /// Constructor.
  NodeModel(const std::string &name="");
  /// Virtual destructor.
  virtual ~NodeModel() {}
  /// Evaluate the density matrix function, returning the value.
  virtual double evaluate(const VArray &r1, const VArray &r2, 
                          const int islice)=0;
  /// Evaluate distance to the node in units of @f$ \sqrt{\tau/2m}@f$.
  /// Assumes that evaluate has already been called on the slice.
  virtual void evaluateDistance(const VArray &r1, const VArray &r2,
                                const int islice, Array &d1, Array &d2)=0;
  /// Evaluate the time-derivative of the distance to the 
  /// node in units of @f$ \sqrt{\tau/2m}@f$.
  /// Assumes that evaluate has already been called on the slice.
  virtual void evaluateDotDistance(const VArray &r1, const VArray &r2,
                                   const int islice, Array &d1, Array &d2) {;}
  /// Evaluate gradient of log of the distance to the node
  /// in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateGradLogDist(const VArray &r1, const VArray &r2,
           const int islice, VMatrix &gradd1, VMatrix &gradd2,
                             const Array &d1, const Array &d2)=0;
  /// Returns true if action depends on other particle coordinates.
  virtual bool dependsOnOtherParticles() {return false;}
  /// Numerically test the gradient of the log of the distance to the node.
  /// This routine can be called from evaluateGradLogDist from subclasses
  /// for debugging when preprocessor variable NODE_DIST_DEBUG set. 
  void testGradLogDist( const VArray &r1, const VArray &r2,
                        const int islice, VMatrix &gradd1, VMatrix &gradd2,
                        const Array &d1, const Array &d2, 
                        const int npart, const int ifirst);
  /// Get a pointer to the update object (null if no updates).
  MatrixUpdate* getUpdateObj() const {return updateObj;}
protected:
#ifdef NODE_DIST_DEBUG
  std::ofstream logFile; 
#endif
  /// Pointer to update object.
  MatrixUpdate *updateObj;
};
#endif
