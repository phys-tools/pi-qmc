//$Id$
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
#ifndef __GroundStateSNode_h_
#define __GroundStateSNode_h_

#include <blitz/array.h>
#include "NodeModel.h"
class Species;

/** Simple nodal model for two S electrons in an atom.
* Simply define the nodal function as 
* @f[ (|r_1-r_{\text{ion}}|-|r_0-r_{\text{ion}}|)
*     (|r'_1-r'_{\text{ion}}|-|r'_0-r'_{\text{ion}}|), @f]
* where 1 and 0 are the two identical fermions.
* The nodal distance is
* @f[ \sqrt{(|r_1-r_{\text{ion}}|-|r_1-r_{\text{ion}}|)^2/2
*          +(|r_1-r_{\text{ion}}|-|r_1-r_{\text{ion}}|)^2/2}, @f]
* and the gradient of the log of the distance to the node is a 
* unit vector pointing radially inward or outward.
* The force on the center ion is determined by Newton's Third Law.
* @author John Shumway */
class GroundStateSNode : public NodeModel {
public:
  /// Constructor.
  GroundStateSNode(const Species& s, const int icenter, const double tau);
  /// Virtual destructor.
  virtual ~GroundStateSNode();
  /// Evaluate the density matrix function, returning value.
  virtual double evaluate(const VArray &r1, const VArray &r2, 
                          const int islice);
  /// Evaluate distance to the node in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateDistance(const VArray &r1, const VArray &r2,
                                const int islice, Array &d1, Array &d2);
  /// Evaluate gradient of log of the distance to the node
  /// in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateGradLogDist(const VArray &r1, const VArray &r2,
           const int islice, VMatrix &gradd1, VMatrix &gradd2,
                             const Array &d1, const Array &d2);
  /// Returns true if action depends on other particle coordinates.
  virtual bool dependsOnOtherParticles() {return true;}
private:
  /// The index of the first particle. 
  const int ifirst;
  /// The number of particles.
  const int npart;
  /// The index of the center particle (nucleus). 
  const int icenter;
  /// The unit for distance, @f$ \sqrt{\tau/2m}@f$.
  const double dunit;
};
#endif
