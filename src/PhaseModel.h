//$Id: PhaseModel.h,v 1.5 2006/10/18 17:08:19 jshumwa Exp $
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
#ifndef __PhaseModel_h_
#define __PhaseModel_h_

#include <blitz/array.h>

/** Base class for phase models used in the fixed phase method. 

Suppose we have density matrix of complex funtion, we can write 
@f$ \rho = |\rho|e^{i\phi}  @f$.
We need the derivative of density matrix,
@f[ \nabla\rho = (\nabla|\rho|)e^{i\phi} + i |\rho|e^{i\phi} \nabla \phi  
= \dfrac{\nabla|\rho|}{|\rho|}\rho+i(\nabla\phi)\rho @f]
Divide by @f$\rho@f$ on each side, we get
@f[ \dfrac{\nabla\rho}{\rho} =  \dfrac{\nabla|\rho|}{|\rho|}+i(\nabla\phi) 
= Re(\dfrac{\nabla\rho}{\rho}) + i Im(\dfrac{\nabla\rho}{\rho})@f] 
Now we can get the derivative of the phase.
@f[\boxed{\nabla\phi=Im \left(\dfrac{\nabla\rho}{\rho} \right)  }  @f]    
How to get @f$\frac{\nabla\rho}{\rho}@f$ in the simulation?
@version $Revision: 1.5 $
@author Daejin Shin and John Shumway */
class PhaseModel {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.
  PhaseModel(const int npart) : phi(0),gradPhi1(npart),gradPhi2(npart),
                                vecPot1(npart),vecPot2(npart) {}
  /// Virtual destructor.
  virtual ~PhaseModel() {}
  /// Evaluate the PhaseModel.
  virtual void evaluate(const VArray&, const VArray&, const int islice)=0;
  /// Get the value of the phase.
  double getPhi() const {return phi;}
  /// Get the value of the gradient of the phase.
  const VArray& getGradPhi(const int i) const {return (i==1)?gradPhi1:gradPhi2;}
  /// Get the value of the  vector potential.
  const VArray& getVecPot(const int i) const {return (i==1)?vecPot1:vecPot2;}
protected:
  /// The value of the phase.
  double phi;
  /// The gradient of the phase.
  VArray gradPhi1,gradPhi2;
  /// The vector potential.
  VArray vecPot1,vecPot2;
};
#endif
