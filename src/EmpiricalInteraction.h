// $Id$
/*  Copyright (C) 2008 John B. Shumway, Jr.

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
#ifndef __EmpriricalInteraction_h_
#define __EmpriricalInteraction_h_
class MultiLevelSampler;
class Paths;
class Species;
class SimulationInfo;
#include "PairAction.h"
#include <blitz/array.h>

/** Class for setting up empirical pair action between particles.
  * Right now it just uses the primative approximation.
  * We may add methods for squaring or parsing formulas later.
  * @author John Shumway. */
class EmpiricalInteraction : public PairAction::EmpiricalPairAction {
public:
  // Base class for empirical potentials.
  class Potential {
  public:
    virtual double operator()(double r) const=0;
  };

  //Empirical short range potential for trapped atoms.
  class Cosh2Potential : public Potential {
  public: 
    Cosh2Potential(double v0, double kappa);
    virtual double operator()(double r) const;
    double v0, kappa;
  };

  virtual double u(double r, int iorder) const;
  virtual double utau(double r, int iorder) const;

  //Construct by providing an emprical potential and timestep.
  EmpiricalInteraction(const Potential& v, const double tau);
private:
  const Potential& v;
  const double tau;
};
#endif
