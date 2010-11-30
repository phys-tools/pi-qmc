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
#ifndef __ModelSampler_h_
#define __ModelSampler_h_
class Action;
class ActionChoiceBase;
class Paths;
class AccRejEstimator;
class MPIManager;

#include "Algorithm.h"
#include <vector>
#include <blitz/array.h>
#include <iostream>

/** Class to sample different action models.
  @author John Shumway
*/
class ModelSampler : public Algorithm {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.
  ModelSampler(Paths&, Action*, ActionChoiceBase*, const MPIManager* mpi);
  /// Destructor.
  virtual ~ModelSampler();
  /// Run method, performs nrepeat samplings with probability ifreq.
  virtual void run();
  /// Get a pointer to the accept/reject statistic estimator.
  /// (You are responsible for deleting this new object.) 
  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);
protected:
  /// A reference to the paths.
  Paths& paths;
  /// The action to be evaluated during the move.
  Action *action;
  /// The action to be evaluated during the move.
  ActionChoiceBase *actionChoice;
  /// The number of action models to choose from.
  const int nmodel;
  /// A pointer to the accept-reject estimator.
  AccRejEstimator* accRejEst;
  /// A pointer to the MPI manager, zero if MPI is not used.
  const MPIManager* mpi;
  /// Method to atempt a Monte Carlo move, return true if accepted.
  virtual bool tryMove();
};
#endif
