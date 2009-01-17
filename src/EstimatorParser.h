// $Id: EstimatorParser.h,v 1.8 2008/02/08 00:50:30 jshumwa Exp $
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
#ifndef __EstimatorParser_h_
#define __EstimatorParser_h_
#include <string>
#include "XMLUnitParser.h"
#include "PairCFEstimator.h"
class EstimatorManager;
class SimulationInfo;
class Action;
class DoubleAction;
class MPIManager;
/// XML Parser for estimators.
/// @version $Revision: 1.8 $
/// @author John Shumway */
class EstimatorParser : public XMLUnitParser {
public:
  /// Constructor.
  EstimatorParser(const SimulationInfo&, const double tau, 
      const Action* action, const DoubleAction* doubleAction,MPIManager *mpi=0);
  /// Virtual destructor.
  ~EstimatorParser();
  /// Parse some xml.
  void parse(const xmlXPathContextPtr& ctxt);
  /// Return the Estimator object.
  EstimatorManager* getEstimatorManager() {return manager;}
private:
  /// The estimator manager.
  EstimatorManager* manager;
  /// General simulation information.
  const SimulationInfo& simInfo;
  /// The timestep.
  const double tau;
  /// The action.
  const Action* action;
  /// The double action.
  const DoubleAction* doubleAction;
  /// The MPI manager.
  MPIManager *mpi;
  /// Parser for pair correlation estimator.
  template<int N> PairCFEstimator<N>* parsePairCF(xmlNodePtr estNode,
                                                  xmlXPathObjectPtr ctxt);
};
#endif
