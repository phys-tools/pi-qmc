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
#ifndef __PIMCParser_h_
#define __PIMCParser_h_
#include "parser/XMLUnitParser.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <vector>
class Paths;
class SuperCell;
class SectionChooser;
class DoubleSectionChooser;
class SimulationInfo;
class Algorithm;
class CompositeAlgorithm;
class Action;
class ActionChoiceBase;
class DoubleAction;
class EstimatorManager;
class ProbDensityGrid;
class MPIManager;
class BeadFactory;
/** XML Parser for PIMC simulation data.
  * @version $Revision$
  * @todo Get rid of explicit use of blitz.
  * @todo Clean up parsing of ConditionalDensityGrid.
  * @author John Shumway */
class PIMCParser : public XMLUnitParser {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  /// Constructor.
  PIMCParser(const SimulationInfo&, Action*, DoubleAction*, ActionChoiceBase*,
      EstimatorManager*, const BeadFactory&, MPIManager* mpi=0);
  /// Virtual destructor.
  virtual ~PIMCParser();
  /// Parse some xml.
  void parse(const xmlXPathContextPtr& ctxt);
  /// Return the Algorithm object.
  Algorithm* getAlgorithm() {return algorithm;}
private:
  /// Holder for the Paths object.
  Paths* paths;
  /// Holder for the Algorithm tree.
  Algorithm* algorithm;
  /// Holder for current the SectionChooser object for parsing.
  SectionChooser* sectionChooser;
  /// Holder for current the DoubleSectionChooser object for parsing.
  DoubleSectionChooser* doubleSectionChooser;
  /// General simulation information.
  const SimulationInfo& simInfo;
  /// Parse an algorithm.
  Algorithm* parseAlgorithm(const xmlXPathContextPtr& ctxt);
  /// Parse the body of an algorithm.
  void parseBody(const xmlXPathContextPtr& ctxt, CompositeAlgorithm*);
  /// Parse the boundary.
  void parseBoundary(const xmlXPathContextPtr& ctxt, Vec& min, Vec& max);
  /// The Action.
  Action* action;
  /// The DoubleAction.
  DoubleAction* doubleAction;
  /// Pointer to ActionChoice.
  ActionChoiceBase* actionChoice;
  /// The estimator manager.
  EstimatorManager* estimators;
  /// The probability density grid.
  ProbDensityGrid *probDensityGrid;
  std::vector<ProbDensityGrid*> condDensityGrid;
  /// The BeadFactory.
  const BeadFactory &beadFactory;
  /// The MPIManager.
  MPIManager *mpi;
  /// Current number of levels in sectionChooser or doubleSectionChooser.
  int nlevel;
  /// Count number of loop iterations surrounding an XML node.
  int getLoopCount(const xmlXPathContextPtr& ctxt);
};
#endif
