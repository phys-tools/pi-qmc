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
#ifndef __ActionParser_h_
#define __ActionParser_h_
#include "XMLUnitParser.h"
class Action;
class DoubleAction;
class ActionChoiceBase;
class SimulationInfo;
class MPIManager;
class CompositeAction;
class CompositeDoubleAction;
#include "AugmentedNodes.h"
#include <blitz/array.h>
#include <vector>
/** XML Parser for actions.
  * @version $Revision$
  * @author John Shumway */
class ActionParser : public XMLUnitParser {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Constructor.
  ActionParser(const SimulationInfo&, const int maxlevel, MPIManager*);
  /// Virtual destructor.
  ~ActionParser(){};
  /// Parse some xml.
  void parse(const xmlXPathContextPtr& ctxt);
  /// Return the Action object.
  Action* getAction() {return action;}
  /// Return the DoubleAction object.
  DoubleAction* getDoubleAction() {return doubleAction;}
  /// Return the ActionChoice object.
  ActionChoiceBase* getActionChoice() {return actionChoice;}
private:
  /// Function to parse all the action xml tags.
  void parseActions(const xmlXPathContextPtr& ctxt, xmlXPathObjectPtr& obj,
    CompositeAction* action, CompositeDoubleAction* doubleComposite);
  /// The Action object.
  Action* action;
  /// The DoubleAction object.
  DoubleAction* doubleAction;
  /// Choice of actions if different models are being sampled.
  ActionChoiceBase* actionChoice;
  /// General simulation information.
  const SimulationInfo& simInfo;
  /// The timestep.
  const double tau;
  /// The maximum number of levels.
  const int maxlevel;
  /// Parse actions for EwaldAction.
  Action* parseEwaldActions(const xmlXPathContextPtr& ctxt);
  /// Parse atomic orbitals for AugmentedNdoes.
  void parseOrbitalDM(std::vector<const AugmentedNodes::AtomicOrbitalDM*>&,
                      const xmlXPathContextPtr& ctxt);
  /// Letters associated with directions in input file.
  static const std::string dimName;
  /// MPIManager.
  MPIManager *mpi;
};
#endif
