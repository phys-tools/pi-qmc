//$Id: SpinChoiceFixedNodeAction.cc 383 2011-04-20 16:56:02Z john.shumwayjr $
/*  Copyright (C) 2004-2009, 2011 John B. Shumway, Jr.

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "SpinChoiceFixedNodeAction.h"
#include "SpinModelState.h"

SpinChoiceFixedNodeAction::SpinChoiceFixedNodeAction(
  const SimulationInfo &simInfo,
  const Species &species, NodeModel *nodeModel, bool withNodalAction,
  bool useDistDerivative, int maxlevel, bool useManyBodyDistance) 
  : FixedNodeAction(simInfo,species,nodeModel,withNodalAction,
      useDistDerivative,maxlevel,useManyBodyDistance) {
  std::cout << "nSpeciesPart" << nSpeciesPart << std::endl;
  modelState = new SpinModelState(nSpeciesPart);
}

SpinChoiceFixedNodeAction::~SpinChoiceFixedNodeAction() {
}
