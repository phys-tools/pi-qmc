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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "MainSpinParser.h"
#include "SimulationInfo.h"
#include "Beads.h"
#include "BeadFactory.h"
#include <vector>
#include <cstdlib>
#include <blitz/tinyvec-et.h>

MainSpinParser::~MainSpinParser() {
//  delete simInfo;
}

void MainSpinParser::parse(const xmlXPathContextPtr& ctxt) {
  // Parse the species information.
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Spin",ctxt);

  if (obj->nodesetval->nodeNr>0) {
    xmlNodePtr node=obj->nodesetval->nodeTab[0]; ctxt->node=node;
    std::cout << "Using spin" << std::endl;

    simInfo.getBeadFactory().addAuxBeads(new Beads<4>(1,1),"spin");
    simInfo.spinOmega=getDoubleAttribute(node,"omega");
    if (simInfo.spinOmega==0) simInfo.spinOmega=1.0;
    std::cout << "Spin omega = " << simInfo.spinOmega << std::endl;
  } else {
    std::cout << "Not using spin" << std::endl;
  }
}
