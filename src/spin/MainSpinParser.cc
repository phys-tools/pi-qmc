#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "MainSpinParser.h"
#include "base/Beads.h"
#include "base/BeadFactory.h"
#include "base/SimulationInfo.h"
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
