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
#include "SimInfoParser.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include <vector>
#include <blitz/tinyvec-et.h>
#include "stats/Units.h"

SimInfoParser::~SimInfoParser() {
  delete simInfo;
}

void SimInfoParser::parse(const xmlXPathContextPtr& ctxt) {
  if (!simInfo) simInfo = new SimulationInfo;
  // Parse information about Units (hard-coded for atomic units for now).
  //xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Units",ctxt);
  //xmlNodePtr node=obj->nodesetval->nodeTab[0]; ctxt->node=node;
  //xmlXPathFreeObject(obj);
  //const std::string eunit=getStringAttribute(node,"energyUnit");
  //const std::string lunit=getStringAttribute(node,"lengthUnit");
  simInfo->units=new Units("Ha","a0");
  this->units=simInfo->units;
  // Parse the species information.
  simInfo->npart=0;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Species",ctxt);
  int nspecies=obj->nodesetval->nodeNr;
  std::vector<Species*>& speciesList(simInfo->speciesList);
  speciesList.resize(nspecies);
  for (int ispecies=0; ispecies<nspecies; ++ispecies) {
    xmlNodePtr specNode=obj->nodesetval->nodeTab[ispecies];    
    Species *species=new Species();
    speciesList[ispecies]=species;
    species->name=getStringAttribute(specNode,"name");
    species->count=getIntAttribute(specNode,"count");
    species->mass=getMassAttribute(specNode,"mass");
    species->charge=getDoubleAttribute(specNode,"charge");
    species->twos=getIntAttribute(specNode,"twos");
    std::string type=getStringAttribute(specNode,"type");
    species->isFermion=(type=="fermion")?true:false;
    species->isStatic=getBoolAttribute(specNode,"isStatic");
    simInfo->npart+=species->count;
    // Look for optional child nodes.
    ctxt->node=specNode;
    xmlXPathObjectPtr obj2 = xmlXPathEval(BAD_CAST"AnisotropicMass",ctxt);
    if (obj2->nodesetval->nodeNr>0) {
      xmlNodePtr node=obj2->nodesetval->nodeTab[0];    
      species->anMass = new Vec(getVecAttribute(node,""));
    }
    xmlXPathFreeObject(obj2);
  }
  simInfo->speciesIndex.resize(simInfo->npart);
  int ipart=0;
  for (int ispecies=0; ispecies<nspecies; ++ispecies) {
    Species *species=speciesList[ispecies];
    species->ifirst=ipart;
    for (int i=0; i<species->count; ++i) {
      simInfo->speciesIndex[ipart++]=species;
    }
  }
  xmlXPathFreeObject(obj);
  // Parse the SuperCell information.
  obj = xmlXPathEval(BAD_CAST"//SuperCell",ctxt);
  xmlNodePtr node=obj->nodesetval->nodeTab[0]; ctxt->node=node;
  xmlXPathFreeObject(obj);
  const double a=getLengthAttribute(node,"a");
  Vec extent = getVecAttribute(node,"");
  simInfo->superCell=new SuperCell(a*extent);
  simInfo->superCell->computeRecipricalVectors();
  std::cout << "Supercell dimensions: "
            << simInfo->superCell->a << std::endl;
  // Parse the simulation temperature.
  obj = xmlXPathEval(BAD_CAST"//Temperature",ctxt);
  simInfo->temperature=getEnergyAttribute(obj->nodesetval->nodeTab[0],"value");
  simInfo->nslice=getIntAttribute(obj->nodesetval->nodeTab[0],"nslice");
  if (simInfo->nslice>0) {
    simInfo->tau=1.0/(simInfo->temperature * simInfo->nslice);
  } else {
    simInfo->tau=getTimeAttribute(obj->nodesetval->nodeTab[0],"tau");
    simInfo->nslice=(int)(1.0/(simInfo->temperature*simInfo->tau)+0.1);
  }
  xmlXPathFreeObject(obj);
}
