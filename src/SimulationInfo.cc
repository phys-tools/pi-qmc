// $Id: SimulationInfo.cc,v 1.5 2006/10/18 17:08:19 jshumwa Exp $
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
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "BeadFactory.h"

SimulationInfo::SimulationInfo()
  : superCell(0),npart(0),speciesList(0),speciesIndex(0),temperature(0),
    tau(1), nslice(0), beadFactory(*new BeadFactory()) {
}

SimulationInfo::SimulationInfo(SuperCell* cell, int npart, 
    std::vector<Species*> &specList, std::vector<Species*> &specIndex ,
    double temperature, double tau, int nslice) 
  : superCell(cell),npart(npart),speciesList(specList), speciesIndex(specIndex),
    temperature(temperature), tau(tau), nslice(nslice),
    beadFactory(*new BeadFactory()) {
}

SimulationInfo::~SimulationInfo() {
  for (unsigned int i=0; i<speciesList.size(); ++i) delete speciesList[i];
  delete superCell;
  delete &beadFactory;
}

const Species& SimulationInfo::getSpecies(const std::string& name) const {
  for (unsigned int ispec=0; ispec<speciesList.size(); ++ispec) {
    if (speciesList[ispec]->name==name) return *speciesList[ispec];
  }
  return *speciesList[0];
}

std::ostream& operator<<(std::ostream& os, const SimulationInfo& si) {
  os << "Simulation Info:" << std::endl;
  for (unsigned int ispec=0; ispec<si.speciesList.size(); ++ispec) {
    os << "  " << *si.speciesList[ispec];
  }
//  os << "  " << *si.superCell;
  os << "  T=" << si.temperature << std::endl;
  return os;
}
