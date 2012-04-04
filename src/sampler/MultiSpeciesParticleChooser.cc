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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MultiSpeciesParticleChooser.h"
#include "util/RandomNumGenerator.h"
#include "SimulationInfo.h"
#include <iostream>

MultiSpeciesParticleChooser::MultiSpeciesParticleChooser( const Species *speciesList, const int nspecies,
  const int nmoving) : ParticleChooser(nmoving) {
  
  count = 0;
  for (int ispec=0; ispec<nspecies; ispec++){
    count +=  speciesList[ispec].count;
  }
 
  speciesContainer.resize(count);int k=0;
  for (int ispec=0; ispec<nspecies; ispec++){
    for (int i=0; i< speciesList[ispec].count; i++){
      speciesContainer(k++) = speciesList[ispec].ifirst + i;
    }
  }
}

MultiSpeciesParticleChooser::~MultiSpeciesParticleChooser() {
}

void MultiSpeciesParticleChooser::chooseParticles() {
  for (int imoving=0; imoving<index.size();imoving++) {
    int i = (int)(count*RandomNumGenerator::getRand());
    index(imoving) = speciesContainer(i) ;
  }
}
