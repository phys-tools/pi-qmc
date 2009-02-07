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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "SpeciesParticleChooser.h"
#include "RandomNumGenerator.h"
#include "SimulationInfo.h"
#include <iostream>

SpeciesParticleChooser::SpeciesParticleChooser(const Species& species,
  const int nmoving) 
  : ParticleChooser(nmoving), npart(species.count), ifirst(species.ifirst) {
  name=species.name;
}

SpeciesParticleChooser::~SpeciesParticleChooser() {}

void SpeciesParticleChooser::chooseParticles() {
  for (int imoving=0; imoving<index.size();) {
    int i = (int)(npart*RandomNumGenerator::getRand()+ifirst);
    if (i==npart+ifirst) i=npart+ifirst-1;
    for (int jmoving=0; jmoving<=imoving; ++jmoving) {
      if (jmoving==imoving) index(imoving++)=i;
      if (index(jmoving)==i) break;
    }
  }
}
