// $Id:$
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
#include "AlternatingParticleChooser.h"
#include "util/RandomNumGenerator.h"
#include "SimulationInfo.h"
#include <iostream>

AlternatingParticleChooser::AlternatingParticleChooser(const Species& species1,
		 const Species& species2, const int nmoving)
                 : ParticleChooser(nmoving), npart1(species1.count),
		   npart2(species2.count), ifirst1(species1.ifirst),
		   ifirst2(species2.ifirst)
{
}

AlternatingParticleChooser::~AlternatingParticleChooser() {}

void AlternatingParticleChooser::chooseParticles() {
  int i = (int)(npart1 * RandomNumGenerator::getRand() + ifirst1);
  if (i==ifirst1+npart1) i=ifirst1+npart1-1;
  index(0) = i;
  i = (int)(npart2 * RandomNumGenerator::getRand() + ifirst2);
  if (i==ifirst2+npart2) i=ifirst2+npart2-1;
  index(1) = i;
  if (RandomNumGenerator::getRand()<0.5) {
    int temp = index(0);
    index(0) = index(1);
    index(1) = temp;
  }
//  std::cout<<"index(0) = "<<index(0)<<", index(1) = "<<index(1)<<std::endl;
//  std::cin.ignore();
}

