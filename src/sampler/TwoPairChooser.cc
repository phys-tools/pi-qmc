// $Id: TwoPairChooser.cc 338 2010-11-30 18:56:16Z john.shumwayjr $
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
#include "TwoPairChooser.h"
#include "Species.h"
#include "RandomNumGenerator.h"
#include "Permutation.h"
#include <iostream>

TwoPairChooser::TwoPairChooser(const Species &s1, const Species &s2)
   : ParticleChooser(4), PermutationChooser(4),
     npart1(s1.count), npart2(s2.count),
     ifirst1(s1.ifirst), ifirst2(s2.ifirst) {
  (*this->permutation)[0]=1;
  (*this->permutation)[1]=0;
  (*this->permutation)[2]=3;
  (*this->permutation)[3]=2;
}

TwoPairChooser::~TwoPairChooser() {}

void TwoPairChooser::chooseParticles() {
  int i=0,j=0,k=0,l=0;
  do {
    i = (int)(npart1*RandomNumGenerator::getRand())+ifirst1;
  } while (i==npart1+ifirst1);
  do {
    j = (int)(npart1*RandomNumGenerator::getRand())+ifirst1;
  } while (j==i || j==npart1+ifirst1);
  do {
    k = (int)(npart2*RandomNumGenerator::getRand())+ifirst2;
  } while (k==i || k==j || k==npart2+ifirst2);
  do {
    l = (int)(npart2*RandomNumGenerator::getRand())+ifirst2;
  } while (l==k || l==i || l==j || l==npart2+ifirst2);
  index(0)=i;
  index(1)=j;
  index(2)=k;
  index(3)=l;
}
