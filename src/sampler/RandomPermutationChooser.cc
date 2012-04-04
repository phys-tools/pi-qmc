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
#include "RandomPermutationChooser.h"
#include "util/RandomNumGenerator.h"
#include "util/Permutation.h"

RandomPermutationChooser::RandomPermutationChooser(const int nsize)
  : PermutationChooser(nsize), nsize(nsize) {
}

bool RandomPermutationChooser::choosePermutation() {
  for (int ifrom=0; ifrom<nsize;) {
    int ito = (int)(nsize*RandomNumGenerator::getRand());
    if (ito==nsize) ito=nsize-1;
    for (int jfrom=0; jfrom<=ifrom; ++jfrom) {
      if (jfrom==ifrom) (*permutation)[ifrom++]=ito;
      if ((*permutation)[jfrom]==ito) break;
    }
  }
  return true;
}
