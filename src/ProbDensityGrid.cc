// $Id: ProbDensityGrid.cc,v 1.8 2006/10/18 17:08:19 jshumwa Exp $
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
#include "ProbDensityGrid.h"
#include "SimulationInfo.h"
#include "Paths.h"

ProbDensityGrid::ProbDensityGrid(const IVec n,
  const double a, const SimulationInfo& simInfo, const Paths* paths)
  : grid(simInfo.getNSpecies()), paths(paths),
    norm(0), n(n), a(a), b(3), index(simInfo.getNPart(),0) {
  for (int i=0; i<NDIM; ++i) {
    b[i]=0.0;
    b[i][i]=1/a;
  }
  for (int ispec=0; ispec<simInfo.getNSpecies(); ++ispec) {
    grid[ispec]=new LIArray(n);
    (*grid[ispec])=0;
    for (int ipart=0; ipart<simInfo.getNPart(); ++ipart) {
      if (&simInfo.getPartSpecies(ipart) == &simInfo.getSpecies(ispec)) {
        index[ipart]=ispec; 
      }
    }
  }
}

ProbDensityGrid::~ProbDensityGrid() {
  for (unsigned int i=0; i<grid.size(); ++i) delete grid[i];
}

void ProbDensityGrid::bin() {
  //Use callback for loop over all links.
  paths->sumOverLinks(*this);
}

void ProbDensityGrid::handleLink(const Vec& start, const Vec& end,
  const int ipart, const int islice, const Paths& paths) {
  IVec i;
  for (int idim=0; idim<NDIM; ++idim) {
    i[idim]=(int)floor(dot(start,b[idim])+n[idim]/2);
    if (i[idim]<0||i[idim]>=n[idim]) return;
  }
  ++((*grid[index[ipart]])(i));
}

void ProbDensityGrid::endCalc(const int nslice) {norm+=nslice;}
