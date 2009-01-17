// $Id: ConditionalDensityGrid.cc,v 1.6 2006/10/18 17:08:18 jshumwa Exp $
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
#include "ConditionalDensityGrid.h"
#include "Paths.h"
#include "Species.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

ConditionalDensityGrid::ConditionalDensityGrid(const IVec n,
  const double a, const SimulationInfo& simInfo,
  const Paths* paths, const Species& species, Vec center, const double radius)
  : ProbDensityGrid(n,a,simInfo,paths), ifirst(species.ifirst),
    ilast(species.ifirst+species.count), center(center), radius(radius) {
}

void ConditionalDensityGrid::handleLink(const Vec& start, const Vec& end,
  const int ipart, const int islice, const Paths& paths) {
  // Bin the density for this slice if we meet the condition.
  if (ipart>=ifirst && ipart<ilast) { 
    Vec delta=start-center;
    if (dot(delta,delta)<radius*radius) {
      // Good, but make sure we are the first particle meeting the condition.
      for (int jpart=ifirst; jpart<ipart;++jpart) {
        Vec delta=paths(jpart,islice)-center;
        if (dot(delta,delta)<radius*radius) return;
      }
      // Now bin the density.
      for (int jpart=0; jpart<(int)index.size(); ++jpart) {
        if (ipart==jpart) continue;
        Vec point=paths(jpart,islice);
        IVec i;
        for (int idim=0; idim<NDIM; ++idim) {
          i[idim]=(int)floor(dot(point,b[idim])+n[idim]/2);
          if (i[idim]<0||i[idim]>=n[idim]) break;
          // Increment the bin if we made it through all dimensions.
          if (idim==NDIM-1) ++((*grid[index[jpart]])(i));
        }
      }
      ++norm; 
    }
  }
}

void ConditionalDensityGrid::endCalc(const int nslice) {;}
