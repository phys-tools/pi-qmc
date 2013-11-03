#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ConditionalDensityGrid.h"
#include "base/Paths.h"
#include "base/Species.h"
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
