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
#include "BondLengthEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Paths.h"
#include "stats/Units.h"
#include <blitz/tinyvec-et.h>
#include "SuperCell.h"

BondLengthEstimator::BondLengthEstimator(const SimulationInfo& simInfo,
  const Species& species1, const Species& species2, const std::string& unitName)
  : ScalarEstimator("bond_length_"+species1.name+"-"+species2.name,
      unitName,1./simInfo.getUnits()->getLengthScaleIn(unitName),0.),
    length(0), norm1(0), norm2(0), avlength(0), names(simInfo.getNPart()),
    spec1(species1.name), spec2(species2.name), cell(*simInfo.getSuperCell()) {
  for (int i=0; i<names.size(); ++i)
    names(i)=simInfo.getPartSpecies(i).name;
  std::cout << "Bond Length Estimator "
            << spec1 << " " << spec2 << std::endl;
}

void BondLengthEstimator::initCalc(const int nslice, const int firstSlice) {
  length=norm1=0;
}

void BondLengthEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  if(names(ipart)==spec1) {
    for (int jpart=0; jpart<names.size(); ++jpart) {
      if(names(jpart)==spec2 && jpart!=ipart) {
        Vec delta=end-paths(jpart,islice);
        cell.pbc(delta);
        length+=sqrt(dot(delta,delta));
        norm1++;
      }
    }
  }
}

void BondLengthEstimator::endCalc(const int nslice) {
   avlength+=(length/norm1);
   norm2++;
}
