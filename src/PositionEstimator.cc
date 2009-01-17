//$Id:
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
#include "PositionEstimator.h"

PositionEstimator::PositionEstimator(const SimulationInfo& simInfo, const Species& species,
                                     const int index)
  : ScalarEstimator(index==0 ? "position_x_"+species.name : 
                    index==1 ? "position_y_"+species.name : "position_z_"+species.name),
    position(0), ptot(0), pnorm(0), index(index),
    names(simInfo.getNPart()), spec(species.name) {
    for (int i=0; i<names.size(); ++i)
      names(i)=simInfo.getPartSpecies(i).name;
    std::cout << "Position Estimator " << spec << " component = "
              << (index==0 ? "x" : index==1 ? "y" : "z") << std::endl;
}

void PositionEstimator::initCalc(const int nslice, const int firstSlice) {
  position=0;
}

void PositionEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  if(names(ipart)==spec)
    position+=end(index);
}

void PositionEstimator::endCalc(const int nslice) {
  position/=nslice;
  ptot+=position; pnorm+=1;
}
