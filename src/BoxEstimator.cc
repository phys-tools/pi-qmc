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
#include "BoxEstimator.h"

BoxEstimator::BoxEstimator(const SimulationInfo& simInfo, const Species& species,
                           const double lower, const double upper, const int index)
  : ScalarEstimator(index==0 ? "box_x_"+species.name : 
                    index==1 ? "box_y_"+species.name : "box_z_"+species.name),
    in(0.), intot(0.), innorm(0), low(lower), up(upper), index(index),
    names(simInfo.getNPart()), spec(species.name) {
    for (int i=0; i<names.size(); ++i)
      names(i)=simInfo.getPartSpecies(i).name;
    std::cout << "Box Estimator " << spec << " bounds: "
              << (index==0 ? "x = " : index==1 ? "y = " : "z = ")
              << "[" << low << ", " << up <<"]" << std::endl;
}

void BoxEstimator::initCalc(const int nslice, const int firstSlice) {
  in=0.;
}

void BoxEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  if(names(ipart)==spec)
    if(end(index)>=low && end(index)<=up)
      in++;
}

void BoxEstimator::endCalc(const int nslice) {
  in/=nslice;
  intot+=in; innorm+=1;
}
