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
#include "QuadrupoleMomentEstimator.h"

QuadrupoleMomentEstimator::QuadrupoleMomentEstimator(const SimulationInfo& simInfo, const int index)
  : ScalarEstimator(index==0 ? "quadrupole_moment_xx" : index==1 ? "quadrupole_moment_xy" :
                    index==2 ? "quadrupole_moment_xz" : index==3 ? "quadrupole_moment_yy" :
                    index==4 ? "quadrupole_moment_yz" : index==5 ? "quadrupole_moment_zz" : ""),
    quadrupole(0), Qtot(0), Qnorm(0), q(simInfo.getNPart()), index(index) {
    for (int i=0; i<q.size(); ++i) q(i)=simInfo.getPartSpecies(i).charge;
    std::cout << "Quadrupole Moment Estimator component = "
              << (index==0 ? "xx" : index==1 ? "xy" : index==2 ? "xz" : 
                  index==3 ? "yy" : index==4 ? "yz" : "zz") << std::endl;
}

void QuadrupoleMomentEstimator::initCalc(const int nslice, const int firstSlice) {
  quadrupole=0;
}

void QuadrupoleMomentEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  switch(index)
  {
    case 0: //Q11
      quadrupole+=q(ipart)*(2.*end(0)*end(0)-end(1)*end(1)-end(2)*end(2));
    break;
    case 1: //Q12
      quadrupole+=q(ipart)*3.*end(0)*end(1);
    break;
    case 2: //Q13
      quadrupole+=q(ipart)*3.*end(0)*end(2);
    break;
    case 3: //Q22
      quadrupole+=q(ipart)*(2.*end(1)*end(1)-end(0)*end(0)-end(2)*end(2));
    break;
    case 4: //Q23
      quadrupole+=q(ipart)*3.*end(1)*end(2);
    break;
    default:
    case 5: //Q33
      quadrupole+=q(ipart)*(2.*end(2)*end(2)-end(0)*end(0)-end(1)*end(1));
    break;
  }

}

void QuadrupoleMomentEstimator::endCalc(const int nslice) {
  quadrupole/=nslice;
  Qtot+=quadrupole; Qnorm+=1;
}
