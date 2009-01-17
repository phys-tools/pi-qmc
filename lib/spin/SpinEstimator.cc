// $Id: SpinEstimator.cc,v 1.4 2007/10/03 12:53:56 jshumwa Exp $
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
#include "SpinEstimator.h"
#include "SimulationInfo.h"
#include "Action.h"
#include "DoubleAction.h"
#include <blitz/tinyvec.h>

SpinEstimator::SpinEstimator(
  const SimulationInfo& simInfo, const int idim, const double gc)
  : ScalarEstimator("spin"+((idim==0)?namex:(idim==1)?namey:namez)),
    idim(idim), gc(gc), invgc(pow((1.0+gc),3)), energy(0), etot(0), enorm(0),
    l2(1/(simInfo.getSpecies(0).mass*simInfo.spinOmega)) {
}

void SpinEstimator::initCalc(const int nslice, const int firstSlice) {
  energy=0;
}

void SpinEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
#if NDIM==4 
  const Vec s=end;
#else
  const blitz::TinyVector<double,4> s
            =*reinterpret_cast<const blitz::TinyVector<double,4>*>(
               paths.getAuxBead(ipart,islice,1));
#endif
	double dots=dot(s,s)/l2;
//  energy+=dots;
	double f=1.5*invgc*exp(-gc*dots)/dots;
  switch (idim) {
    case 0:
      energy+=2*f*(-s[0]*s[2]+s[1]*s[3])/l2;
      break;
    case 1:
      energy+=2*f*(s[0]*s[1]+s[2]*s[3])/l2;
      break;
    case 2:
      energy+=f*(s[0]*s[0]-s[1]*s[1]-s[2]*s[2]+s[3]*s[3])/l2;
      break;
  }
}

void SpinEstimator::endCalc(const int nslice) {
  energy/=nslice;
  etot+=energy; enorm+=1;
}

const std::string  SpinEstimator::namex("_x"); 
const std::string  SpinEstimator::namey("_y"); 
const std::string  SpinEstimator::namez("_z"); 
