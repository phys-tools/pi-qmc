// $Id$
/*  Copyright (C) 2007-8 John B. Shumway, Jr.

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
#include "StillWebAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "sampler/MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"
#include "SimulationInfo.h"

StillWebAction::StillWebAction(const SimulationInfo &simInfo, 
  const std::string &filename) : 
    tau(simInfo.getTau()),
    npart(simInfo.getNPart()),
    A(7.049556277),
    B(0.6022245584), p(4), q(0), a(1.80), gamma(1.20), param(3) {
  // Si-Si:
  param[0].lambda=21.0;
  param[0].sigmainv=1./((2.347772155/0.529177)/pow(2,1./6));
  param[0].epsilon=(2.17/27.211);
  // Si-Ge, Ge-Si:
  param[1].lambda=sqrt(21.0*31.0);
  param[1].sigmainv=1./((2.396885/0.529177)/pow(2,1./6));
  param[1].epsilon=(2.0427/27.211);
  // Ge-Ge
  param[2].lambda=31.0;
  param[2].sigmainv=1./((2.44598/0.529177)/pow(2,1./6));
  param[2].epsilon=(1.93/27.211);
}

double StillWebAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  double taueff = tau*nStride;
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      for (int jpart=0; jpart<npart; ++jpart) {
        if (jpart==i) continue;
        Vec deltaj=movingBeads(iMoving,islice)-sectionBeads(jpart,islice);
        cell.pbc(deltaj);
        double rij=sqrt(dot(deltaj,deltaj));
        deltaAction+=f2(rij,0,0)*taueff;
        for (int kpart=0; kpart<jpart; ++kpart) {
          if (kpart==i) continue;
          // Add action for moving beads.
          Vec deltak=movingBeads(iMoving,islice)-sectionBeads(kpart,islice);
          cell.pbc(deltak);
          double rik=sqrt(dot(deltak,deltak));
          double costheta=dot(deltaj,deltak)/(rij*rik);
          deltaAction+=h(rij,rik,costheta,0,0,0)*taueff;
        }
        deltaj=sectionBeads(i,islice)-sectionBeads(jpart,islice);
        cell.pbc(deltaj);
        rij=sqrt(dot(deltaj,deltaj));
        deltaAction-=f2(rij,0,0)*taueff;
        for (int kpart=0; kpart<jpart; ++kpart) {
          if (kpart==i) continue;
          // Subtract action for old beads.
          Vec deltak=sectionBeads(i,islice)-sectionBeads(kpart,islice);
          cell.pbc(deltak);
          double rik=sqrt(dot(deltak,deltak));
          double costheta=dot(deltaj,deltak)/(rij*rik);
          deltaAction-=h(rij,rik,costheta,0,0,0)*taueff;
        }
      }
    }
  }
  return deltaAction;
}

double StillWebAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void StillWebAction::getBeadAction(const Paths& paths, int ipart, int islice,
     double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  Vec ri=paths(ipart,islice);
  const SuperCell &cell = paths.getSuperCell();
  utau=0.;
  for (int jpart=0; jpart<npart; ++jpart) {
    if (jpart==ipart) continue;
    Vec deltaj=ri-paths(jpart,islice);
    cell.pbc(deltaj);
    double rij=sqrt(dot(deltaj,deltaj));
    utau+=0.5*f2(rij,0,0);
    for (int kpart=0; kpart<jpart; ++kpart) {
      if (kpart==ipart) continue;
      Vec deltak=ri-paths(kpart,islice);
      cell.pbc(deltak);
      double rik=sqrt(dot(deltak,deltak));
      double costheta=dot(deltaj,deltak)/(rij*rik);
      utau+=h(rij,rik,costheta,0,0,0);
    }
  }
  u=utau*tau;
}


double StillWebAction::h(const double r1, const double r2, 
                   const double costheta,
                   const int i, const int j, const int k ) const {
  const SWParam& p1=param[0];//ind(i,j)];
  const SWParam& p2=param[0];//ind(i,k)];
  const double onethird=1./3.;
  double y1=r1*p1.sigmainv;
  double y2=r2*p2.sigmainv;
  return (y1<a && y2<a)
          ? sqrt(p1.epsilon*p2.epsilon*p1.lambda*p2.lambda)
//            *exp(gamma/(y1-a)+gamma/(y2-a))*pow((costheta+onethird),2)
            *exp(gamma/(y1-a)+gamma/(y2-a))
            *(costheta+onethird)*(costheta+onethird)
          : 0;
}

double StillWebAction::f2(const double r, const int i, const int j) const {
  const SWParam& p1=param[0];//param[ind(i,j)];
  double y=r*p1.sigmainv;
  return (y<a) ? p1.epsilon*A*(B*pow(y,-p)-1)*exp(1/(y-a)) : 0;
  //return (y<a) ? p1.epsilon*(A*(B/(y*y*y*y)-1)*exp(1/(y-a))+1) : p1.epsilon;
}

