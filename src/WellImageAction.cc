// $Id$
/*  Copyright (C) 2008 John B. Shumway, Jr.

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
#include "WellImageAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "sampler/MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include <fstream>

WellImageAction::WellImageAction(const SimulationInfo& simInfo,
  double epsIn, double epsOut, double width, double z0, double delta) 
  : tau(simInfo.getTau()), npart(simInfo.getNPart()),
    epsIn(epsIn), epsOut(epsOut), width(width), z0(z0), delta(delta),
    p1((epsIn-epsOut)/(epsIn+epsOut)), p2(2.*epsIn/(epsIn+epsOut)),
    q(npart), cell(*simInfo.getSuperCell()) {
  q=1.0;
}

double WellImageAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  double deltaAction=0;
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double tauStride = tau*nStride;
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Add action for moving beads.
      Vec pos=movingBeads(iMoving,islice);
      pos[NDIM-1] -= z0;  
      cell.pbc(pos);
      double z = fabs(pos[NDIM-1]);
      if (z<0.5*width) {
        deltaAction+=vSelfC(z);
      } else {
        deltaAction+=vSelfR(z);
      }
      // Subtract action for old beads.
      pos=sectionBeads(i,islice);
      pos[NDIM-1] -= z0;  
      cell.pbc(pos);
      z = fabs(pos[NDIM-1]);
      if (z<0.5*width) {
        deltaAction-=vSelfC(z);
      } else {
        deltaAction-=vSelfR(z);
      }
    }
  }
  for (int iMoving=0; iMoving<nMoving; ++iMoving) {
    const int i=index(iMoving);
    continue; //Particle not in this interaction.
    for (int j=0; j<npart; ++j) {
      bool isMoving=false; int jMoving=0;
      for (int k=0;k<nMoving;++k) {
        if (j==index(k)) {isMoving=true; jMoving=k; break;}
      }
      if (isMoving && i<=j) continue; //Don't double count moving interactions.
      for (int islice=nStride; islice<nSlice; islice+=nStride) {
        // Add action for moving beads.
        Vec pos1=movingBeads(iMoving,islice);
        Vec pos2=(isMoving)?movingBeads(jMoving,islice):sectionBeads(j,islice);
        pos1[NDIM-1]-=z0; cell.pbc(pos1);
        pos2[NDIM-1]-=z0; cell.pbc(pos2);
        Vec delta=pos1-pos2; cell.pbc(delta);
        double z1=pos1[NDIM-1];
        double z2=pos2[NDIM-1];
        double rho2=0;
        for (int idim=0; idim<NDIM-1; ++idim) rho2+=delta[idim]*delta[idim];
       	if (fabs(z1)<0.5*width && fabs(z2)<0.5*width) {
          deltaAction+=vCC(z1,z2,rho2);
       	} else if (fabs(z1)<0.5*width) {
          if (z2>0.) {
            deltaAction+=vRC(z2,z1,rho2);
          } else {
            deltaAction+=vRC(-z2,-z1,rho2);
          }
       	} else if (fabs(z2)<0.5*width) {
          if (z1>0.) {
            deltaAction+=vRC(z1,z2,rho2);
          } else {
            deltaAction+=vRC(-z1,-z2,rho2);
          }
        } else if (z1*z2>0) {
          deltaAction+=vRR(fabs(z1),fabs(z2),rho2);
        } else {
          deltaAction+=vLR(-fabs(z1),fabs(z2),rho2);
        } 
        // Subtract action for old beads.
        pos1=sectionBeads(i,islice);
        pos2=sectionBeads(j,islice);
        pos1[NDIM-1]-=z0; cell.pbc(pos1);
        pos2[NDIM-1]-=z0; cell.pbc(pos2);
        delta=pos1-pos2; cell.pbc(delta);
        z1=pos1[NDIM-1];
        z2=pos2[NDIM-1];
        rho2=0;
        for (int idim=0; idim<NDIM-1; ++idim) rho2+=delta[idim]*delta[idim];
       	if (fabs(z1)<0.5*width && fabs(z2)<0.5*width) {
          deltaAction-=vCC(z1,z2,rho2);
       	} else if (fabs(z1)<0.5*width) {
          if (z2>0.) {
            deltaAction-=vRC(z2,z1,rho2);
          } else {
            deltaAction-=vRC(-z2,-z1,rho2);
          }
       	} else if (fabs(z2)<0.5*width) {
          if (z1>0.) {
            deltaAction-=vRC(z1,z2,rho2);
          } else {
            deltaAction-=vRC(-z1,-z2,rho2);
          }
        } else if (z1*z2>0) {
          deltaAction-=vRR(fabs(z1),fabs(z2),rho2);
        } else {
          deltaAction-=vLR(-fabs(z1),fabs(z2),rho2);
        }
      }
    }
  }
  return tauStride*deltaAction;
}

double WellImageAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void WellImageAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  Vec pos1=paths(ipart,islice);
  pos1[NDIM-1] -= z0;  
  cell.pbc(pos1);
  double z1 = pos1[NDIM-1];
  if (fabs(z1)<0.5*width) {
    utau+=vSelfC(z1);
  } else {
    utau+=vSelfR(fabs(z1));
  }
  for (int j=0; j<npart; ++j) {
    if (j!=ipart) {
      Vec pos2=paths(j,islice);
      pos2[NDIM-1] -= z0;  
      cell.pbc(pos2);
      double z2 = pos2[NDIM-1];
      Vec delta=pos1-pos2; cell.pbc(delta);
      double rho2=0;
      for (int idim=0; idim<NDIM-1; ++idim) rho2+=delta[idim]*delta[idim];
      if (fabs(z1)<0.5*width && fabs(z2)<0.5*width) {
        utau+=0.5*vCC(z1,z2,rho2);
      } else if (fabs(z1)<0.5*width) {
        if (z2>0.) {
          utau+=0.5*vRC(z2,z1,rho2);
        } else {
          utau+=0.5*vRC(-z2,-z1,rho2);
        }
      } else if (fabs(z2)<0.5*width) {
        if (z1>0.) {
          utau+=0.5*vRC(z1,z2,rho2);
        } else {
          utau+=0.5*vRC(-z1,-z2,rho2);
        }
      } else if (z1*z2>0) {
        utau+=0.5*vRR(fabs(z1),fabs(z2),rho2);
      } else {
        utau+=0.5*vLR(-fabs(z1),fabs(z2),rho2);
      }
    }
  }
  u=utau*tau;
}

inline double WellImageAction::vSelfC(double z) const {
  double sum=p1/(2.*epsIn*fabs(2.*z+width+2.*delta))
            +p1/(2.*epsIn*fabs(2.*z-width-2.*delta));
  double pn=p1;
  int sign=-1;
  for (int n=2; n<4; ++n) {
    pn *= p1;
    sign *= -1;
    sum += pn/(2.*epsIn*fabs(z-sign*z+n*width));
    sum += pn/(2.*epsIn*fabs(z-sign*z-n*width));
  } 
  return sum;
}

inline double WellImageAction::vSelfR(double z) const {
  double sum = -p1/(2*epsOut*fabs(2.*z-width+2.*delta));
  double pn=p1;
  for (int n=1; n<6; n+=2) {
    sum += p2*pn/((epsIn+epsOut)*fabs(2.*z+n*width));
    pn *= p1*p1;
  }
  return sum;
}

inline double WellImageAction::vRC(double z1, double z2, double rho2) const {
  double sum=-1./(epsOut*sqrt((z1-z2)*(z1-z2)+rho2));
  double pn=1.;
  int sign=1;
  for (int n=0; n<4; ++n) {
    sum += p2*pn/(epsIn
                 *sqrt((z1-sign*z2+n*width)*(z1-sign*z2+n*width)+rho2));
    pn *= p1;
    sign *= -1;
  }
  return sum;
}

inline double WellImageAction::vRR(double z1, double z2, double rho2) const {
  double sum=-p1/(epsOut*sqrt((z1+z2-width)*(z1+z2-width)+rho2));
  double pn=p1;
  for (int n=1; n<4; ++n) {
    sum += p2*p2*pn/(epsIn
                 *sqrt((z1+z2+n*width)*(z1+z2+n*width)+rho2));
    pn *= p1*p1;
  }
  return sum;
}

inline double WellImageAction::vLR(double z1, double z2, double rho2) const {
  double sum=-1./(epsOut*sqrt((z1-z2)*(z1-z2)+rho2));
  double p2n=1.;
  for (int n=0; n<2; ++n) {
    sum += p2*p2*p2n/(epsIn*sqrt((z1-z2-n*width)*(z1-z2-n*width)+rho2));
    p2n *= p1*p1;
  }
  return sum;
}

inline double WellImageAction::vCC(double z1, double z2, double rho2) const {
  double sum=(1./epsIn-1./epsOut)*1./(sqrt((z1-z2)*(z1-z2)+rho2));
  double pn=p1;
  int sign=-1;
  for (int n=1; n<4; ++n) {
    sum += pn/(epsIn
                 *sqrt((z1-sign*z2+n*width)*(z1-sign*z2+n*width)+rho2));
    sum += pn/(epsIn
                 *sqrt((z1-sign*z2-n*width)*(z1-sign*z2-n*width)+rho2));
    pn *= p1;
    sign *= -1;
  }
  return sum;
}
