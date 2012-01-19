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
#include "EwaldAction.h"
#include <cmath>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Paths.h"
#include "Beads.h"
#include "sampler/SectionSamplerInterface.h"
#include "sampler/DisplaceMoveSampler.h"

#if NDIM==3
EwaldAction::EwaldAction(const SimulationInfo& simInfo, 
    const double rcut, const double kcut, const int nreserve) 
  : cell(*simInfo.getSuperCell()), actions(0), rcut(rcut), kcut(kcut),
    kappa(2.5/rcut), tau(simInfo.getTau()), npart(simInfo.getNPart()),
    q(npart), pos(npart), selfEnergy(0),


    ikmax((int)(cell.a[0]*kcut/(2*PI)),(int)(cell.a[1]*kcut/(2*PI)),
          (int)(cell.a[2]*kcut/(2*PI))),
    deltak(2*PI*cell.b[0],2*PI*cell.b[1],2*PI*cell.b[2]),
    // allocate eikx(0,npart-1,0:kmax)
    eikx(blitz::shape(npart,ikmax[0]+1)),
    // allocate eiky(0,npart-1,-kmax:kmax)
    eiky(blitz::shape(0,-ikmax[1]), blitz::shape(npart,2*ikmax[1]+1)),
    // allocate eikz(0,npart-1,-kmax:kmax)
    eikz(blitz::shape(0,-ikmax[2]), blitz::shape(npart,2*ikmax[2]+1)),
    twoPiOverV(2*PI/(cell.a[0]*cell.a[1]*cell.a[2]))
{
  if (nreserve>0) actions.reserve(nreserve);
  for (int i=0; i<npart; ++i) q(i)=simInfo.getPartSpecies(i).charge;
  // Calculate self-energy contribution.
  for (int jpart=0; jpart<npart; ++jpart) selfEnergy+=q(jpart)*q(jpart);
  selfEnergy*=kappa/sqrt(PI);
  // Calculate k-vectors.
  totk=0;
  for (int kx=0; kx<=ikmax[0]; ++kx) {
    double kx2=kx*kx*deltak[0]*deltak[0];
    for (int ky=-ikmax[1]; ky<=ikmax[1]; ++ky) {
      double ky2=ky*ky*deltak[1]*deltak[1];
      for (int kz=-ikmax[2]; kz<=ikmax[2]; ++kz) {
        double k2=kx2+ky2+kz*kz*deltak[2]*deltak[2];
        if (k2<kcut*kcut && k2!=0) ++totk;
      }
    }
  }
  std::cout << "Ewald: totk=" << totk << std::endl;
  expoverk2.resize(totk);
  int ikvec=0;
  for (int kx=0; kx<=ikmax[0]; ++kx) {
    double kx2=kx*kx*deltak[0]*deltak[0];
    for (int ky=-ikmax[1]; ky<=ikmax[1]; ++ky) {
      double ky2=ky*ky*deltak[1]*deltak[1];
      for (int kz=-ikmax[2]; kz<=ikmax[2]; ++kz) {
        double k2=kx2+ky2+kz*kz*deltak[2]*deltak[2];
        if (k2<kcut*kcut && k2!=0) {
          expoverk2(ikvec++)=exp(-k2/(4*kappa*kappa))/k2;
        }
      }
    }
  }
}

EwaldAction::~EwaldAction() {
  for (SRActIter a=actions.begin(); a<actions.end(); ++a) delete *a;
}

double EwaldAction::getActionDifference(
    const SectionSamplerInterface& sampler, const int level) {
  double diff=0;
return diff;
  for (ConstSRActIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action) diff+=(*action)->getActionDifference(sampler,level);
  }
  if (level==0) { // Calculate long range action.
    const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
    const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
    const int nSlice=sectionBeads.getNSlice();
    const IArray& index=sampler.getMovingIndex(); 
    const int nMoving=index.size();
    for (int islice=1;islice<nSlice; ++islice) { 
      for (int i=0; i<npart; ++i) {
        pos(i)=sectionBeads(i,islice);
        cell.pbc(pos(i)-=sectionBeads(i,islice-1));
        pos(i)*=-0.5;
        pos(i)+=sectionBeads(i,islice);
      }
      double utauLR=calcLongRangeUtau(pos);
      diff-=utauLR*tau;
      for (int imoving=0; imoving<nMoving; ++imoving) {
        const int i=index(imoving);
        pos(i)=movingBeads(imoving,islice);
        cell.pbc(pos(i)-=movingBeads(imoving,islice-1));
        pos(i)*=-0.5;
        pos(i)+=movingBeads(imoving,islice);
      }
      utauLR=calcLongRangeUtau(pos);
      diff+=utauLR*tau;
    }
  }
  return diff;
}

//displace move
double EwaldAction::getActionDifference( const DisplaceMoveSampler& sampler, const int nMoving) {
  double diff=0;
return diff;
/*
  for (ConstSRActIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action) diff+=(*action)->getActionDifference(sampler, nMoving);
  }
  // if (level==0) { // Calculate long range action.
    const Beads<NDIM>& pathsBeads=sampler.getPathsBeads();
    const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
    const int nSlice=pathsBeads.getNSlice();
    const IArray& index=sampler.getMovingIndex(); 
    //  const int nMoving=index.size();
    for (int islice=1;islice<nSlice; ++islice) { 
      for (int i=0; i<npart; ++i) {
        pos(i)=pathsBeads(i,islice);
        cell.pbc(pos(i)-=pathsBeads(i,islice-1));
        pos(i)*=-0.5;
        pos(i)+=pathsBeads(i,islice);
      }
      double utauLR=calcLongRangeUtau(pos);
      diff-=utauLR*tau;
      for (int imoving=0; imoving<nMoving; ++imoving) {
        const int i=index(imoving);
        pos(i)=movingBeads(imoving,islice);
        cell.pbc(pos(i)-=movingBeads(imoving,islice-1));
        pos(i)*=-0.5;
        pos(i)+=movingBeads(imoving,islice);
      }
      utauLR=calcLongRangeUtau(pos);
      diff+=utauLR*tau;
    }
    //}
  return diff;
*/
}

double EwaldAction::getTotalAction(const Paths& paths, int level) const {
  double total=0;
  for (ConstSRActIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action) total+=(*action)->getTotalAction(paths,level);
  }
  return total;
}

void EwaldAction::getBeadAction(const Paths& paths, int ipart, int islice,
       double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  for (ConstSRActIter action=actions.begin(); action<actions.end(); ++action) {
    double ui=0, utaui=0, ulambdai=0; 
    Vec fmi,fpi;
    if (*action) (*action)->getBeadAction(paths,ipart,islice,
                                          ui,utaui,ulambdai,fmi,fpi);
    u+=ui; utau+=utaui; ulambda+=ulambdai; fm+=fmi; fp+=fpi;
  }
  if (ipart==0) {
    u-=selfEnergy*tau; utau-=selfEnergy;
    paths.getSlice(islice,pos);
    for (int i=0; i<npart; ++i) {
      Vec delta(paths.delta(i,islice,-1));
      pos(i)-=(cell.pbc(delta)*=0.5);
    }
    double utauLR=calcLongRangeUtau(pos);
    u+=utauLR*tau; utau+=utauLR;
  } 
}

void EwaldAction::initialize(const SectionChooser& sectionChooser) {
  for (ConstSRActIter action=actions.begin(); action<actions.end(); ++action) {
    (*action)->initialize(sectionChooser);
  }
}

void EwaldAction::acceptLastMove() {
  for (ConstSRActIter action=actions.begin(); action<actions.end(); ++action) {
    (*action)->acceptLastMove();
  }
}

void EwaldAction::setup() {
  for (SRActIter actptr=actions.begin(); actptr<actions.end(); ++actptr) {
    PairAction &action=**actptr;
    /// Remove long-range 1/r tail from PairAction.
    int n=action.ngpts;
    double coef=q(action.ifirst1)*q(action.ifirst2);
    //double coef=action.ugrid(n-1,0)
    //         *exp((n-1)/action.logrratioinv)/action.rgridinv; 
    for (int i=0; i<n; ++i) {
      double r=exp(i/action.logrratioinv)/action.rgridinv;
      action.ugrid(i,0)-=coef*tau*erf(kappa*r)/r;
      action.ugrid(i,action.norder+1)-=coef*erf(kappa*r)/r;
    }
  }
}

double EwaldAction::calcLongRangeUtau(VArray1& r) const {
  // set up the exponential tables
  for (int ipart=0; ipart<npart; ++ipart) {
    const Complex I(0,1);
    eikx(ipart,0)=eiky(ipart,0)=eikz(ipart,0)=1.0;
    eikx(ipart,1)=exp(I*deltak[0]*r(ipart)[0]);
    eiky(ipart,1)=exp(I*deltak[1]*r(ipart)[1]);
    eikz(ipart,1)=exp(I*deltak[2]*r(ipart)[2]);
    eiky(ipart,-1)=conj(eiky(ipart,1));
    eikz(ipart,-1)=conj(eikz(ipart,1));
    for (int kx=2; kx<=ikmax[0]; ++kx) {
      eikx(ipart,kx)=eikx(ipart,kx-1)*eikx(ipart,1);
    }
    for (int ky=2; ky<=ikmax[1]; ++ky) {
      eiky(ipart, ky)=eiky(ipart, ky-1)*eiky(ipart, 1);
      eiky(ipart,-ky)=eiky(ipart,-ky+1)*eiky(ipart,-1);
    }
    for (int kz=2; kz<=ikmax[2]; ++kz) {
      eikz(ipart, kz)=eikz(ipart, kz-1)*eikz(ipart, 1);
      eikz(ipart,-kz)=eikz(ipart,-kz+1)*eikz(ipart,-1);
    }
  }
  // Sum long range action over all k-vectors.
  double utau=0;
  int ikvec=0;
  for (int kx=0; kx<=ikmax[0]; ++kx) {
    double factor=(kx==0)?1.0:2.0;
    double kx2=kx*kx*deltak[0]*deltak[0];
    for (int ky=-ikmax[1]; ky<=ikmax[1]; ++ky) {
      double ky2=ky*ky*deltak[1]*deltak[1];
      for (int kz=-ikmax[2]; kz<=ikmax[2]; ++kz) {
        double k2=kx2+ky2+kz*kz*deltak[2]*deltak[2];
        if (k2<kcut*kcut && k2!=0) {
          Complex sum=0;
          for (int jpart=0; jpart<npart; ++jpart) {
            sum+=q(jpart)*eikx(jpart,kx)*eiky(jpart,ky)*eikz(jpart,kz);
          }
          utau+=factor*expoverk2(ikvec++)*norm(sum);
        }
      }
    }
  }
  return utau*twoPiOverV;
}

const double EwaldAction::PI=3.14159265358979;
#endif
