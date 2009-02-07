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
#include "PairAction.h"
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include <vector>
#include <fstream>
#include <blitz/tinyvec-et.h>

PairAction::PairAction(const Species& s1, const Species& s2,
            const std::string& filename, const SimulationInfo& simInfo, 
            const int norder, const bool isDMD) 
  : tau(simInfo.getTau()), species1(s1), species2(s2),
    ifirst1(s1.ifirst), ifirst2(s2.ifirst),
    npart1(s1.count), npart2(s2.count), norder(norder), isDMD(isDMD) {
std::cout << "constructing PairAction" << std::endl;
std::cout << "species1= " <<  s1 << "species2= " <<  s2 << std::endl;
std::cout << "filename=" << filename << std::endl;
  std::vector<Array> buffer; buffer.reserve(300);
  std::ifstream dmfile((filename+".dmu").c_str()); std::string temp;
  double r;
  double junk;
  Array u(norder+1);
  dmfile >> r; if (isDMD) dmfile >> junk; 
  for (int k=0; k<=norder; ++k) {
    dmfile >> u(k);
    if (isDMD && k>0) u(k)/=r*r;
  }
  getline(dmfile,temp); buffer.push_back(u.copy());
  rgridinv=1./r;
  dmfile >> r; if (isDMD) dmfile >> junk;
  for (int k=0; k<=norder; ++k) {
    dmfile >> u(k); 
    if (isDMD && k>0) u(k)/=r*r;
  }
  getline(dmfile,temp); buffer.push_back(u.copy());
  logrratioinv = 1/log(r*rgridinv);
  while (dmfile) {
    dmfile >> r; if (isDMD) dmfile >> junk;
    for (int k=0; k<=norder; ++k) {
      dmfile >> u(k); 
      if (isDMD && k>0) u(k)/=r*r;
    }
    getline(dmfile,temp); if (dmfile) buffer.push_back(u.copy());
  }
  ngpts=buffer.size();
  ugrid.resize(ngpts, (norder+1)*2);
  for (int i=0; i<ngpts; ++i) {
    for (int k=0; k<=norder; ++k) ugrid(i,k)=buffer[i](k);
  }
  std::ifstream dmefile((filename+".dme").c_str());
  for (int i=0; i<ngpts; ++i) {
    dmefile >> r; if (isDMD) dmefile >> junk;
    for (int k=0; k<=norder; ++k) {
      dmefile >> ugrid(i,k+norder+1);
      if (isDMD && k>0) ugrid(i,k+norder+1)/=r*r;
    }
    getline(dmefile,temp);
  }
  if (isDMD) { // DMD drops potential contribution to tau derivative.
    std::ifstream potfile((filename+".pot").c_str());
    for (int i=0; i<ngpts; ++i) {
      double v;
      potfile >> r >> v; getline(potfile,temp);
      ugrid(i,norder+1)+=v;
    }
  }
  std::ofstream file((filename+".dat").c_str());
  for (int i=0;i<100; ++i) file << i*0.1 << " " << u00(i*0.1) << std::endl;
}

PairAction::PairAction(const Species& s1, const Species& s2,
            const EmpiricalPairAction &action, const SimulationInfo& simInfo, 
            const int norder, const double rmin, const double rmax,
            const int ngpts) 
  : tau(simInfo.getTau()), ngpts(ngpts), rgridinv(1.0/rmin),
    logrratioinv((ngpts-1)/log(rmax/rmin)), ugrid(ngpts,(norder+1)*2),
    species1(s1), species2(s2), ifirst1(s1.ifirst), ifirst2(s2.ifirst),
    npart1(s1.count), npart2(s2.count), norder(norder) {
std::cout << "constructing PairAction" << std::endl;
std::cout << "species1= " <<  s1 << "species2= " <<  s2 << std::endl;
  ugrid=0;
//  std::cout << 1./rgridinv << ", " << logrratioinv << std::endl;
  for (int i=0; i<ngpts; ++i) {
    double r=(1./rgridinv)*exp(i/logrratioinv);
    for (int iorder=0; iorder<norder+1; ++iorder) {
      ugrid(i,iorder)=action.u(r,iorder);
      ugrid(i,iorder+norder+1)=action.utau(r,iorder);
//std::cout << r << " " << ugrid(i,iorder) 
//               << " " << ugrid(i,iorder+norder+1) << std::endl;
    }
  }
}

double PairAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  for (int iMoving=0; iMoving<nMoving; ++iMoving) {
    const int i=index(iMoving);
    int jbegin,jend;
    if (i>=ifirst1 && i<ifirst1+npart1) {
      jbegin=ifirst2; jend=ifirst2+npart2;
    } else
    if (i>=ifirst2 && i<ifirst2+npart2) {
      jbegin=ifirst1; jend=ifirst1+npart1;
    } else 
    continue; //Particle not in this interaction.
    for (int j=jbegin; j<jend; ++j) {
      bool isMoving=false; int jMoving=0;
      for (int k=0;k<nMoving;++k) {
        if (j==index(k)) {isMoving=true; jMoving=k; break;}
      }
      if (isMoving && i<=j) continue; //Don't double count moving interactions.
      Vec prevDelta =sectionBeads(i,0);
      prevDelta-=sectionBeads(j,0); cell.pbc(prevDelta);
      Vec prevMovingDelta=prevDelta;
      double prevR=sqrt(dot(prevDelta,prevDelta));
      double prevMovingR=prevR;
      for (int islice=nStride; islice<nSlice; islice+=nStride) {
        // Add action for moving beads.
        Vec delta=movingBeads(iMoving,islice);
        delta-=(isMoving)?movingBeads(jMoving,islice):sectionBeads(j,islice);
        cell.pbc(delta);
        double r=sqrt(dot(delta,delta));
        double q=0.5*(r+prevMovingR);
       	if (level<3) {
          Vec svec=delta-prevMovingDelta; double s2=dot(svec,svec)/(q*q);
          deltaAction+=uk0(q,s2);
        } else { 
          deltaAction+=u00(q);
        }
        prevMovingDelta=delta; prevMovingR=r;
        // Subtract action for old beads.
        delta=sectionBeads(i,islice);
        delta-=sectionBeads(j,islice);
        cell.pbc(delta);
        r=sqrt(dot(delta,delta));
        q=0.5*(r+prevR);
        if (level<3) {
          Vec svec=delta-prevDelta; double s2=dot(svec,svec)/(q*q);
          deltaAction-=uk0(q,s2);
        } else {
          deltaAction-=u00(q);
        }
        prevDelta=delta; prevR=r;
      }
    }
  }
  return deltaAction*nStride;
}

double PairAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void PairAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  int jbegin,jend;
  if (ipart>=ifirst1 && ipart<ifirst1+npart1) {
    jbegin=ifirst2; jend=ifirst2+npart2;
  } else
  if (ipart>=ifirst2 && ipart<ifirst2+npart2) {
    jbegin=ifirst1; jend=ifirst1+npart1;
  } else 
  return; //Particle not in this interaction.
  for (int j=jbegin; j<jend; ++j) {
    if (ipart==j) continue;
    Vec delta=paths(ipart,islice);
    delta-=paths(j,islice);
    paths.getSuperCell().pbc(delta);
    double r=sqrt(dot(delta,delta));
    Vec prevDelta=paths(ipart,islice,-1);
    prevDelta-=paths(j,islice,-1);
    paths.getSuperCell().pbc(prevDelta);
    double q=0.5*(r+sqrt(dot(prevDelta,prevDelta)));
    Vec svec=delta-prevDelta; double s2=dot(svec,svec)/(q*q);
    double v,vtau,vq,vs2;
    uk0CalcDerivatives(q,s2,v,vtau,vq,vs2);
    u += 0.5*v;
    utau += 0.5*vtau;
    fm -= vq*delta/(2*r) + vs2*(2*svec/(q*q) - s2*delta/(q*r));
    // And force contribution from next slice.
    Vec nextDelta=paths(ipart,islice,+1);
    nextDelta-=paths(j,islice,+1);
    paths.getSuperCell().pbc(nextDelta);
    q=0.5*(r+sqrt(dot(nextDelta,nextDelta)));
    svec=delta-nextDelta; s2=dot(svec,svec)/(q*q);
    uk0CalcDerivatives(q,s2,v,vtau,vq,vs2);
    fp -= vq*delta/(2*r) + vs2*(2*svec/(q*q) - s2*delta/(q*r));
  }
}

double PairAction::u00(double r) const {
  r*=rgridinv; double x=log(r)*logrratioinv;
  int i = (int)x;
  if (r<1) {i=0; x=0;}
  else if (i>ngpts-2) {i=ngpts-2; x=1;}
  else x-=i;
  return (1-x)*ugrid(i,0)+x*ugrid(i+1,0);
}

double PairAction::uk0(double q, double s2) const {
  q*=rgridinv; double x=log(q)*logrratioinv;
  int i = (int)x;
  if (q<1) {i=0; x=0;}
  else if (i>ngpts-2) {i=ngpts-2; x=1;}
  else x-=i;
  double action=0;
  for (int k=norder; k>=0; k--) {
    action*=s2; action+=(1-x)*ugrid(i,k)+x*ugrid(i+1,k);
  }
  return action;
}

void PairAction::uk0CalcDerivatives(double q, double s2, double &u,
            double &utau, double &uq, double &us2) const {
  q*=rgridinv; double x=log(q)*logrratioinv;
  int i = (int)x;
  if (q<1) {i=0; x=0;}
  else if (i>ngpts-2) {i=ngpts-2; x=1;}
  else x-=i;
  int koffset=norder+1;
  u=utau=uq=us2=0;
  for (int k=norder; k>=0; k--) {
    u*=s2; u+=(1-x)*ugrid(i,k)+x*ugrid(i+1,k);
    utau*=s2; utau+=(1-x)*ugrid(i,k+koffset)+x*ugrid(i+1,k+koffset);
    uq*=s2; uq+=ugrid(i+1,k)-ugrid(i,k);
  }
  uq*=logrratioinv*rgridinv/q;
  for (int k=norder; k>0; k--) {
    us2*=s2; us2+=k*((1-x)*ugrid(i,k)+x*ugrid(i+1,k));
  }
}

void PairAction::write(const std::string &fname) const {
  std::string filename=(fname=="")?(species1.name+species2.name):fname;
  std::ofstream ufile((filename+".dmu").c_str()); ufile.precision(6);
  ufile.setf(ufile.scientific,ufile.floatfield); ufile.setf(ufile.showpos);
  std::ofstream efile((filename+".dme").c_str()); efile.precision(6);
  efile.setf(efile.scientific,efile.floatfield); efile.setf(efile.showpos);
  for (int i=0; i<ngpts; ++i) {
    double r=(1./rgridinv)*exp(i/logrratioinv);
    ufile << r; efile << r;
    for (int iorder=0; iorder<norder+1; ++iorder) {
      ufile << " " << ugrid(i,iorder);
      efile << " " << ugrid(i,iorder+norder+1);
    }
    ufile << std::endl; efile << std::endl;
  }
}
