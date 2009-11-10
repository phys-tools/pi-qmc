// $Id$
/*  Copyright (C) 2004-2009 John B. Shumway, Jr.

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
#include "PairIntegrator.h"
#include <vector>
#include <fstream>
#include <blitz/tinyvec-et.h>
#include "DisplaceMoveSampler.h"

PairAction::PairAction(const Species& s1, const Species& s2,
            const std::string& filename, const SimulationInfo& simInfo, 
            const int norder, const bool hasZ, const bool isDMD) 
  : tau(simInfo.getTau()), species1(s1), species2(s2),
    ifirst1(s1.ifirst), ifirst2(s2.ifirst),
    npart1(s1.count), npart2(s2.count), norder(norder), isDMD(isDMD),
    hasZ(hasZ) {
  //std::cout << "constructing PairAction" << std::endl;
  //std::cout << "species1= " <<  s1 << "species2= " <<  s2 << std::endl;
  //std::cout << "filename=" << filename << std::endl;
  std::vector<Array> buffer; buffer.reserve(300);
  std::ifstream dmfile((filename+".dmu").c_str()); std::string temp;
  double r;
  double junk;
  int ndata = (hasZ) ? (norder+1)*(norder+2)/2 : norder+1;
  Array u(ndata);
  dmfile >> r; if (isDMD) dmfile >> junk; 
  for (int k=0; k<ndata; ++k) {
    dmfile >> u(k);
    if (isDMD && k>0) u(k)/=r*r;
  }
  getline(dmfile,temp); buffer.push_back(u.copy());
  rgridinv=1./r;
  dmfile >> r; if (isDMD) dmfile >> junk;
  for (int k=0; k<ndata; ++k) {
    dmfile >> u(k); 
    if (isDMD && k>0) u(k)/=r*r;
  }
  getline(dmfile,temp); buffer.push_back(u.copy());
  logrratioinv = 1/log(r*rgridinv);
  while (dmfile) {
    dmfile >> r; if (isDMD) dmfile >> junk;
    for (int k=0; k<ndata; ++k) {
      dmfile >> u(k); 
      if (isDMD && k>0) u(k)/=r*r;
    }
    getline(dmfile,temp); if (dmfile) buffer.push_back(u.copy());
  }
  ngpts=buffer.size();
  ugrid.resize(ngpts,2,ndata);
  for (int i=0; i<ngpts; ++i) {
    for (int k=0; k<ndata; ++k) ugrid(i,0,k)=buffer[i](k);
  }
  std::ifstream dmefile((filename+".dme").c_str());
  for (int i=0; i<ngpts; ++i) {
    dmefile >> r; if (isDMD) dmefile >> junk;
    for (int k=0; k<ndata; ++k) {
      dmefile >> ugrid(i,1,k);
      if (isDMD && k>0) ugrid(i,1,k)/=r*r;
    }
    getline(dmefile,temp);
  }
  if (isDMD) { // DMD drops potential contribution to tau derivative.
    std::ifstream potfile((filename+".pot").c_str());
    for (int i=0; i<ngpts; ++i) {
      double v;
      potfile >> r >> v; getline(potfile,temp);
      ugrid(i,1,0)+=v;
    }
  }
  std::ofstream file((filename+".dat").c_str());
  for (int i=0;i<100; ++i) file << i*0.1 << " " << u00(i*0.1) << std::endl;
}

PairAction::PairAction(const Species& s1, const Species& s2,
            const EmpiricalPairAction &action, const SimulationInfo& simInfo, 
            const int norder, const double rmin, const double rmax,
            const int ngpts, const bool hasZ) 
  : tau(simInfo.getTau()), ngpts(ngpts), rgridinv(1.0/rmin),
    logrratioinv((ngpts-1)/log(rmax/rmin)), 
    ugrid(ngpts,2,(hasZ?(norder+1)*(norder+2)/2:norder+1)),
    species1(s1), species2(s2), ifirst1(s1.ifirst), ifirst2(s2.ifirst),
    npart1(s1.count), npart2(s2.count), norder(norder), hasZ(hasZ) {
  //std::cout << "constructing PairAction" << std::endl;
  //std::cout << "species1= " <<  s1 << "species2= " <<  s2 << std::endl;
  ugrid=0;
  for (int i=0; i<ngpts; ++i) {
    double r=(1./rgridinv)*exp(i/logrratioinv);
    int idata=0;
    for (int iorder=0; iorder<norder+1; ++iorder) {
      for (int ioff=0; ioff<=(hasZ?iorder:0); ++ioff) {
        ugrid(i,0,idata)=action.u(r,idata);
        ugrid(i,1,idata)=action.utau(r,idata);
        ++idata;
      }
    }
  }
}

PairAction::PairAction(const Species& s1, const Species& s2,
            PairIntegrator &integrator, const SimulationInfo& simInfo, 
            const int norder, const double rmin, const double rmax,
            const int ngpts) 
  : tau(simInfo.getTau()), ngpts(ngpts), rgridinv(1.0/rmin),
    logrratioinv((ngpts-1)/log(rmax/rmin)), 
    ugrid(ngpts,2,(hasZ?(norder+1)*(norder+2)/2:norder+1)),
    species1(s1), species2(s2), ifirst1(s1.ifirst), ifirst2(s2.ifirst),
    npart1(s1.count), npart2(s2.count), norder(norder), hasZ(hasZ) {
std::cout << "constructing PairAction" << std::endl;
std::cout << "species1= " <<  s1 << "species2= " <<  s2 << std::endl;
  ugrid=0;
  int ndata = (norder+1)*(norder+2)/2;
  for (int i=0; i<ngpts; ++i) {
    double r=(1./rgridinv)*exp(i/logrratioinv);
    integrator.integrate(r);
    Array u = integrator.getU();
    for (int idata=0; idata<ndata; ++idata) ugrid(i,0,idata) = u(idata);
    integrator.integrate(r,1.01);
    u = integrator.getU();
    for (int idata=0; idata<ndata; ++idata) ugrid(i,1,idata) = u(idata);
    integrator.integrate(r,0.99);
    u = integrator.getU();
    for (int idata=0; idata<ndata; ++idata) {
      ugrid(i,1,idata) -= u(idata);
      ugrid(i,1,idata) /= 0.02*tau;
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
       	if (level<3 && norder>0) {
          Vec svec=delta-prevMovingDelta; double s2=dot(svec,svec)/(q*q);
          if (hasZ) {
            double z=(r-prevMovingR)/q;
            deltaAction+=uk0(q,s2,z*z);
          } else {
            deltaAction+=uk0(q,s2);
          }
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
        if (level<3 && norder>0) {
          Vec svec=delta-prevDelta; double s2=dot(svec,svec)/(q*q);
          if (hasZ) {
            double z=(r-prevR)/q;
            deltaAction-=uk0(q,s2,z*z);
          } else {
            deltaAction-=uk0(q,s2);
          }
        } else {
          deltaAction-=u00(q);
        }
        prevDelta=delta; prevR=r;
      }
    }
  }
  return deltaAction*nStride;
}

double PairAction::getActionDifference(const Paths &paths, 
    const VArray &displacement, int nmoving, const IArray &movingIndex, 
    int iFirstSlice, int nslice) {
  const SuperCell& cell=paths.getSuperCell();
  double deltaAction=0;
  for (int iMoving=0; iMoving<nmoving; ++iMoving) {
    const int i=movingIndex(iMoving);
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
      for (int k=0;k<nmoving;++k) {
        if (j==movingIndex(k)) {isMoving=true; jMoving=k; break;}
      }
      if (isMoving && i<=j) continue; //Don't double count moving interactions.
      Vec prevDelta =paths(i,iFirstSlice);
      prevDelta-=paths(j,iFirstSlice); cell.pbc(prevDelta);
      Vec prevMovingDelta=prevDelta;
      double prevR=sqrt(dot(prevDelta,prevDelta));
      double prevMovingR=prevR;
      //     for (int islice=iFirstSlice+1; islice<nslice; islice++) {
      for (int islice=iFirstSlice; islice<nslice; islice++) {
        // Add action for moving beads.
        Vec delta=paths(i,islice);
        delta+=displacement(iMoving);
        delta-=paths(j,islice);
        if (isMoving) delta-=displacement(jMoving);
        cell.pbc(delta);
        double r=sqrt(dot(delta,delta));
        double q=0.5*(r+prevMovingR);
      
          Vec svec=delta-prevMovingDelta; double s2=dot(svec,svec)/(q*q);
          deltaAction+=uk0(q,s2);
      
        prevMovingDelta=delta; prevMovingR=r;
        // Subtract action for old beads.
        delta=paths(i,islice);
        delta-=paths(j,islice);
        cell.pbc(delta);
        r=sqrt(dot(delta,delta));
        q=0.5*(r+prevR);
      
          svec=delta-prevDelta;  s2=dot(svec,svec)/(q*q);
          deltaAction-=uk0(q,s2);
      
        prevDelta=delta; prevR=r;
      }
    }
  }
  return deltaAction;
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
    double prevR=sqrt(dot(prevDelta,prevDelta));
    double q=0.5*(r+prevR);
    Vec svec=delta-prevDelta; double s2=dot(svec,svec)/(q*q);
    double v,vtau,vq,vs2,vz2,z;
    if (hasZ) { 
      z=(r-prevR)/q;
      uk0CalcDerivatives(q,s2,z*z,v,vtau,vq,vs2,vz2);
    } else {
      uk0CalcDerivatives(q,s2,v,vtau,vq,vs2);
      vz2=0.; z=0;
    }
    u += 0.5*v;
    utau += 0.5*vtau;
    fm -= vq*delta/(2*r) + vs2*(2*svec/(q*q) - s2*delta/(q*r))
         +vz2*z*delta*(2-z)/(q*r);
    // And force contribution from next slice.
    Vec nextDelta=paths(ipart,islice,+1);
    nextDelta-=paths(j,islice,+1);
    paths.getSuperCell().pbc(nextDelta);
    double nextR=sqrt(dot(nextDelta,nextDelta));
    q=0.5*(r+nextR);
    svec=delta-nextDelta; s2=dot(svec,svec)/(q*q);
    if (hasZ) { 
      z=(r-nextR)/q;
      uk0CalcDerivatives(q,s2,z*z,v,vtau,vq,vs2,vz2);
    } else {
      uk0CalcDerivatives(q,s2,v,vtau,vq,vs2);
      vz2=0.; z=0.;
    }
    fp -= vq*delta/(2*r) + vs2*(2*svec/(q*q) - s2*delta/(q*r))
         +vz2*z*delta*(2-z)/(q*r);
  }
 
}

double PairAction::u00(double r) const {
  r*=rgridinv; double x=log(r)*logrratioinv;
  int i = (int)x;
  if (r<1) {i=0; x=0;}
  else if (i>ngpts-2) {i=ngpts-2; x=1;}
  else x-=i;
  return (1-x)*ugrid(i,0,0)+x*ugrid(i+1,0,0);
}

double PairAction::uk0(double q, double s2) const {
  q*=rgridinv; double x=log(q)*logrratioinv;
  int i = (int)x;
  if (q<1) {i=0; x=0;}
  else if (i>ngpts-2) {i=ngpts-2; x=1;}
  else x-=i;
  double action=0;
  for (int k=norder; k>=0; k--) {
    action*=s2; action+=(1-x)*ugrid(i,0,k)+x*ugrid(i+1,0,k);
  }
  return action;
}

double PairAction::uk0(double q, double s2, double z2) const {
  q*=rgridinv; double x=log(q)*logrratioinv;
  int i = (int)x;
  if (q<1) {i=0; x=0;}
  else if (i>ngpts-2) {i=ngpts-2; x=1;}
  else x-=i;
  double action=0;
  for (int l=norder; l>=0; l--) {
    int index=norder*(norder+1)/2+l;
    double temp=0;
    for (int k=norder-l; k>=0; k--) {
      temp*=s2; temp+=(1-x)*ugrid(i,0,index)+x*ugrid(i+1,0,index);
      index -= k+l;
    }
    action*=z2; action+=temp;
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
  u=utau=uq=us2=0;
  for (int k=norder; k>=0; k--) {
    u*=s2; u+=(1-x)*ugrid(i,0,k)+x*ugrid(i+1,0,k);
    utau*=s2; utau+=(1-x)*ugrid(i,1,k)+x*ugrid(i+1,1,k);
    uq*=s2; uq+=ugrid(i+1,0,k)-ugrid(i,0,k);
  }
  uq*=logrratioinv*rgridinv/q;
  for (int k=norder; k>0; k--) {
    us2*=s2; us2+=k*((1-x)*ugrid(i,0,k)+x*ugrid(i+1,0,k));
  }
}

void PairAction::uk0CalcDerivatives(double q, double s2, double z2, double &u,
            double &utau, double &uq, double &us2, double& uz2) const {
  q*=rgridinv; double x=log(q)*logrratioinv;
  int i = (int)x;
  if (q<1) {i=0; x=0;}
  else if (i>ngpts-2) {i=ngpts-2; x=1;}
  else x-=i;
  u=utau=uq=us2=uz2=0.;
  for (int l=norder; l>=0; l--) {
    int index=norder*(norder+1)/2+l;
    double a=0, atau=0, aq=0;
    for (int k=norder-l; k>=0; k--) {
      a*=s2; a+=(1-x)*ugrid(i,0,index)+x*ugrid(i+1,0,index);
      atau*=s2; atau+=(1-x)*ugrid(i,1,index)+x*ugrid(i+1,1,index);
      aq*=s2; aq += ugrid(i+1,0,index)-ugrid(i,0,index);
      index -= k+l;
    }
    u*=z2; u+=a;
    utau*=z2; utau+=atau;
    uq*=z2; uq+=aq;
  }
  uq*=logrratioinv*rgridinv/q;
  for (int l=norder; l>=0; l--) {
    int index=norder*(norder+1)/2+l;
    double as2=0;
    for (int k=norder-l; k>0; k--) {
      as2*=s2; as2+=k*(1-x)*ugrid(i,0,index)+x*ugrid(i+1,0,index);
      index -= k+l;
    }
    us2*=z2; us2+=as2;
  }
  for (int l=norder; l>0; l--) {
    int index=norder*(norder+1)/2+l;
    double az2=0;
    for (int k=norder-l; k>=0; k--) {
      az2*=s2; az2+=(1-x)*ugrid(i,0,index)+x*ugrid(i+1,0,index);
    }
    uz2*=z2; uz2+=l*az2;
  }
}

void PairAction::write(const std::string &fname) const {
  std::string filename=(fname=="")?(species1.name+species2.name):fname;
  std::ofstream ufile((filename+".dmu").c_str()); ufile.precision(6);
  ufile.setf(ufile.scientific,ufile.floatfield); ufile.setf(ufile.showpos);
  std::ofstream efile((filename+".dme").c_str()); efile.precision(6);
  efile.setf(efile.scientific,efile.floatfield); efile.setf(efile.showpos);
  int ndata = (hasZ) ? (norder+1)*(norder+2)/2 : norder+1;
  for (int i=0; i<ngpts; ++i) {
    double r=(1./rgridinv)*exp(i/logrratioinv);
    ufile << r; efile << r;
    for (int idata=0; idata<ndata; ++idata) {
      ufile << " " << ugrid(i,0,idata);
      efile << " " << ugrid(i,1,idata);
    }
    ufile << std::endl; efile << std::endl;
  }
}
