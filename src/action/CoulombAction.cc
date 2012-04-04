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
#include "CoulombAction.h"
#include "sampler/SectionSamplerInterface.h"
#include "sampler/DisplaceMoveSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PairAction.h"
#include "ImagePairAction.h"
#include "EwaldImagePairAction.h"
#include "util/TradEwaldSum.h"
#include "util/OptEwaldSum.h"

CoulombAction::CoulombAction(const double epsilon, 
  const SimulationInfo& simInfo, const int norder, double rmin, double rmax,
  int ngpts, const bool dumpFiles, bool useEwald, int ewaldNDim, 
  double ewaldRcut, double ewaldKcut, double screenDist, const double kappa, 
  const int nImages, const std::string ewaldType, int exLevel)
  : epsilon(epsilon), tau(simInfo.getTau()), npart(simInfo.getNPart()),
    pairActionArray(0), screenDist(screenDist), ewaldSum(0) {
  typedef blitz::TinyVector<int, NDIM> IVec;
  const int nspecies=simInfo.getNSpecies();
  if (useEwald && (NDIM==3||NDIM==2) && ewaldNDim==NDIM) {
    SuperCell &cell(*simInfo.getSuperCell());
    if (ewaldRcut==0.) ewaldRcut = cell.a[0]/2.;
    std::cout << "EwaldRcut = " << ewaldRcut << std::endl;
  
    if (ewaldType=="tradEwald") ewaldSum = new TradEwaldSum(cell,npart,ewaldRcut,ewaldKcut, kappa);
    if (ewaldType=="optEwald")  ewaldSum = new OptEwaldSum(cell,npart,ewaldRcut,ewaldKcut,4*ewaldKcut,8);
    rewald.resize(npart);
    EwaldSum::Array &q=ewaldSum->getQArray();  
    for (int i=0; i<npart; ++i) q(i)=simInfo.getPartSpecies(i).charge;
    ewaldSum->evalSelfEnergy();
  }
  for (int i=0; i<nspecies; ++i) {
    for (int j=i; j<nspecies; ++j) {
      const Species &s1=simInfo.getSpecies(i);
      const Species &s2=simInfo.getSpecies(j);
      q1q2=s1.charge*s2.charge/epsilon;
      if (q1q2!=0 && ((s1.name!=s2.name) || s1.count>1) ) {
        mu=1./(1./s1.mass+1./s2.mass);
        displace2 = s1.displace-s2.displace; displace2 *= displace2;
        /// Hack to handle ion-ion interaction properly.
        if (rmin==0) rmin=0.005/(2*((mu>500)?1:mu)*fabs(q1q2)); 
        IVec nimage=0; bool needImages=false;
        if (useEwald && ewaldNDim<NDIM) {
          SuperCell::Vec a=simInfo.getSuperCell()->a;
          for (int idim=0; idim<ewaldNDim; ++idim) {
            nimage[idim]=(int)(rmax/a[idim]);
            if (nimage[idim]>1) nimage[idim]-=2;
            if (nimage[idim]>0) needImages=true;
          }
        } else {
          if (rmax==0) {
            SuperCell::Vec a=simInfo.getSuperCell()->a;
            rmax=sqrt(dot(a,a));
          }
        }
        if (ngpts==0) ngpts=500;
        if (needImages) {
          pairActionArray.push_back(
            // Hack to handle ion-ion interaction properly.
            new ImagePairAction(s1,s2, *this, simInfo, (mu>500)?0:norder,
                  nimage, rmin, rmax, ngpts,(i==j)?exLevel:-1));
        } else {
	  if (nImages > 1 && ewaldType=="tradEwald"){
	    pairActionArray.push_back(
              new EwaldImagePairAction(s1,s2, *this, simInfo, (mu>500)?0:norder,
                    rmin, rmax, ngpts, nImages, (i==j)?exLevel:-1));
	  } else {
	    pairActionArray.push_back(
              // Hack to handle ion-ion interaction properly.
            new PairAction(s1,s2, *this, simInfo, (mu>500)?0:norder,
                           rmin, rmax, ngpts, false, (i==j)?exLevel:-1));
	  }
        }
        if (dumpFiles) (*(pairActionArray.end()-1))->write("",0);
      }
    }
  }
}

CoulombAction::~CoulombAction() {
  for (unsigned int i=0; i<pairActionArray.size(); ++i) {
     delete pairActionArray[i];
  }
  delete ewaldSum;
}

double CoulombAction::getActionDifference(const SectionSamplerInterface& sampler,
                                         const int level) {
  double u=0;
  for (unsigned int i=0; i<pairActionArray.size(); ++i) {
    u += pairActionArray[i]->getActionDifference(sampler,level);
  }
  // Compute long range Ewald action at lowest level.
  if (ewaldSum && level==0) {
    const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
    const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
    const int nSlice=sectionBeads.getNSlice();
    const IArray& index=sampler.getMovingIndex(); 
    const int nMoving=index.size();
    for (int islice=1; islice<nSlice-1; ++islice) {
      for (int i=0; i<npart; ++i) rewald(i)=sectionBeads(i,islice); 
      u -= ewaldSum->evalLongRange(rewald)*tau/epsilon; 
      for (int i=0; i<nMoving; ++i) rewald(index(i))=movingBeads(i,islice);
      u += ewaldSum->evalLongRange(rewald)*tau/epsilon; 
    }
  } 
  return u;
}

double CoulombAction::getActionDifference(const Paths &paths, 
    const VArray &displacement, int nmoving, const IArray &movingIndex, 
    int iFirstSlice, int iLastSlice) {

  SuperCell cell=paths.getSuperCell();
  double u=0;
  for (unsigned int i=0; i<pairActionArray.size(); ++i) {
    u += pairActionArray[i]->getActionDifference(paths, displacement, nmoving,
                               movingIndex, iFirstSlice, iLastSlice); 
  }
  // Compute long range Ewald action at lowest level.
  if (ewaldSum) {
     for (int islice=iFirstSlice; islice<=iLastSlice; ++islice) {
      for (int i=0; i<npart; ++i)  rewald(i)=paths(i,islice);
      u -= ewaldSum->evalLongRange(rewald)*tau/epsilon;
      for (int i=0; i<nmoving; ++i) {
        rewald(movingIndex(i))+=displacement(i);
        rewald(movingIndex(i))=cell.pbc(rewald(movingIndex(i)));
      }
      u += ewaldSum->evalLongRange(rewald)*tau/epsilon;
    }
  }
  return u;
}



double CoulombAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void CoulombAction::getBeadAction(const Paths& paths, int ipart, int islice,
         double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  for (unsigned int i=0; i<pairActionArray.size(); ++i) {
    double ui=0, utaui=0, ulambdai=0; 
    Vec fmi=0.0, fpi=0.0; 
    pairActionArray[i]->
      getBeadAction(paths,ipart,islice,ui,utaui,ulambdai,fmi,fpi);
    u+=ui; utau+=utaui; ulambda+=ulambdai; fm+=fmi; fp+=fpi;
  }
  // Compute long range action.
  if (ewaldSum && ipart==0) {
    paths.getSlice(islice,rewald);
    double longRange = ewaldSum->evalLongRange(rewald)/epsilon;
    u += longRange*tau;
    utau += longRange;
  }
  //std :: cout << "CA :: "<<ipart<<" "<<islice<<"  "<<utau<<"  "<<u<<std ::endl;
}

double CoulombAction::getAction(const Paths& paths, int islice) const {
  double u=0;
  for (int ipart=0; ipart<npart; ++ipart) {
    for (unsigned int i=0; i<pairActionArray.size(); ++i) {
      double ui=0, utaui=0, ulambdai=0; 
      Vec fmi=0.0, fpi=0.0; 
      pairActionArray[i]->
        getBeadAction(paths,ipart,islice,ui,utaui,ulambdai,fmi,fpi);
      u+=ui;
    }
  }
  // Compute long range action.
  if (ewaldSum) {
    paths.getSlice(islice,rewald);
    double longRange = ewaldSum->evalLongRange(rewald)/epsilon;
    u += longRange*tau;
  }
  return u;
}

double CoulombAction::getEField(const Paths& paths, int ipart, 
                                int islice) const {
  double u=0,utau=0,ulambda=0;
  Vec fm=0.,fp=0.;
  for (unsigned int i=0; i<pairActionArray.size(); ++i) {
    double ui=0, utaui=0, ulambdai=0; 
    Vec fmi=0.0, fpi=0.0; 
    pairActionArray[i]->
      getBeadAction(paths,ipart,islice,ui,utaui,ulambdai,fmi,fpi);
    u+=ui; utau+=utaui; ulambda+=ulambdai; fm+=fmi; fp+=fpi;
  }
  return fp[0]/tau;
}

double CoulombAction::u(double r, int order) const {
  r = sqrt(r*r+displace2);
  double taueff = 2.0*mu*q1q2*q1q2*tau;
  double u=0;
  double reff = 2.0*mu*q1q2*r;
  double stau = q1q2*sqrt(2.0*mu*tau);
  //First compute the 1D actions.
  double u0 = stau*(1.772453851+stau*(-0.074137740+stau*(0.005834805
             +stau*(-0.000382686+stau*(0.000008738+stau*0.000002138)))));
  double u1_0=0., u1_1=0., u1_2=0., u1_3=0., u1_4=0.;
  { // u0 for 1D
    double a= 0.25300593*pow(stau,-1)+0.01432126;
    double b= 0.07936898*pow(stau,-2)-0.01634421/stau;
    double c= 0.07263383*pow(stau,-3);
    double d= 0.00940013*pow(stau,-4);
    double e= 0.03160181*pow(stau,-5);
    double f= 0.18976335*pow(stau,-1);
    double g= 0.00343053*pow(stau,-2);
    u1_0 = (u0 + reff*( (u0*a-1.) + reff*(f+reff*(g+reff*e*taueff)))) 
         / (1.+reff*(a+reff*(b+reff*(c+reff*(d+reff*e)))));
  }
  { // u1 for 1D
    double a= 0.07776906*pow(stau,-1);
    double b=-0.04695226*pow(stau,-2);
    double c= 0.01291994*pow(stau,-3);
    double d= 0.60942687*pow(stau,-1);
    double e=-1.09588064*pow(stau,-2);
    double f= 1.15591565*pow(stau,-3);
    double g=-0.57485278*pow(stau,-4);
    double h= 0.15212765*pow(stau,-5);
    u1_1 = pow(reff,2)*(a+reff*(b+reff*c))
           /(1+reff*(d+reff*(e+reff*(f+reff*(g+reff*h)))));
  }
  { // u2 for 1D
    double a= 9.52130652e-04*pow(stau,-3);
    double b=-6.39946071e-04*pow(stau,-4);
    double c= 1.14359558e-04*pow(stau,-5);
    double d=-7.59575215e-01*pow(stau,-3);
    double e= 5.53798980e-01*pow(stau,-4);
    double f=-7.15626226e-02*pow(stau,-5);
    double g=-3.49081618e-02*pow(stau,-6);
    double h= 8.40084881e-03*pow(stau,-7);
    u1_2 = pow(reff,4)*(a+reff*(b+reff*c))
          /(1+pow(reff,3)*(d+reff*(e+reff*(f+reff*(g+reff*h)))));
  }
  { // u3 for 1D
    double a=-2.06984256e-1*pow(stau,-5); 
    double b=-9.51423947e-3*pow(stau,-6); 
    double c= 8.97532561e-2*pow(stau,-7); 
    double d= 6.10876797e+4*pow(stau,-3); 
    double e=-1.19016292e+5*pow(stau,-4); 
    double f= 9.29836599e+4*pow(stau,-5); 
    double g=-3.49919130e+4*pow(stau,-6); 
    double h= 6.42490539e+3*pow(stau,-7); 
    double i=-6.02729064e+2*pow(stau,-8); 
    double j= 6.01285550e+1*pow(stau,-9); 
    u1_3 = pow(reff,6)*(a+reff*(b+reff*c))
       /(1+pow(reff,3)*(d+reff*(e+reff*(f+reff*(g+reff*(h+reff*(i+reff*j)))))));
  }
  { // u4 for 1D
    double a=-2.84102296e-6*pow(stau,-5);
    double b= 6.26961672e-7*pow(stau,-7); 
    double c=-2.20166473e-1*pow(stau,-3); 
    double d= 1.07981903e-1*pow(stau,-5); 
    double e=-1.94671111e-2*pow(stau,-7); 
    double f= 1.56217930e-3*pow(stau,-9); 
    u1_4 = pow(reff,6)*(a+b*reff*reff)
          /(1+reff*reff*reff*(c+reff*reff*(d+reff*reff*(e+reff*reff*f))));
  }
#if NDIM==1
  if (order==0) {
    u = u1_0;
  } else if (order==1) {
    u = u1_1;
  } else if (order==2) {
    u = u1_2;
  } else if (order==3) {
    u = u1_3;
  } else if (order==4) {
    u = u1_4;
  }
#endif
#if NDIM==3
  double a=2.*taueff/(reff*reff);
  double b=exp(-2./a);
  double c=1./(1+2*a*(1-b)*u1_1);
  // Equations from Mathematica's power series expansion.
  if (order==0) {
    u = u1_0 + log(c);
    if (ewaldSum) {
      u -= tau*q1q2*ewaldSum->evalFR(r);
    }
  } else if (order==1) {
    u = u1_1 + c*(b*u1_1 - 4*a*(1-b)*u1_2);
  } else if (order==2) {
    u = u1_2 +0.25*c/a*(c* ( b*u1_1*(1+2*a*u1_1) + 8*a*b*u1_2
                          +32*a*a*a*(1-b)*(1-b)*u1_2*u1_2) - 24*a*a*(1-b)*u1_3);
  } else if (order==3) {
    u = u1_3 + c*c*c/(24.*a*a)*(
         128*a*a*a*a*a*b*b*b*(4*u1_2*u1_2*u1_2-9*u1_1*u1_2*u1_3
                            +6*u1_1*u1_1*u1_4)
         +64*a*a*a*(a*u1_2*(-8*a*u1_2*u1_2+9*(1+2*a*u1_1)*u1_3)
                  -3*(1+2*a*u1_1)*(1+2*a*u1_1)*u1_4)
         +2*a*b*b*(2*a*u1_1*u1_1*u1_1
                 +96*a*a*u1_2*(u1_2-8*a*a*u1_2*u1_2+3*a*u1_3)
                 +12*a*u1_1*(u1_2-6*a*u1_3+144*a*a*a*u1_2*u1_3-32*a*a*u1_4)
                 +u1_1*u1_1*(1-1152*a*a*a*a*u1_4))
         +b*(4*a*a*u1_1*u1_1*u1_1+12*a*(u1_2-16*a*a*u1_2*u1_2
                                      +128*a*a*a*a*u1_2*u1_2*u1_2+6*a*u1_3
                                      -96*a*a*a*u1_2*u1_3+16*a*a*u1_4)
            +4*u1_1*u1_1*(a+576*a*a*a*a*a*u1_4)
            +u1_1*(1+24*a*a*(u1_2+6*a*u1_3-144*a*a*a*u1_2*u1_3+64*a*a*u1_4))));

  } else if (order==4) {
    u = u1_4; //Not exact
  }
#endif
#if NDIM==2
  //double ureff = 2.0*mu*fabs(q1q2)*r;
  //double sign= (q1q2>0)?1:-1;
  double a=2.*taueff/(reff*reff);
  double b=exp(-2./a);
  double c=1./(1+2*a*(1-b)*u1_1);
  double u3_0 = u1_0 + log(c);
  double u3_1 = u1_1 + c*(b*u1_1 - 4*a*(1-b)*u1_2);
  if (order==0) {
    //double reff = 2.0*mu*q1q2*r;
    //double stau = q1q2*sqrt(2.0*mu*tau);
    //double u0 = stau*(2.784163998+stau*(-0.331414576+stau*(0.069113708
    //           +stau*(-0.007701621+stau*(-0.000062512+stau*0.000335724)))));
    //double a = 0.81644994/stau + 0.20750388;
    //double b = 0.73884476/taueff + 0.22874030/stau;
    //double c = 0.48135991*pow(stau,-3) + 0.15828689*pow(stau,-2);
    //double d = 0.15769661*pow(stau,-4) + 0.01512610*pow(stau,-3);
    //double e = 0.25472475*pow(stau,-5) + 0.16019595*pow(stau,-4);
    //double f = 0.51100692/stau + 0.15237602;
    //double g = 0.15705841/taueff + 0.01503638/stau;
    //u = (u0 + reff*( (-2.+u0*a) + reff*(f + reff*(g + e*taueff*reff)))) /
    //    (1. + reff*(a + reff*(b + reff*(c + reff*(d + e*reff)))));
    u = u3_0 + 0.5*log(1.+4*stau*stau*u3_1/(reff*reff));
    if (ewaldSum) {
      u -= tau*q1q2*ewaldSum->evalFR(r);
    }
    if (screenDist!=0) {
      //ureff=2.0*mu*fabs(q1q2)*sqrt(r*r+screenDist*screenDist);
      //u -= (u0 + ureff*( (-2.+u0*a) + ureff*(f + ureff*(g + e*taueff*ureff)))) /
      //     (1. + ureff*(a + ureff*(b + ureff*(c + ureff*(d + e*ureff)))));
      u -= q1q2*tau/sqrt(r*r+screenDist*screenDist);
    }
  } else if (order==1) {
    u = u3_1;
  }
#endif
 
  return u;
}

double CoulombAction::utau(double r, int order) const {
  double tausave=tau;
  const double eps=1e-7;
  tau = tausave*(1.0+eps);
  double up=u(r,order);
  tau = tausave*(1.0-eps);
  double um=u(r,order);
  tau = tausave;
  return (up-um)/(2.0*eps*tau);
}
