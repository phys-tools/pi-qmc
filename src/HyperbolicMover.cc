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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "HyperbolicMover.h"
#include "Beads.h"
#include "MultiLevelSampler.h"
#include "RandomNumGenerator.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <fstream>

HyperbolicMover::HyperbolicMover(const SimulationInfo& simInfo,
                                 const int maxlevel)
  : tau(simInfo.getTau()), mass(0.067), alpha(27.211396/1.525),
    normG2(maxlevel+1) {
  normG2=1;
  setFTable(maxlevel);
}

HyperbolicMover::~HyperbolicMover() {
}

HyperbolicMover::Vec HyperbolicMover::sampleDelta(Vec& center,
    const double teff, const int ilevel) {
  blitz::Array<Vec,1> delta(1);
  double r=sqrt(dot(center,center));
  double gr=g(r,teff,ilevel);
  double scaleG1G2=1.0*gr*gr*normG2(ilevel)/g1g2(center,center,teff,ilevel);
  if (scaleG1G2>25.0) {
    /// Go ahead and use squared envelope.
    double p, pEnvelope;
    do {
      RandomNumGenerator::makeGaussRand(delta);
      double r=sqrt(dot(delta(0),delta(0)));
      double scale=fG2(r,ilevel);
      delta(0)*=scale; r*=scale;
      double gr=g(r,teff,ilevel);
      pEnvelope=normG2(ilevel)*gr*gr;
      p=scaleG1G2*g1g2(center,delta(0),teff,ilevel);
//std::cout << ilevel << ": " << p/pEnvelope << ", " << scale << std::endl;
    } while (p/pEnvelope<RandomNumGenerator::getRand());
    if (p/pEnvelope>1.0) std::cout << "!!!!" << p/pEnvelope << std::endl;
  } else {
    /// Use sum envelope instead.
    double p, pEnvelope;
    double dist=sqrt(dot(center,center));
    double scaleGsum=2*g(dist,teff,ilevel)/g1g2(center,0.0,teff,ilevel);
    do {
      RandomNumGenerator::makeGaussRand(delta);
      int dir=(RandomNumGenerator::getRand()>0.5)?+1:-1;
      double r=sqrt(dot(delta(0),delta(0)));
      double scale=fG(r,ilevel);
      delta(0)*=scale; r*=scale;
      pEnvelope=g(r,teff,ilevel);
      delta(0)+=2*dir*center;
      r=sqrt(dot(delta(0),delta(0)));
      pEnvelope+=g(r,teff,ilevel);
      delta(0)-=dir*center;
      p=scaleGsum*g1g2(center,delta(0),teff,ilevel);
//std::cout << ilevel << ": " << p/pEnvelope << ", " << scale << std::endl;
    } while (p/pEnvelope<RandomNumGenerator::getRand());
    if (p/pEnvelope>1.0) {
      std::cout << "!?!!" << p/pEnvelope << center << delta << std::endl;
    }
  }
  return delta(0);
}

double HyperbolicMover::makeMove(MultiLevelSampler& sampler, const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const double teff=tau*nStride;
  const int nSlice=sectionBeads.getNSlice();
  const blitz::Array<int,1>& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  blitz::Array<Vec,1> gaussRand(nMoving);
  double toldOverTnew=1;
  for (int islice=nStride; islice<nSlice-nStride; islice+=2*nStride) {
    RandomNumGenerator::makeGaussRand(gaussRand);
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Calculate the new position.
      Vec midpoint=movingBeads.delta(iMoving,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      Vec delta=sampleDelta(midpoint,teff,level);
      toldOverTnew/=g1g2(midpoint,delta,teff,level);
      midpoint+=movingBeads(iMoving,islice-nStride);
      cell.pbc(midpoint);
      //Vec delta = gaussRand(iMoving);
      //double r=sqrt(dot(delta,delta));
      //delta*=f(r,level);
      (movingBeads(iMoving,islice)=midpoint)+=delta;
      cell.pbc(movingBeads(iMoving,islice));
      // Add transition probability for move.
      //cell.pbc(delta);
      //double r=sqrt(dot(delta,delta));
      //double gr=g(r,teff,level);
      //toldOverTnew/=gr*gr;
      // Calculate and add reverse transition probability.
      midpoint=sectionBeads.delta(i,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      Vec center=midpoint;
      midpoint+=sectionBeads(i,islice-nStride);
      cell.pbc(midpoint);
      delta=sectionBeads(i,islice); 
      delta-=midpoint;
      cell.pbc(delta);
      toldOverTnew*=g1g2(center,delta,teff,level);
      //r=sqrt(dot(delta,delta));
      //gr=g(r,teff,level);
      //toldOverTnew*=gr*gr;
    }
  }
  return toldOverTnew;
}


double HyperbolicMover::fG(const double r, const int ilevel) const {
  int i=(int)(r*drinv); if (i>4998) i=4998;
  double x=r*drinv-i;
  return x*fGtable(i+1,ilevel)+(1-x)*fGtable(i,ilevel); 
}

double HyperbolicMover::fG2(const double r, const int ilevel) const {
  int i=(int)(r*drinv); if (i>4998) i=4998;
  double x=r*drinv-i;
  return x*fG2table(i+1,ilevel)+(1-x)*fG2table(i,ilevel); 
}

double HyperbolicMover::g(const double r, const double t, const int ilevel)
    const {
  double z=t/(2*alpha)*sqrt(1.0+2*mass*alpha*r*r/(t*t));
  double k2=gsl_sf_bessel_Kn(2,z);
  return exp(t/(2*alpha))*t/(32*pi*pi*alpha*z*z)*pow(2*mass/alpha,1.5)*k2;
}

double HyperbolicMover::g0(const double r, const double t) const {
  return pow(2*pi,-1.5)*exp(-0.5*r*r);
}

void HyperbolicMover::setFTable(int maxlevel) {
  std::cout << "Setting FTable" << std::endl;
  fGtable.resize(5000,maxlevel+1);
  fG2table.resize(5000,maxlevel+1);
  std::ofstream file("ftable.dat");
  file << "# Sampling functions" << std::endl;
  file << "# mass=" << mass 
       << ", alpha=" << alpha/27.211 << " eV-1" << std::endl;;
  const double dr=0.001;
  drinv=1.0/dr;
  for (int ilevel=0; ilevel<maxlevel+1; ++ilevel) {
    double sumg=0,sumg0=0;
    double fr=0;
    double teff=tau*pow(2,ilevel);
    /// Set normalization for g and g2.
    normG2(ilevel)=1./g(0,2*teff,ilevel+1);
    sumg=0;
    for (int i=0; i<5000; ++i) {
      double r=i*dr;
      double g0r=g0(r,teff);
      double gfr=g(fr,teff,ilevel);
      gfr*=gfr*normG2(ilevel);
      sumg0+=g0r*2*pi*r*r*dr;
      if (i>0) {
        double oldfr=fr, oldgfr=gfr;
        do {
          fr+=dr;
          gfr=g(fr,teff,ilevel);
          gfr*=gfr*normG2(ilevel);
        } while (2*pi*(oldgfr*oldfr*oldfr+gfr*fr*fr)*(fr-oldfr)<sumg0-sumg);
        double x1=oldfr;
        double x2=fr;
        double y1=2*pi*(oldgfr*oldfr*oldfr+oldgfr*x1*x1)
                  *(x1-oldfr)+sumg-sumg0;
        double y2=2*pi*(oldgfr*oldfr*oldfr+gfr*x2*x2)
                      *(x2-oldfr)+sumg-sumg0;
        double g1=oldgfr;
        double g2=gfr;
        double newfr=x2, newgfr=gfr;
        for (int j=0; j<20; ++j) {
          double x3=0.5*(x1+x2);
          double g3=g(x3,teff,ilevel);
          g3*=g3*normG2(ilevel);
          double y3=2*pi*(oldgfr*oldfr*oldfr+g3*x3*x3)
                        *(x3-oldfr)+sumg-sumg0;
          if (y1*y3<0) {
            x2=x3; y2=y3; g2=g3;
            newfr=x2; newgfr=g2;
          } else { 
            x1=x3; y1=y3; g1=g3;
            newfr=x1; newgfr=g1;
          }
        }
        sumg+=2*pi*(oldgfr*oldfr*oldfr+newgfr*newfr*newfr)*(newfr-oldfr);
        fr=newfr; gfr=newgfr;
        file << r << " " << fr/r << " " << g0r << " " << gfr
                  << " " << sumg0 << " " << sumg << std::endl;
        fG2table(i,ilevel)=fr/r;
      }
      sumg0+=g0r*2*pi*r*r*dr;
    }
    fG2table(0,ilevel)=2*fG2table(2,ilevel)-fG2table(1,ilevel);
    file << std::endl;
    sumg0=sumg=fr=0;
    for (int i=0; i<5000; ++i) {
      double r=i*dr;
      double g0r=g0(r,teff);
      double gfr=g(fr,teff,ilevel);
      sumg0+=g0r*2*pi*r*r*dr;
      if (i>0) {
        double oldfr=fr, oldgfr=gfr;
        do {
          fr+=dr;
          gfr=g(fr,teff,ilevel);
        } while (2*pi*(oldgfr*oldfr*oldfr+gfr*fr*fr)*(fr-oldfr)<sumg0-sumg);
        double x1=oldfr;
        double x2=fr;
        double y1=2*pi*(oldgfr*oldfr*oldfr+oldgfr*x1*x1)
                  *(x1-oldfr)+sumg-sumg0;
        double y2=2*pi*(oldgfr*oldfr*oldfr+gfr*x2*x2)
                      *(x2-oldfr)+sumg-sumg0;
        double g1=oldgfr;
        double g2=gfr;
        double newfr=x2, newgfr=gfr;
        for (int j=0; j<20; ++j) {
          double x3=0.5*(x1+x2);
          double g3=g(x3,teff,ilevel);
          double y3=2*pi*(oldgfr*oldfr*oldfr+g3*x3*x3)
                        *(x3-oldfr)+sumg-sumg0;
          if (y1*y3<0) {
            x2=x3; y2=y3; g2=g3;
            newfr=x2; newgfr=g2;
          } else { 
            x1=x3; y1=y3; g1=g3;
            newfr=x1; newgfr=g1;
          }
        }
        sumg+=2*pi*(oldgfr*oldfr*oldfr+newgfr*newfr*newfr)*(newfr-oldfr);
        fr=newfr; gfr=newgfr;
        fGtable(i,ilevel)=fr/r;
      }
      sumg0+=g0r*2*pi*r*r*dr;
    }
    fGtable(0,ilevel)=2*fGtable(2,ilevel)-fGtable(1,ilevel);
  }
}

const double HyperbolicMover::pi(3.14159265358979);

double HyperbolicMover::g1g2(const Vec center, const Vec delta,
                             const double teff, const int level) {
  double r1=sqrt(dot(center+delta,center+delta));
  double r2=sqrt(dot(center-delta,center-delta));
  double r=sqrt(dot(center,center));
  return g(r1,teff,level)*g(r2,teff,level)/g(2*r,teff*2,level+1);
}
