// $Id$
/*  Copyright (C) 2007 John B. Shumway, Jr.

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
#include "EIndEstimator.h"
#include "SimulationInfo.h"
#include "action/CoulombAction.h"
#include "SuperCell.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include "stats/MPIManager.h"

EIndEstimator::EIndEstimator(const SimulationInfo& simInfo,
  const CoulombAction* coulAction, const int nfreq, 
  const int nbin, const int ndbin, const int nstride, MPIManager *mpi)
  : BlitzArrayBlkdEst<3>("eind","dynamic-array/induced-e-field",
                         IVecN(2*ndbin-1,nbin,nfreq),true), 
    coulAction(coulAction), npart(simInfo.getNPart()), 
    nslice(simInfo.getNSlice()), nfreq(nfreq), nbin(nbin), ndbin(ndbin),
    nstride(nstride), tau(simInfo.getTau()), tauinv(1./tau), 
    massinv(1./simInfo.getSpecies(0).mass),
    dx(simInfo.getSuperCell()->a[0]/nbin), dxinv(1/dx),q(npart),
    temp(2*ndbin-1,nbin,nslice/nstride), 
    temp2(2*ndbin-1,nbin,nslice/nstride), mpi(mpi) {
  fftw_complex *ptr = (fftw_complex*)temp.data();
  int nsliceEff=nslice/nstride;
  fwd = fftw_plan_many_dft(1,&nsliceEff,nbin,
                           ptr,0,1,nsliceEff,
                           ptr,0,1,nsliceEff,
                           FFTW_FORWARD,FFTW_MEASURE);
  rev = fftw_plan_many_dft(1,&nsliceEff,nbin*(2*ndbin-1),
                           ptr,0,1,nsliceEff,
                           ptr,0,1,nsliceEff,
                           FFTW_BACKWARD,FFTW_MEASURE);
  for (int i=0; i<npart; ++i) q(i)=simInfo.getPartSpecies(i).charge; 
}

EIndEstimator::~EIndEstimator() {
  fftw_destroy_plan(rev);
  fftw_destroy_plan(fwd);
}

void EIndEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0; temp2=0;
}

void EIndEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  //Vec delta=paths.delta(ipart,islice,-1);
  //int ibin = (int)((end[0]-0.5*delta[0]+1.5*nbin*dx)*dxinv)%nbin;
  //temp(0,ibin,islice) += q(ipart) * delta[0] * tauinv;
  //++ninbin(ibin);
  //if (islice%nstride==0) {
    int ibin=((int)(end[0]*dxinv+nbin))%nbin;
    int jbin=((int)(start[0]*dxinv+nbin))%nbin;
    Vec efield=coulAction->getEField(paths,ipart,islice);
  //temp(0,ibin,islice)+=1.0;  
    temp2(0,ibin,islice/nstride)+=efield[0];
    if (ibin!=jbin) {
      //Vec efield=coulAction->getEField(paths,ipart,islice);
      int nstep = ((ibin-jbin+3*nbin/2)%nbin)-nbin/2;
      int idir = (nstep>0)?1:-1;
      if (idir>0) {
        for (int i=1; i<=nstep; i++) {
          temp(0,(jbin+i)%nbin,islice/nstride)+=idir;
        }
        //for (int i=1; i<=nstep; i++) 
        //  temp2(0,(jbin+i)%nbin,islice)+=efield[0]*idir;
      } else {
        for (int i=1; i<=-nstep; i++) {
          temp(0,(ibin+i)%nbin,islice/nstride)+=idir;
        }
        //for (int i=1; i<=-nstep; i++) 
        //  temp2(0,(ibin+i)%nbin,islice)+=efield[0]*idir;
      }
    }
  //}
}

void EIndEstimator::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
  blitz::Range allBin = blitz::Range::all();
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Reduce(&temp(0,0,0),&temp(1,0,0),
                                2*nslice/nstride*nbin,MPI::DOUBLE,MPI::SUM,0);
    temp(0,allBin,allSlice)=temp(1,allBin,allSlice); 
    mpi->getWorkerComm().Reduce(&temp2(0,0,0),&temp2(1,0,0),
                                2*nslice/nstride*nbin,MPI::DOUBLE,MPI::SUM,0);
    temp2(0,allBin,allSlice)=temp2(1,allBin,allSlice); 
  }
#endif
  // Calculate autocorrelation function using FFT's.
  if (workerID==0) {
    //temp/=nslice;
    //temp2/=nslice;
    fftw_execute(fwd);
    fftw_complex *ptr2 = (fftw_complex*)temp2.data();
    fftw_execute_dft(fwd,ptr2,ptr2);
    temp2*=tau;
    for (int ibin=0; ibin<nbin; ++ibin) {
      for (int jdbin=1; jdbin<ndbin; ++jdbin) {
        int jbin=(ibin+jdbin)%nbin;
        temp(jdbin,ibin,allSlice)=conj(temp2(0,ibin,allSlice))
                                      *temp(0,jbin,allSlice);
        jbin=(ibin-jdbin+nbin)%nbin;
        temp(2*ndbin-1-jdbin,ibin,allSlice)=conj(temp2(0,ibin,allSlice))
                                                *temp(0,jbin,allSlice);
      }
    }
    for (int ibin=0; ibin<nbin; ++ibin) {
      temp(0,ibin,allSlice)*=conj(temp2(0,ibin,allSlice));
    }
    //fftw_execute(rev);
    value += imag(temp(allBin,allBin,blitz::Range(0,nfreq-1)))
             /(tau*nslice);
    norm+=1;
  }
}
