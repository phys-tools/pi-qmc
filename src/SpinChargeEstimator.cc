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
#include "SpinChargeEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include <blitz/tinyvec.h>
#include "stats/MPIManager.h"

SpinChargeEstimator::SpinChargeEstimator(const SimulationInfo& simInfo,
  const Species &sup, const Species &sdn, const int nfreq,
  const int nbin, const int ndbin, const int nstride, MPIManager *mpi)
  : BlitzArrayBlkdEst<5>("spincharge",
                         IVecN(2*ndbin-1,nbin,2,2,nfreq),true), 
    npart(simInfo.getNPart()), nslice(simInfo.getNSlice()),
    nfreq(nfreq), nbin(nbin), ndbin(ndbin), nstride(nstride),
    nup(sup.count), ndn(sdn.count), ifirstup(sup.ifirst), ifirstdn(sdn.ifirst),
    beta(1./simInfo.getTemperature()), tau(simInfo.getTau()),
    tauinv(1./tau), massinv(1./simInfo.getSpecies(0).mass),
    dx(simInfo.getSuperCell()->a[0]/nbin), dxinv(1/dx),q(npart),
    temp(2*ndbin-1,nbin,2,2,nslice/nstride), mpi(mpi) {
  fftw_complex *ptr = (fftw_complex*)temp.data();
  int nsliceEff=nslice/nstride;
  fwd = fftw_plan_many_dft(1,&nsliceEff,4*nbin,
                           ptr,0,1,nsliceEff,
                           ptr,0,1,nsliceEff,
                           FFTW_FORWARD,FFTW_MEASURE);
  rev = fftw_plan_many_dft(1,&nsliceEff,4*nbin*(2*ndbin-1),
                           ptr,0,1,nsliceEff,
                           ptr,0,1,nsliceEff,
                           FFTW_BACKWARD,FFTW_MEASURE);
  for (int i=0; i<npart; ++i) q(i)=simInfo.getPartSpecies(i).charge; 
}

SpinChargeEstimator::~SpinChargeEstimator() {
  fftw_destroy_plan(rev);
  fftw_destroy_plan(fwd);
}

void SpinChargeEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0;
}

void SpinChargeEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  //Vec delta=paths.delta(ipart,islice,-1);
  //int ibin = (int)((end[0]-0.5*delta[0]+1.5*nbin*dx)*dxinv)%nbin;
  //temp(0,ibin,islice) += q(ipart) * delta[0] * tauinv;
  //++ninbin(ibin);
  int ispin=0;
  if (ipart>=ifirstup && ipart<ifirstup+nup) {
    ispin=0;
  } else if (ipart>=ifirstdn && ipart<ifirstdn+ndn) {
    ispin=1;
  } else {
    return;
  }
  int ibin=((int)(end[0]*dxinv+nbin))%nbin;
  int jbin=((int)(start[0]*dxinv+nbin))%nbin;
  if (ibin!=jbin) {
    int nstep = ((ibin-jbin+3*nbin/2)%nbin)-nbin/2;
    int idir = (nstep>0)?1:-1;
    if (idir>0) {
      for (int i=1; i<=nstep; i++) temp(0,(jbin+i)%nbin,0,ispin,islice/nstride)+=idir;
    } else {
      for (int i=1; i<=-nstep; i++) temp(0,(ibin+i)%nbin,0,ispin,islice/nstride)+=idir;
    }
  }
}

void SpinChargeEstimator::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
  blitz::Range allBin = blitz::Range::all();
  blitz::Range allSpin = blitz::Range::all();
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Reduce(&temp(0,0,0,0,0),&temp(1,0,0,0,0),
                                8*nslice/nstride*nbin,MPI::DOUBLE,MPI::SUM,0);
    temp(0,allBin,allSpin,allSpin,allSlice)
     =temp(1,allBin,allSpin,allSpin,allSlice); 
  }
#endif
  // Calculate autocorrelation function using FFT's.
  if (workerID==0) {
    //temp/=nslice;
    fftw_execute(fwd);
    for (int ibin=0; ibin<nbin; ++ibin) {
      for (int jdbin=1; jdbin<ndbin; ++jdbin) {
        // up-up
        int jbin=(ibin+jdbin)%nbin;
        temp(jdbin,ibin,0,0,allSlice)
         =conj(temp(0,ibin,0,0,allSlice))
              *temp(0,jbin,0,0,allSlice);
        jbin=(ibin-jdbin+nbin)%nbin;
        temp(2*ndbin-1-jdbin,ibin,0,0,allSlice)
          =conj(temp(0,ibin,0,0,allSlice))
               *temp(0,jbin,0,0,allSlice);
        // up-down
        jbin=(ibin+jdbin)%nbin;
        temp(jdbin,ibin,0,1,allSlice)
         =conj(temp(0,ibin,0,0,allSlice))
          *temp(0,jbin,0,1,allSlice);
        jbin=(ibin-jdbin+nbin)%nbin;
        temp(2*ndbin-1-jdbin,ibin,0,1,allSlice)
          =conj(temp(0,ibin,0,0,allSlice))
               *temp(0,jbin,0,1,allSlice);
        // down-up
        jbin=(ibin+jdbin)%nbin;
        temp(jdbin,ibin,1,0,allSlice)
         =conj(temp(0,ibin,0,1,allSlice))
              *temp(0,jbin,0,0,allSlice);
        jbin=(ibin-jdbin+nbin)%nbin;
        temp(2*ndbin-1-jdbin,ibin,1,0,allSlice)
          =conj(temp(0,ibin,0,1,allSlice))
               *temp(0,jbin,0,0,allSlice);
        // down-down
        jbin=(ibin+jdbin)%nbin;
        temp(jdbin,ibin,1,1,allSlice)
         =conj(temp(0,ibin,0,1,allSlice))
              *temp(0,jbin,0,1,allSlice);
        jbin=(ibin-jdbin+nbin)%nbin;
        temp(2*ndbin-1-jdbin,ibin,1,1,allSlice)
          =conj(temp(0,ibin,0,1,allSlice))
               *temp(0,jbin,0,1,allSlice);
      }
    }
    for (int ibin=0; ibin<nbin; ++ibin) {
      // down-up
      temp(0,ibin,1,0,allSlice)
        =conj(temp(0,ibin,0,1,allSlice))
             *temp(0,ibin,0,0,allSlice);
      // down-down
      temp(0,ibin,1,1,allSlice)
        =conj(temp(0,ibin,0,1,allSlice))
             *temp(0,ibin,0,1,allSlice);
      // up-up
      temp(0,ibin,0,0,allSlice)*=conj(temp(0,ibin,0,0,allSlice));
      // up-down
      temp(0,ibin,0,1,allSlice)=conj(temp(0,ibin,1,0,allSlice));
    }
    //fftw_execute(rev);
    value += real(temp(allBin,allBin,allSpin,allSpin,blitz::Range(0,nfreq-1)))
             /beta;
    norm+=1;
  }
}
