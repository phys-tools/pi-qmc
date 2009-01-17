// $Id: DensDensEstimator.cc,v 1.7 2007/08/08 22:14:04 jshumwa Exp $
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
#include "DensDensEstimator.h"
#include "SimulationInfo.h"
#include "Action.h"
#include "DoubleAction.h"
#include "SuperCell.h"
#include <blitz/tinyvec.h>
#include "stats/MPIManager.h"

DensDensEstimator::DensDensEstimator(const SimulationInfo& simInfo,
  const Action* action, const DoubleAction* doubleAction,
  const int nbin, const int ndbin, MPIManager *mpi)
  : BlitzArrayBlkdEst<3>("densityDensity",
                         IVecN(2*ndbin-1,nbin,simInfo.getNSlice()),true), 
    action(action), doubleAction(doubleAction),
    npart(simInfo.getNPart()), nslice(n[2]), nbin(nbin), ndbin(ndbin),
    tauinv(1./simInfo.getTau()), massinv(1./simInfo.getSpecies(0).mass),
    dx(simInfo.getSuperCell()->a[0]/nbin), dxinv(1/dx),
    temp(2*ndbin-1,nbin,nslice), ninbin(nbin), ninbinbuff(nbin), mpi(mpi) {
  fftw_complex *ptr = (fftw_complex*)temp.data();
  fwd = fftw_plan_many_dft(1,&nslice,nbin,
                           ptr,0,1,nslice,
                           ptr,0,1,nslice,
                           FFTW_FORWARD,FFTW_MEASURE);
  rev = fftw_plan_many_dft(1,&nslice,nbin*(2*ndbin-1),
                           ptr,0,1,nslice,
                           ptr,0,1,nslice,
                           FFTW_BACKWARD,FFTW_MEASURE);
}

DensDensEstimator::~DensDensEstimator() {
  fftw_destroy_plan(rev);
  fftw_destroy_plan(fwd);
}

void DensDensEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0; ninbin=0;
}

void DensDensEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  int ibin = (int)((end[0]+0.5*nbin*dx)*dxinv);
  if (ibin>=0 && ibin<nbin) temp(0,ibin,islice)+=1;
}

void DensDensEstimator::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
  blitz::Range allBin = blitz::Range::all();
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Reduce(&temp(0,0,0),&temp(1,0,0),
                                2*nslice*nbin,MPI::DOUBLE,MPI::SUM,0);
    temp(0,allBin,allSlice)=temp(1,allBin,allSlice); 
    //mpi->getWorkerComm().Reduce(ninbin.data(),ninbinbuff.data(),
    //                            nbin,MPI::INTEGER,MPI::SUM,0);
    //ninbin=ninbinbuff;
  }
#endif
  // Calculate autocorrelation function using FFT's.
  if (workerID==0) {
    temp/=nslice;
    fftw_execute(fwd);
    for (int ibin=0; ibin<nbin; ++ibin) {
      for (int jdbin=1; jdbin<ndbin; ++jdbin) {
        int jbin=(ibin+jdbin)%nbin;
        temp(jdbin,ibin,allSlice)=conj(temp(0,ibin,allSlice))
                                      *temp(0,jbin,allSlice);
        jbin=(ibin-jdbin+nbin)%nbin;
        temp(2*ndbin-1-jdbin,ibin,allSlice)=conj(temp(0,ibin,allSlice))
                                                *temp(0,jbin,allSlice);
      }
    }
    for (int ibin=0; ibin<nbin; ++ibin) {
      temp(0,ibin,allSlice)*=conj(temp(0,ibin,allSlice));
    }
    fftw_execute(rev);
    //for (int ibin=0; ibin<nbin; ++ibin) {
    //  temp(0,ibin,0) -= ninbin(ibin)*tauinv*massinv/nslice;
    //}
    value -= real(temp);
    norm+=1;
  }
}
