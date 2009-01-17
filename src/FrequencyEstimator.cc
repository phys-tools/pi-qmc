// $Id: FrequencyEstimator.cc,v 1.7 2007/10/03 12:53:56 jshumwa Exp $
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
#include "FrequencyEstimator.h"
#include "SimulationInfo.h"
#include "Action.h"
#include "DoubleAction.h"
#include "SuperCell.h"
#include "Species.h"
#include "Paths.h"
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"

FrequencyEstimator::FrequencyEstimator(const SimulationInfo& simInfo,
  const Species& species1, const Species& species2, MPIManager *mpi)
  : BlitzArrayBlkdEst<1>("frequency", IVecN(simInfo.getNSlice() ),true),
    npart(simInfo.getNPart()), nslice(simInfo.getNSlice()), 
    tau(simInfo.getTau()),
    names(simInfo.getNPart()), spec1(species1.name), spec2(species2.name),
    temp(nslice), mpi(mpi){
  for (int i=0; i<names.size(); ++i)
    names(i)=simInfo.getPartSpecies(i).name;
  std::cout << "Frequency Estimator " << spec1 << " " << spec2 << std::endl;
  fftw_complex *ptr = (fftw_complex*)temp.data();
  fwd = fftw_plan_dft_1d(nslice, ptr, ptr, FFTW_FORWARD, FFTW_MEASURE);
  rev = fftw_plan_dft_1d(nslice, ptr, ptr, FFTW_BACKWARD, FFTW_MEASURE);
}

FrequencyEstimator::~FrequencyEstimator() {
  fftw_destroy_plan(fwd);
  fftw_destroy_plan(rev);
}

void FrequencyEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0 ;
}

void FrequencyEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  if(names(ipart)==spec1) {
     for(int jpart=0; jpart<names.size();++jpart){
        if(names(jpart)==spec2 && jpart!=ipart) {
           Vec delta=paths(ipart,islice)-paths(jpart,islice);
  	   double delta2=dot(delta, delta);
           temp(islice)= sqrt(delta2);
        }
     }
  }
}

void FrequencyEstimator::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
 // First move all data to 1st worker. 
 // int workerID=(mpi)?mpi->getWorkerID():0;
 #ifdef ENABLE_MPI
  if (mpi) {
//    mpi->getWorkerComm().Reduce(&temp(0),&temp(0),
//                                nslice*nbin,MPI::DOUBLE,MPI::SUM,0);
//    temp(0,allBin,allSlice)=temp(1,allBin,allSlice); 
//    mpi->getWorkerComm().Reduce(ninbin.data(),ninbinbuff.data(),
//                                nbin,MPI::INTEGER,MPI::SUM,0);
//    ninbin=ninbinbuff;
  }
 #endif 
  // Calculate autocorrelation function using FFT's.
 // if (workerID==0) {
    temp/=nslice;
    fftw_execute(fwd);
    temp(allSlice)=conj(temp(allSlice))*temp(allSlice);
    fftw_execute(rev);
    value -= real(temp);
    norm+=1;
 // }
}
