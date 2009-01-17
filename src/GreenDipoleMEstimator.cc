//$Id: GreenDipoleMEstimator.cc,v 1.7 2007/10/03 12:53:56 jshumwa Exp $
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
#include "GreenDipoleMEstimator.h"
#include "SimulationInfo.h"
#include "SuperCell.h"
#include "Species.h"
#include "Paths.h"
#include <blitz/tinyvec-et.h>
#include <iostream>
#include "stats/MPIManager.h"

GreenDipoleMEstimator::GreenDipoleMEstimator(const SimulationInfo& simInfo, const int index,
					      MPIManager *mpi)
  : BlitzArrayBlkdEst<1>(index==0?"Gdipole_moment_x" :
                         index==1?"Gdipole_moment_y" :  
		"Gdipole_moment_z" , IVecN(simInfo.getNSlice() ),true),
    npart(simInfo.getNPart()), nslice(simInfo.getNSlice()),
    tau(simInfo.getTau()), q(simInfo.getNPart()),
    index(index), names(simInfo.getNPart()), 
    temp(nslice), mpi(mpi) {
  for (int i=0; i<names.size(); ++i)
  names(i)=simInfo.getPartSpecies(i).name;
  for (int i=0; i<q.size(); ++i) q(i)=simInfo.getPartSpecies(i).charge;
  std::cout << "Green's Dipole Moment Estimator component = "
            << (index==0 ? "x" : index==1 ? "y" : "z") << std::endl;
  fftw_complex *ptr = (fftw_complex*)temp.data();
  fwd = fftw_plan_dft_1d(nslice, ptr, ptr, FFTW_FORWARD, FFTW_MEASURE);
  rev = fftw_plan_dft_1d(nslice, ptr, ptr, FFTW_BACKWARD, FFTW_MEASURE);
}

GreenDipoleMEstimator::~GreenDipoleMEstimator() {
  fftw_destroy_plan(fwd);
  fftw_destroy_plan(rev);
}

void GreenDipoleMEstimator::initCalc(const int nslice, const int firstSlice) {
  dipole=0, temp=0;
}

void GreenDipoleMEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
          temp(islice)+=q(ipart)*paths(ipart,islice)[index];
}

void GreenDipoleMEstimator::endCalc(const int nslice) {
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
 //  }
}
