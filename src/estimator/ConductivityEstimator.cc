#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "ConductivityEstimator.h"
#include "action/Action.h"
#include "action/DoubleAction.h"
#include "base/SimulationInfo.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include <cstdlib>
#include <blitz/tinyvec.h>

ConductivityEstimator::ConductivityEstimator(const SimulationInfo& simInfo,
  const int nfreq, const int nbin, const int ndbin, const int nstride,
  MPIManager *mpi)
  : BlitzArrayBlkdEst<3>("conductivity","dynamic-array/conductivity",
                         IVecN(2*ndbin-1,nbin,nfreq),true), 
    npart(simInfo.getNPart()), nslice(simInfo.getNSlice()), 
    nfreq(nfreq), nbin(nbin), ndbin(ndbin), nstride(nstride),
    tau(simInfo.getTau()),
    tauinv(1./tau), massinv(1./simInfo.getSpecies(0).mass),
    dx(simInfo.getSuperCell()->a[0]/nbin), dxinv(1/dx),q(npart),
    temp(2*ndbin-1,nbin,nslice/nstride), mpi(mpi) {
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

ConductivityEstimator::~ConductivityEstimator() {
  fftw_destroy_plan(rev);
  fftw_destroy_plan(fwd);
}

void ConductivityEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0;
}

void ConductivityEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  //Vec delta=paths.delta(ipart,islice,-1);
  //int ibin = (int)((end[0]-0.5*delta[0]+1.5*nbin*dx)*dxinv)%nbin;
  //temp(0,ibin,islice) += q(ipart) * delta[0] * tauinv;
  //++ninbin(ibin);
  //if (islice%nstride==0) {
    int ibin=((int)(end[0]*dxinv+nbin))%nbin;
    int jbin=((int)(start[0]*dxinv+nbin))%nbin;
    if (ibin!=jbin) {
      int nstep = ((ibin-jbin+3*nbin/2)%nbin)-nbin/2;
      int idir = (nstep>0)?1:-1;
      if (idir>0) {
        for (int i=1; i<=nstep; i++) {
          temp(0,(jbin+i)%nbin,islice/nstride)+=idir*q(ipart);
        }
      } else {
        for (int i=1; i<=-nstep; i++) {
          temp(0,(ibin+i)%nbin,islice/nstride)+=idir*q(ipart);
        }
      }
    }
  //}
}

void ConductivityEstimator::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
  blitz::Range allBin = blitz::Range::all();
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Reduce(&temp(0,0,0),&temp(1,0,0),
                                2*nslice/nstride*nbin,MPI::DOUBLE,MPI::SUM,0);
    temp(0,allBin,allSlice)=temp(1,allBin,allSlice); 
  }
#endif
  // Calculate autocorrelation function using FFT's.
  if (workerID==0) {
    //temp/=nslice;
    fftw_execute(fwd);
    //temp/=tau;
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
    //fftw_execute(rev);
    value += real(temp(allBin,allBin,blitz::Range(0,nfreq-1)))
            /(tau*nslice);
    norm+=1;
  }
}
