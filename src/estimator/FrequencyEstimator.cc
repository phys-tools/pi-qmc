#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "FrequencyEstimator.h"
#include "action/Action.h"
#include "action/DoubleAction.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include <blitz/tinyvec-et.h>

FrequencyEstimator::FrequencyEstimator(const SimulationInfo& simInfo,
  const Species& species1, const Species& species2, 
  int nfreq, int nstride, MPIManager *mpi)
  : BlitzArrayBlkdEst<1>("frequency", "dynamic-scalar/frequency",
                         IVecN(nfreq),true),
    npart(simInfo.getNPart()), nslice(simInfo.getNSlice()), 
    nfreq(nfreq), nstride(nstride), tau(simInfo.getTau()),
    ipart(species1.ifirst), jpart(species2.ifirst+species2.count-1),
    temp(nslice/nstride), mpi(mpi) {
  std::cout << "Frequency Estimator " << species1.name 
            << " " << species2.name << std::endl;
  fftw_complex *ptr = (fftw_complex*)temp.data();
  fwd = fftw_plan_dft_1d(nslice/nstride, ptr, ptr, FFTW_FORWARD, FFTW_MEASURE);
#ifdef ENABLE_MPI
  if (mpi) mpiBuffer.resize(nslice/nstride);
#endif
}

FrequencyEstimator::~FrequencyEstimator() {
  fftw_destroy_plan(fwd);
}

void FrequencyEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0 ;
}

void FrequencyEstimator::handleLink(const Vec& start, const Vec& end,
          const int iipart, const int islice, const Paths& paths) {
  if (iipart==ipart) {
    Vec delta = start-paths(jpart,islice);
    double delta2 = dot(delta, delta);
    temp(islice/nstride) += sqrt(delta2);
  }
}

void FrequencyEstimator::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  mpiBuffer = 0.;
  if (mpi) {
    mpi->getWorkerComm().Reduce(temp.data(), mpiBuffer.data(),
                                2*(nslice/nstride),MPI::DOUBLE,MPI::SUM,0);
    temp = mpiBuffer; 
  }
#endif 
  // Calculate autocorrelation function using FFT's.
  if (workerID==0) {
    fftw_execute(fwd);
    temp *= tau;
    temp(allSlice)=conj(temp(allSlice))*temp(allSlice);
    double betainv=1./(tau*nslice);
    value -= real(temp(blitz::Range(0,nfreq-1)))*betainv;
    norm+=1;
  }
}
