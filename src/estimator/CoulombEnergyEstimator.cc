#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "CoulombEnergyEstimator.h"
#include "action/Action.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "base/Paths.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include <blitz/tinyvec-et.h>

CoulombEnergyEstimator::CoulombEnergyEstimator(
  const SimulationInfo& simInfo, const Action* action, const double epsilon,
  MPIManager *mpi, const std::string& unitName, double scale, double shift)
  : ScalarEstimator("coulomb_energy","scalar-energy/coulomb-energy",
                    unitName,scale,shift),
    energy(0), etot(0), enorm(0),
    action(action), epsilon(epsilon),q(simInfo.getNPart()), mpi(mpi) {
  for (int i=0; i<q.size(); ++i) q(i)=simInfo.getPartSpecies(i).charge;

}

void CoulombEnergyEstimator::initCalc(const int nslice, const int firstSlice) {
  energy=0;
}

void CoulombEnergyEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  const SuperCell& cell=paths.getSuperCell();
  for (int jpart=0; jpart<ipart; ++jpart) {
    Vec delta=end-paths(jpart,islice);
    cell.pbc(delta);
    energy+=q(ipart)*q(jpart)/(epsilon*sqrt(dot(delta,delta)));
  }
}

void CoulombEnergyEstimator::endCalc(const int lnslice) {
  int nslice=lnslice;
  #ifdef ENABLE_MPI
  if (mpi) {
    double buffer; int ibuffer;
    mpi->getWorkerComm().Reduce(&energy,&buffer,1,MPI::DOUBLE,MPI::SUM,0);
    mpi->getWorkerComm().Reduce(&lnslice,&ibuffer,1,MPI::INT,MPI::SUM,0);
    energy=buffer; nslice=ibuffer;
  }
  #endif
  energy/=nslice;
  etot+=energy; enorm+=1;
}
