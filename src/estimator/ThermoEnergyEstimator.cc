#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "ThermoEnergyEstimator.h"
#include "action/Action.h"
#include "action/DoubleAction.h"
#include "base/SimulationInfo.h"
#include "stats/MPIManager.h"
#include <cstdlib>
#include <blitz/tinyvec.h>

ThermoEnergyEstimator::ThermoEnergyEstimator(
  const SimulationInfo& simInfo, const Action* action,
  const DoubleAction* doubleAction, MPIManager *mpi,
  const std::string& unitName, double scale, double shift)
  : ScalarEstimator("thermo_energy","scalar-energy/thermo-energy",
                    unitName,scale,shift),
    energy(0), etot(0), enorm(0), action(action), doubleAction(doubleAction),
    mpi(mpi) {
}

void ThermoEnergyEstimator::initCalc(const int nslice, const int firstSlice) {
  energy=0;
}

void ThermoEnergyEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {

  if (action) {
    double u(0),utau(0),ulambda(0);
    Vec fm,fp;
    action->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
    energy+=utau;
  }
 
  if (doubleAction) {
    double u(0),utau(0),ulambda(0);
    Vec fm,fp;
    doubleAction->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
    energy+=utau;
  }

}

void ThermoEnergyEstimator::endCalc(const int lnslice) {
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
