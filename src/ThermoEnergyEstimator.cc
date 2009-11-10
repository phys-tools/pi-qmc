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
#include "ThermoEnergyEstimator.h"
#include "SimulationInfo.h"
#include "Action.h"
#include "DoubleAction.h"
#include <blitz/tinyvec.h>
#include "stats/MPIManager.h"

ThermoEnergyEstimator::ThermoEnergyEstimator(
  const SimulationInfo& simInfo, const Action* action,
  const DoubleAction* doubleAction, MPIManager *mpi,
  const std::string& unitName, double scale, double shift)
  : ScalarEstimator("thermo_energy",unitName,scale,shift),
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
    energy+=utau;//std :: cout << "action TE :: "<<energy<<"  "<<ipart<<" "<<islice<<"  "<<utau<<"  "<<u<<std ::endl;
  }
  if (doubleAction) {
    double u(0),utau(0),ulambda(0);
    Vec fm,fp;
    doubleAction->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
    energy+=utau;//std :: cout << "daction TE :: "<<energy<<"  "<<ipart<<" "<<islice<<"  "<<utau<<"  "<<u<<std ::endl;
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
  etot+=energy; enorm+=1; //std :: cout <<"End calc ::"<<energy<<" "<<nslice<<"  "<<etot<<std ::endl;
}
