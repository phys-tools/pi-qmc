// $Id$
/*  Copyright (C) 2004-2007 John B. Shumway, Jr.

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
#include "CoulombEnergyEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Action.h"
#include "Paths.h"
#include "SuperCell.h"
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"


CoulombEnergyEstimator::CoulombEnergyEstimator(
  const SimulationInfo& simInfo, const Action* action, const double epsilon,
  MPIManager *mpi, const std::string& unitName, double scale, double shift)
  : ScalarEstimator("coulomb_energy",unitName,scale,shift),
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
