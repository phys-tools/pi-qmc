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
#include "EwaldCoulombEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Action.h"
#include "Paths.h"
#include "SuperCell.h"
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"
#include "TradEwaldSum.h"
#include "OptEwaldSum.h"


EwaldCoulombEstimator::EwaldCoulombEstimator(
  const SimulationInfo& simInfo, const Action* action, const double epsilon,
  const double rcut, const double kcut,
  MPIManager *mpi, const std::string& unitName, double scale, double shift)
  : ScalarEstimator("coulomb_energy",unitName,scale,shift),
    //ewaldSum(*new TradEwaldSum(*simInfo.getSuperCell(),
    //                            simInfo.getNPart(),rcut,kcut)),
    ewaldSum(*new OptEwaldSum(*simInfo.getSuperCell(),
                               simInfo.getNPart(),rcut,kcut,4*kcut,8)),
    cell(*simInfo.getSuperCell()),
    energy(0), etot(0), enorm(0), vgrid(1001), nradial(1001),
    rcut(rcut), dr(rcut/1000), drinv(1./dr),
    action(action), epsilon(epsilon),q(simInfo.getNPart()),
    r(simInfo.getNPart()), mpi(mpi) {
  for (int i=0; i<q.size(); ++i) q(i)=simInfo.getPartSpecies(i).charge;
  ewaldSum.getQArray() = q;
  ewaldSum.evalSelfEnergy();
  vgrid(0) = -ewaldSum.evalFR0()/epsilon;
  for (int i=1; i<nradial; ++i) {
    vgrid(i) = -ewaldSum.evalFR(i*dr)/epsilon;
  }

}

EwaldCoulombEstimator::~EwaldCoulombEstimator() {
  delete &ewaldSum;
}

void EwaldCoulombEstimator::initCalc(const int nslice, const int firstSlice) {
  energy=0;
}

void EwaldCoulombEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {


  for (int jpart=0; jpart<ipart; ++jpart) {
    Vec delta=end-paths(jpart,islice);
    cell.pbc(delta);
    double r=sqrt(dot(delta,delta));
    if (r<rcut) {
      int igrid=(int)(r*drinv);
      double x=r-igrid*dr;
      energy+=q(ipart)*q(jpart)
             *(1./(r*epsilon) + (1-x)*vgrid(igrid)+x*vgrid(igrid+1));
    }
  }

 
  // Add long range contribution.
  if (ipart==0) {
    paths.getSlice(islice,r);
    energy += ewaldSum.evalLongRange(r)/epsilon;
   }

}

void EwaldCoulombEstimator::endCalc(const int lnslice) {
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
