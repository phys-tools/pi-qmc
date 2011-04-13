// $Id: CollectiveMover.cc $
/*  Copyright (C) 2011 John B. Shumway, Jr.

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
#include "stats/MPIManager.h"
#include "CollectiveMover.h"
#include "Paths.h"
#include "Beads.h"
#include "RandomNumGenerator.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include <cmath>

#include <sstream>
#include <string>

#define DGESV_F77 F77_FUNC(dgesv,DGESV)
extern "C" void DGESV_F77(const int *n, const int *nrhs,
                          const double *a, const int *lda, int *ipiv,
                          double *b, const int *ldb, int *info);


CollectiveMover::CollectiveMover(Paths& paths, const int iFirstSlice,
                 const Vec &kvec, const Vec &amp, const Vec &center,
                 const IVec &randomize, const MPIManager* mpi)
  : UniformMover(amp,mpi), paths(paths), kvec(kvec), amp(amp),
    iFirstSlice(iFirstSlice), center(center), randomize(randomize),
    amplitude(amp), phase(0.), value(0.) {
  std::cout<<"CollectiveMover: Amplitude = " << amp
           <<"\n  k = "<<kvec<<std::endl;
  std::cout << "center=" << center << std::endl; 
  std::cout << "phase=" << center << std::endl; 
  /*  Test during development.
  Vec r(0.5,1.0);
  std::cout << r << std::endl;
  calcShift(r); 
  calcJacobian(r);
  double jacob = jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0);
  std::cout << "delta, Jacobian=" << value << jacob << std::endl;
  calcInverseShift(r); 
  jacob = 1./(jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0));
  std::cout << "back delta, Jacobian=" << value << jacob << std::endl;

  r = Vec(0.,2.0);
  std::cout << r << std::endl;
  calcShift(r); 
  calcJacobian(r);
  jacob = jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0);
  std::cout << "delta, Jacobian=" << value << jacob << std::endl;
  calcInverseShift(r); 
  jacob = 1./(jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0));
  std::cout << "back delta, Jacobian=" << value << jacob << std::endl; */
}

CollectiveMover::~CollectiveMover() {
 }

double CollectiveMover::makeMove(VArray& displacement, 
    const IArray& movingIndex) const {
  double tranProb = 1.;
  for (int idim=0; idim<NDIM; ++idim){
    amplitude[idim] = 2*amp[idim]*(RandomNumGenerator::getRand()-0.5);
  }
  phase = 0.;
  for (int idim=0; idim<NDIM; ++idim){
    if (randomize[idim]) phase[idim] = 6.2831853071795862
                                  * RandomNumGenerator::getRand();
  }
  bool forward = (RandomNumGenerator::getRand()>0.5);
  int npart = displacement.size();
  for (int ipart=0;ipart<npart;++ipart) {
    Vec r = paths(ipart,iFirstSlice);
    if (forward) calcShift(r); else calcInverseShift(r);
    displacement(ipart) = value;
#if NDIM==1
    double jacob = jacobian(0,0);
#elif NDIM==2
    double jacob = jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0);
#elif NDIM==3
    double jacob = jacobian(0,0)
                     *(jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1))
                  -jacobian(0,1)
                     *(jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))
                  +jacobian(0,2)
                     *(jacobian(1,0)*jacobian(2,1)-jacobian(1,1)*jacobian(2,0));
#endif
    tranProb *= (forward ? jacob : 1./jacob); 
  }
  return log(tranProb); 
}

void CollectiveMover::calcShift(const Vec &r) const {
  double v = 1.;
  for (int idim=0; idim<NDIM; ++idim) {
    v *= cos(kvec[idim] * (r[idim]-center[idim]) + phase[idim]);
  }
  value = amplitude*v;
}

void CollectiveMover::calcJacobian(const Vec &r) const {
  jacobian = 0.;
  for (int idim=0; idim<NDIM; ++idim) {
    for (int jdim=0; jdim<NDIM; ++jdim) {
      jacobian(jdim,idim) = amplitude[idim];
      for (int kdim=0; kdim<NDIM; ++kdim) {
        jacobian(jdim,idim) *= (kdim==jdim) 
         ? (-kvec[kdim] * sin(kvec[kdim]*(r[kdim]-center[kdim])+phase[kdim]))
         : (              cos(kvec[kdim]*(r[kdim]-center[kdim])+phase[kdim]));
      }
    }
  }
  for (int idim=0; idim<NDIM; ++idim) jacobian(idim,idim) += 1.;
}

void CollectiveMover::calcInverseShift(const Vec &r) const {
  calcShift(r);
  Vec delta = -value;
  calcShift(r+delta);
  calcJacobian(r+delta);
  double error2 = dot(value+delta,value+delta);
  int niter = 0;
  while (error2 > 1e-26) { 
    // Solve linear equations for Newton's method. 
    ++niter;
    if (niter>20) {
      std::cout << "ERROR: Too many Newton iterations in CollectiveMover!"
                << std::endl;
      std::cout<<"CollectiveMover: Amplitude = " << amp
               <<"\n  k = "<<kvec<<std::endl;
      std::cout << "center=" << center << std::endl; 
      std::cout << "phase=" << phase << std::endl; 
      std::cout << "amplitude=" << amplitude << std::endl; 
      std::cout << "delta=" << delta << std::endl; 
      exit(-1);
    }
    const int N=NDIM, ONE=1;
    int info;
    IVec ipiv;
    value += delta;
    DGESV_F77(&N,&ONE,jacobian.data(),&N,ipiv.data(),value.data(),&N,&info);
    delta -= value;
    calcShift(r+delta);
    calcJacobian(r+delta);
    error2 = dot(value+delta,value+delta);
  }
  value = delta; 
}
