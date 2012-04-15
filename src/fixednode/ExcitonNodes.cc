//$Id$
/*  Copyright (C) 2008-2010 John B. Shumway, Jr.

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
#include <cstdlib>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include "ExcitonNodes.h"
#include "util/PeriodicGaussian.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "util/SuperCell.h"
#include "Beads.h"
#include "sampler/DoubleMLSampler.h"
#include "util/Hungarian.h"

#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int*, const int*, double*, const int*,
                           const int*, int*);
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
extern "C" void DGETRI_F77(const int*, double*, const int*, const int*,
                           double*, const int*, int*);

ExcitonNodes::ExcitonNodes(const SimulationInfo &simInfo,
  const Species &species1, const Species &species2,
  const double temperature, const int maxlevel,
  const double radius, const bool useUpdates, const int maxMovers)
  : NodeModel("_"+species1.name+species2.name),
    alpha(1./radius),
    tau(simInfo.getTau()),mass(species1.mass),npart(species1.count),
    ifirst(species1.ifirst), jfirst(species2.ifirst),
    matrix((int)(pow(2,maxlevel)+0.1)+1),
    ipiv(npart),lwork(npart*npart),work(lwork),
    cell(*simInfo.getSuperCell()), pg(NDIM), pgp(NDIM), pgm(NDIM),
    notMySpecies(false),
    gradArray1(npart), gradArray2(npart), 
    temp1(simInfo.getNPart()), temp2(simInfo.getNPart()),
    uarray(npart,npart,ColMajor()), 
    kindex((int)(pow(2,maxlevel)+0.1)+1,npart), nerror(0),
    scale(1.0),
    hungarian(new Hungarian(npart)) {
    //gradArray(npart), mat2(npart,npart),
    //gradMatrix(npart,npart), grad2Matrix(npart,npart) {
  for (unsigned int i=0; i<matrix.size(); ++i)  {
    matrix[i] = new Matrix(npart,npart,ColMajor());
  }
  std::cout << "ExcitonNodes with radius = "
            << radius << std::endl;
  double tempp=temperature/(1.0+EPSILON); //Larger beta (plus).
  double tempm=temperature/(1.0-EPSILON); //Smaller beta (minus).
  for (int idim=0; idim<NDIM; ++idim) {
    pg[idim]=new PeriodicGaussian(mass*temperature,cell.a[idim],
                            (int)(100*cell.a[idim]*sqrt(mass*temperature)));
    pgm[idim]=new PeriodicGaussian(mass*tempm,cell.a[idim],
                             (int)(100*cell.a[idim]*sqrt(mass*tempm)));
    pgp[idim]=new PeriodicGaussian(mass*tempp,cell.a[idim],
                             (int)(100*cell.a[idim]*sqrt(mass*tempp)));
  }
}

ExcitonNodes::~ExcitonNodes() {
  for (int idim=0; idim<NDIM; ++idim) {
    delete pg[idim]; delete pgm[idim]; delete pgp[idim];
  }
  delete updateObj;
  delete hungarian;
}

NodeModel::DetWithFlag 
ExcitonNodes::evaluate(const VArray &r1, const VArray &r2, 
                       const int islice, bool scaleMagnitude) {
  DetWithFlag result; result.err=false;
  do { // Loop if scale is wrong to avoid overflow/underflow.
    Matrix& mat(*matrix[islice]);
    mat=0;
    for(int jpart=0; jpart<npart; ++jpart) {
      for(int ipart=0; ipart<npart; ++ipart) {
        Vec delta(r1(jpart+ifirst)-r1(ipart+jfirst));
        cell.pbc(delta);
        double ear=exp(-alpha*sqrt(dot(delta,delta)));
        mat(ipart,jpart)=scale*ear+1e-100;
        uarray(ipart,jpart)=-log(fabs(mat(ipart,jpart))+1e-100);
      }
    }
    // Find dominant contribution to determinant (distroys uarray).
    hungarian->solve(uarray.data());
    for (int ipart=0; ipart<npart; ++ipart) {
        kindex(islice,ipart) = (*hungarian)[ipart];
    }
    // Note: u(ipart,jpart=kindex(islice,ipart)) makes maximum contribution
    // or lowest total action.
    // Next calculate determinant and inverse.
    int info=0;//LU decomposition
    DGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
    if (info!=0) {
      result.err=true;
      std::cout << "BAD RETURN FROM ZGETRF!!!!" << std::endl;
      nerror++;
      if (nerror>1000) {
          std::cout << "too many errors!!!!" << std::endl;
        std::exit(-1);
      }
    }
    double det = 1;
    for (int i=0; i<npart; ++i) {
      det*= mat(i,i); 
      det *= (i+1==ipiv(i))?1:-1;
    }
    DGETRI_F77(&npart,mat.data(),&npart,ipiv.data(),work.data(),&lwork,&info);
    if (info!=0) { 
      result.err = true;
      std::cout << "BAD RETURN FROM ZGETRI!!!!" << std::endl;
      nerror++;
      if (nerror>1000) {
        std::cout << "too many errors!!!!" << std::endl;
        std::exit(-1);
      }
    }
    // watch for overflow or underflow.
    if (scaleMagnitude && (!result.err 
                           && (fabs(det)<1e-50 || fabs(det)>1e50) )) {
      if (fabs(det)<1e-250) {
        scale *= 2;
      } else {
        scale *= pow(fabs(det),-1./npart);
      }
      std::cout << "Slater determinant rescaled: scale = " 
                << scale << std::endl;
    }
    result.det = det;
  }
  while (scaleMagnitude && (result.err || 
         (fabs(result.det)<1e-50 || fabs(result.det)>1e50)));
  return result;
}

void ExcitonNodes::evaluateDotDistance(const VArray &r1, const VArray &r2,
         const int islice, Array &d1, Array &d2) {
  return;
}

void ExcitonNodes::evaluateDistance(const VArray& r1, const VArray& r2,
                              const int islice, Array& d1, Array& d2) {
  Matrix& mat(*matrix[islice]);
  // Calculate log gradients to estimate distance.
  d1=200; d2=200; // Initialize distances to a very large value.
  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0, fgrad=0.0;
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec delta(r1(jpart+ifirst)-r1(ipart+jfirst));
      cell.pbc(delta);
      double d=sqrt(dot(delta,delta));
      Vec grad = -alpha*delta/(d+1e-200);
      if (ipart==kindex(islice,jpart)) fgrad=grad;
      grad *= exp(-alpha*d);
      logGrad += mat(jpart,ipart)*scale*grad;
    }
    gradArray1(jpart)=logGrad-fgrad;
    d1(jpart+ifirst)
      =sqrt(2*mass/((dot(gradArray1(jpart),gradArray1(jpart))+1e-15)*tau));
  }
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0, fgrad=0.0;
    for (int jpart=0; jpart<npart; ++jpart) {
      Vec delta(r1(ipart+jfirst)-r1(jpart+ifirst));
      cell.pbc(delta);
      double d=sqrt(dot(delta,delta));
      Vec grad = -alpha*delta/(d+1e-200);
      if (ipart==kindex(islice,jpart)) fgrad=grad;
      grad *= exp(-alpha*d);
      logGrad += mat(ipart,jpart)*scale*grad;
    }
    gradArray1(ipart)=logGrad-fgrad;
    d1(ipart+jfirst)
      =sqrt(2*mass/((dot(gradArray1(ipart),gradArray1(ipart))+1e-15)*tau));
  }
}

void ExcitonNodes::evaluateGradLogDist(const VArray &r1, const VArray &r2,
       const int islice, VMatrix &gradd1, VMatrix &gradd2, 
       const Array& dist1, const Array& dist2) {
  return;
}

const double ExcitonNodes::EPSILON=1e-6;
