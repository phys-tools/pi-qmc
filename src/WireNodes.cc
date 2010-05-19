//$Id$
/*  Copyright (C) 2004-2009 John B. Shumway, Jr.

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
#include "WireNodes.h"
#include "PeriodicGaussian.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Beads.h"
#include "DoubleMLSampler.h"

#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int*, const int*, double*, const int*,
                           const int*, int*);
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
extern "C" void DGETRI_F77(const int*, double*, const int*, const int*,
                           double*, const int*, int*);
#define ASSNDX_F77 F77_FUNC(assndx,ASSNDX)
extern "C" void ASSNDX_F77(const int *mode, double *a, const int *n, 
  const int *m, const int *ida, int *k, double *sum, int *iw, const int *idw);

WireNodes::WireNodes(const SimulationInfo &simInfo, const Species &species,
    const double omega, const double temperature, const int maxlevel,
    const bool useUpdates, const int maxMovers)
  : NodeModel("_"+species.name),
    tau(simInfo.getTau()),mass(species.mass),npart(species.count),
    ifirst(species.ifirst), 
    matrix((int)(pow(2,maxlevel)+0.1)+1),
    ipiv(npart),lwork(npart*npart),work(lwork),
    cell(*simInfo.getSuperCell()), pg(0), pgp(0), pgm(0),
    omega(omega), coshwt(cosh(omega*0.5/temperature)),
    sinhwt(sinh(omega*0.5/temperature)), c(mass*0.5*omega/sinhwt),
    notMySpecies(false),
    gradArray1(npart), gradArray2(npart), 
    temp1(simInfo.getNPart()), temp2(simInfo.getNPart()), 
    uarray(npart,npart,ColMajor()),
    kindex((int)(pow(2,maxlevel)+0.1)+1,npart), kwork(npart*6), nerror(0),
    scale(1.0) {
  for (unsigned int islice=0; islice<matrix.size(); ++islice)  {
    matrix[islice] = new Matrix(npart,npart,ColMajor());
  }
  for (int islice=0; islice<kindex.shape()(0); ++islice)  {
    for (int ipart=0; ipart<npart; ++ipart) kindex(islice,ipart)=ipart;
  }
  std::cout << "WireNodes with temperature = "
            << temperature << std::endl;
  double tempp=temperature/(1.0+EPSILON); //Larger beta (plus).
  double tempm=temperature/(1.0-EPSILON); //Smaller beta (minus).
  pg=new PeriodicGaussian(mass*temperature,cell.a[0],
                            (int)(100*cell.a[0]*sqrt(mass*temperature)));
  pgm=new PeriodicGaussian(mass*tempm,cell.a[0],
                             (int)(100*cell.a[0]*sqrt(mass*tempm)));
  pgp=new PeriodicGaussian(mass*tempp,cell.a[0],
                             (int)(100*cell.a[0]*sqrt(mass*tempp)));
  if (useUpdates) {
    updateObj = new MatrixUpdate(maxMovers,maxlevel,npart,matrix,*this);
  }
}

WireNodes::~WireNodes() {
  delete pg; delete pgm; delete pgp;
  delete updateObj;
}

NodeModel::DetWithFlag
WireNodes::evaluate(const VArray &r1, const VArray &r2, const int islice,
    bool scaleMagnitude) {
  DetWithFlag result; result.err=false;
  do { // Loop if scale is wrong to avoid overflow/underflow.
    Matrix& mat(*matrix[islice]);
    mat=0;
    for(int jpart=0; jpart<npart; ++jpart) {
      Vec rj=r1(jpart+ifirst);
      double rj2=0; for (int i=1; i<NDIM; ++i) rj2+=rj[i]*rj[i];
      for(int ipart=0; ipart<npart; ++ipart) {
        // Free particle in idim=0.
        Vec delta(r1(jpart+ifirst)-r2(ipart+ifirst));
        cell.pbc(delta);
        double ear2=(*pg)(fabs(delta[0]));
        mat(ipart,jpart)=scale*ear2+1e-100;
        // SHO in idim>0.
        Vec ri=r2(ipart+ifirst);
        double ri2=0, rirj=0;
        for (int i=1; i<NDIM; ++i) {
          ri2+=ri[i]*ri[i];
          rirj+=ri[i]*rj[i];
        }
        mat(ipart,jpart)*= exp(-c*((ri2+rj2)*coshwt-2.0*rirj));
        uarray(ipart,jpart)=-log(fabs(mat(ipart,jpart))+1e-100);
      }
    }
    // Find dominant contribution to determinant (distroys uarray).
    const int MODE=1;
    double usum=0;
    ASSNDX_F77(&MODE,uarray.data(),&npart,&npart,&npart,&kindex(islice,0),
               &usum,kwork.data(),&npart);
    for(int ipart=0; ipart<npart; ++ipart) kindex(islice,ipart)-=1;
    // Note: u(ipart,jpart=kindex(islice,ipart)) makes maximum contribution
    // or lowest total action.
    // Next calculate determinant and inverse of slater matrix.
    int info=0;//LU decomposition
    DGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
    if (info!=0) {
      result.err = true;
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
    result.det =  det;
  }
  while (scaleMagnitude && (result.err || 
         (fabs(result.det)<1e-50 || fabs(result.det)>1e50)));
  return result;
}

void WireNodes::evaluateDotDistance(const VArray &r1, const VArray &r2,
         const int islice, Array &d1, Array &d2) {
  PeriodicGaussian* pgSave;
  pgSave=pg;
  double tauSave=tau;
  // Calculate gradient with finite differences.
  // First calculate distance at smaller beta.
  tau = tauSave*(1.-EPSILON);
  pg=pgm;
  evaluate(r1, r2, islice, false);
  evaluateDistance(r1, r2, islice, temp1, temp2);
  // Next calculate distance at larger beta.
  tau = tauSave*(1.+EPSILON);
  pg=pgp;
  evaluate(r1, r2, islice, false);
  evaluateDistance(r1, r2, islice, d1, d2);
  // Now use finite differences.
  double denom=1./(2*tau*EPSILON);
  d1 = denom*(d1-temp1);
  d2 = denom*(d2-temp2);
  // Restore pg and tau to proper temperature.
  pg=pgSave;
  tau=tauSave;
}

void WireNodes::evaluateDistance(const VArray& r1, const VArray& r2,
                              const int islice, Array& d1, Array& d2) {
  // Calculate log gradients to estimate distance.
  d1=200; d2=200; // Initialize distances to a very large value.
  const Matrix& mat(*matrix[islice]);
  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0, fgrad=0.0;
    Vec rj=r1(jpart+ifirst);
    double rj2=0; for (int i=1; i<NDIM; ++i) rj2+=rj[i]*rj[i];
    for (int ipart=0; ipart<npart; ++ipart) {
      // Free particle in idim=0.
      Vec delta=r1(jpart+ifirst)-r2(ipart+ifirst);
      cell.pbc(delta);
      Vec grad; grad=0.0;
      grad[0]=(*pg).grad(fabs(delta[0]));
      if (delta[0]<0) grad[0]=-grad[0];
      // SHO in idim>0.
      Vec ri=r2(ipart+ifirst);
      double ri2=0, rirj=0;
      for (int i=1; i<NDIM; ++i) {
        ri2+=ri[i]*ri[i];
        rirj+=ri[i]*rj[i];
      }
      double ear2=(*pg)(fabs(delta[0]));
      for (int i=1; i<NDIM; ++i) {
        grad[i]=2*c*(ri[i]-coshwt*rj[i])*ear2;
      }
      grad*=exp(-c*((ri2+rj2)*coshwt-2*rirj));
      if (ipart==kindex(islice,jpart)) {
        fgrad=grad/(ear2*exp(-c*((ri2+rj2)*coshwt-2.0*rirj)));
      }
      logGrad+=mat(jpart,ipart)*scale*grad;
    }
    gradArray1(jpart)=logGrad-fgrad;
    d1(jpart+ifirst)=sqrt(2*mass/((dot(logGrad,logGrad)+1e-15)*tau));
  }
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0, fgrad=0.0;
    Vec ri=r2(ipart+ifirst);
    double ri2=0; for (int i=1; i<NDIM; ++i) ri2+=ri[i]*ri[i];
    for(int jpart=0; jpart<npart; ++jpart) {
      // Free particle in idim=0.
      Vec delta=r2(ipart+ifirst)-r1(jpart+ifirst);
      cell.pbc(delta);
      Vec grad; grad=0.0;
      grad[0]=(*pg).grad(fabs(delta[0]));
      if (delta[0]<0) grad[0]=-grad[0];
      // SHO in idim>0.
      Vec rj=r1(jpart+ifirst);
      double rj2=0, rirj=0;
      for (int i=1; i<NDIM; ++i) {
        rj2+=rj[i]*rj[i];
        rirj+=ri[i]*rj[i];
      }
      double ear2=(*pg)(fabs(delta[0]));
      for (int i=1; i<NDIM; ++i) {
        grad[i]=2*c*(rj[i]-coshwt*ri[i])*ear2;
      }
      grad*=exp(-c*((ri2+rj2)*coshwt-2*rirj));
      if (ipart==kindex(islice,jpart)) {
        fgrad=grad/(ear2*exp(-c*((ri2+rj2)*coshwt-2.0*rirj)));
      }
      logGrad+=mat(jpart,ipart)*scale*grad;
    }
    gradArray2(ipart)=logGrad-fgrad;
    d2(ipart+ifirst)=sqrt(2*mass/((dot(logGrad,logGrad)+1e-15)*tau));
  }
}

void WireNodes::evaluateGradLogDist(const VArray &r1, const VArray &r2,
       const int islice, VMatrix &gradd1, VMatrix &gradd2, 
       const Array& dist1, const Array& dist2) {
  gradd1=0.; gradd2=0.; 
}

const double WireNodes::EPSILON=1e-6;

WireNodes::MatrixUpdate::MatrixUpdate(int maxMovers, int maxlevel, int npart,
    std::vector<Matrix*> &matrix, const WireNodes &wireNodes)
  : wireNodes(wireNodes), maxMovers(maxMovers), npart(npart),
    newMatrix((int)(pow(2,maxlevel)+0.1)+1),
    phi((int)(pow(2,maxlevel)+0.1)+1),
    bvec((int)(pow(2,maxlevel)+0.1)+1),
    smallDet(maxMovers,maxMovers),
    matrix(matrix),
    ipiv(maxMovers),lwork(maxMovers*maxMovers),work(lwork),
    isNewMatrixUpdated(false),index1(maxMovers),mindex1(maxMovers),nMoving(0) {
  for (unsigned int i=0; i<matrix.size(); ++i)  {
    newMatrix[i] = new Matrix(npart,npart,ColMajor());
    phi[i] = new Matrix(npart,maxMovers,ColMajor());
    bvec[i] = new Matrix(npart,maxMovers,ColMajor());
  }
}

double WireNodes::MatrixUpdate::evaluateChange(
    const DoubleMLSampler &sampler, int islice) {
  isNewMatrixUpdated=false;
  const Beads<NDIM>& sectionBeads2=sampler.getSectionBeads(2);
  const Beads<NDIM>& movingBeads1=sampler.getMovingBeads(1);
  // Get info on moving paths of this species.
  const IArray& movingIndex(sampler.getMovingIndex(1));
  nMoving=0;
  for (int i=0; i<movingIndex.size(); ++i) {
    if (movingIndex(i) >= wireNodes.ifirst &&
        movingIndex(i) < wireNodes.npart+wireNodes.ifirst) {
      index1(nMoving) = movingIndex(i);
      mindex1(nMoving) = i;
      ++nMoving;
    }
  }
  // Compute new Slater matrix elements for moving particles.
  for (int jmoving=0; jmoving<nMoving; ++jmoving) {
    Vec rj=movingBeads1(mindex1(jmoving),islice);
    double rj2=0; for (int i=1; i<NDIM; ++i) rj2+=rj[i]*rj[i];
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec ri=sectionBeads2(ipart+wireNodes.ifirst,islice);
      // Free particle in idim=0.
      Vec delta(rj-ri);
      wireNodes.cell.pbc(delta);
      double ear2=(*wireNodes.pg)(fabs(delta[0]));
      (*phi[islice])(ipart,jmoving)=ear2;
      // SHO in idim>0.
      double ri2=0, rirj=0;
      for (int i=1; i<NDIM; ++i) {
        ri2+=ri[i]*ri[i];
        rirj+=ri[i]*rj[i];
      }
      (*phi[islice])(ipart,jmoving)*= exp(-wireNodes.c*((ri2+rj2)
                               *wireNodes.coshwt-2.0*rirj));
    }
    for (int imoving=0; imoving<nMoving; ++imoving) {
      int ipart=index1(imoving)-wireNodes.ifirst;
      (*bvec[islice])(ipart,jmoving)=0;
      for (int k=0; k<npart; ++k) {
        (*bvec[islice])(ipart,jmoving)
          += (*phi[islice])(k,jmoving)*(*matrix[islice])(ipart,k);
      }
    }
  }
  // Compute change in the slater determinant.
  for (int jmoving=0; jmoving<nMoving; ++jmoving) {
    for (int imoving=0; imoving<nMoving; ++imoving) {
      int ipart=index1(imoving)-wireNodes.ifirst;
      smallDet(imoving,jmoving)=(*bvec[islice])(ipart,jmoving);
    }
  }
  // Calculate determinant and inverse.
  int info=0;//LU decomposition
  DGETRF_F77(&nMoving,&nMoving,smallDet.data(),&maxMovers,ipiv.data(),&info);
  if (info!=0) std::cout << "BAD RETURN FROM ZGETRF!!!!" << std::endl;
  double det = 1;
  for (int i=0; i<nMoving; ++i) {
    det *= smallDet(i,i); 
    det *= (i+1==ipiv(i))?1:-1;
  }
  return det;
}

void WireNodes::MatrixUpdate::evaluateNewInverse(const int islice) {
  *newMatrix[islice]=*matrix[islice];
  for (int jmoving=0; jmoving<nMoving; ++jmoving) {
    int jpart=index1(jmoving)-wireNodes.ifirst;
    double bjjinv=0;
    for (int k=0; k<npart; ++k) {
      bjjinv += (*phi[islice])(k,jmoving)*(*newMatrix[islice])(jpart,k);
    }
    bjjinv=1./bjjinv;
    for (int k=0; k<npart; ++k) {
      (*newMatrix[islice])(jpart,k) *= bjjinv;
    } 
    for (int ipart=0; ipart<npart; ++ipart) {
      if (ipart==jpart) continue;
      double bij=0;
      for (int k=0; k<npart; ++k) {
        bij += (*phi[islice])(k,jmoving)*(*newMatrix[islice])(ipart,k);
      }
      for (int k=0; k<npart; ++k) {
        (*newMatrix[islice])(ipart,k) -= (*newMatrix[islice])(jpart,k)*bij;
      }
    }
  }
  isNewMatrixUpdated=true;
}

void WireNodes::MatrixUpdate::evaluateNewDistance(
    const VArray& r1, const VArray& r2,
    const int islice, Array &d1, Array &d2) {
  // Then calculate log gradients to estimate distance.
  d1=200; d2=200; // Initialize distances to a very large value.
  const Matrix& mat(*newMatrix[islice]);
  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0, fgrad=0.0;
    Vec rj=r1(jpart+wireNodes.ifirst);
    double rj2=0; for (int i=1; i<NDIM; ++i) rj2+=rj[i]*rj[i];
    for (int ipart=0; ipart<npart; ++ipart) {
      // Free particle in idim=0.
      Vec delta=r1(jpart+wireNodes.ifirst)-r2(ipart+wireNodes.ifirst);
      wireNodes.cell.pbc(delta);
      Vec grad; grad=0.0;
      grad[0]=(*wireNodes.pg).grad(fabs(delta[0]));
      if (delta[0]<0) grad[0]=-grad[0];
      // SHO in idim>0.
      Vec ri=r2(ipart+wireNodes.ifirst);
      double ri2=0, rirj=0;
      for (int i=1; i<NDIM; ++i) {
        ri2+=ri[i]*ri[i];
        rirj+=ri[i]*rj[i];
      }
      double ear2=(*wireNodes.pg)(fabs(delta[0]));
      for (int i=1; i<NDIM; ++i) {
        grad[i]=2*wireNodes.c*(ri[i]-wireNodes.coshwt*rj[i])*ear2;
      }
      grad*=exp(-wireNodes.c*((ri2+rj2)*wireNodes.coshwt-2*rirj));
      if (ipart==wireNodes.kindex(islice,jpart)) {
        fgrad=grad/(ear2*exp(-wireNodes.c*
                               ((ri2+rj2)*wireNodes.coshwt-2.0*rirj)));
      }
      logGrad+=mat(jpart,ipart)*grad;
    }
    wireNodes.gradArray1(jpart)=logGrad-fgrad;
    d1(jpart+wireNodes.ifirst)=sqrt(2*wireNodes.mass
         /((dot(logGrad,logGrad)+1e-15)*wireNodes.tau));
  }
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0, fgrad=0.0;
    Vec ri=r2(ipart+wireNodes.ifirst);
    double ri2=0; for (int i=1; i<NDIM; ++i) ri2+=ri[i]*ri[i];
    for(int jpart=0; jpart<npart; ++jpart) {
      // Free particle in idim=0.
      Vec delta=r2(ipart+wireNodes.ifirst)-r1(jpart+wireNodes.ifirst);
      wireNodes.cell.pbc(delta);
      Vec grad; grad=0.0;
      grad[0]=(*wireNodes.pg).grad(fabs(delta[0]));
      if (delta[0]<0) grad[0]=-grad[0];
      // SHO in idim>0.
      Vec rj=r1(jpart+wireNodes.ifirst);
      double rj2=0, rirj=0;
      for (int i=1; i<NDIM; ++i) {
        rj2+=rj[i]*rj[i];
        rirj+=ri[i]*rj[i];
      }
      double ear2=(*wireNodes.pg)(fabs(delta[0]));
      for (int i=1; i<NDIM; ++i) {
        grad[i]=2*wireNodes.c*(rj[i]-wireNodes.coshwt*ri[i])*ear2;
      }
      grad*=exp(-wireNodes.c*((ri2+rj2)*wireNodes.coshwt-2*rirj));
      if (ipart==wireNodes.kindex(islice,jpart)) {
        fgrad=grad/(ear2*exp(-wireNodes.c*((ri2+rj2)
                             *wireNodes.coshwt-2.0*rirj)));
      }
      logGrad+=mat(jpart,ipart)*grad;
    }
    wireNodes.gradArray2(ipart)=logGrad-fgrad;
    d2(ipart+wireNodes.ifirst)=sqrt(
       2*wireNodes.mass/((dot(logGrad,logGrad)+1e-15)*wireNodes.tau));
  }
}

void WireNodes::MatrixUpdate::acceptLastMove(int nslice) {
  if (isNewMatrixUpdated) {
    for (int islice=1; islice<nslice-1; ++islice) {
      Matrix* ptr=matrix[islice];    
      matrix[islice]=newMatrix[islice];
      newMatrix[islice]=ptr;
    }
  } else {
    for (int islice=1; islice<nslice-1; ++islice) {
      for (int jmoving=0; jmoving<nMoving; ++jmoving) {
        int jpart=index1(jmoving)-wireNodes.ifirst;
        double bjjinv=0;
        for (int k=0; k<npart; ++k) {
          bjjinv += (*phi[islice])(k,jmoving)*(*matrix[islice])(jpart,k);
        }
        bjjinv=1./bjjinv;
        for (int k=0; k<npart; ++k) {
          (*matrix[islice])(jpart,k) *= bjjinv;
        } 
        for (int ipart=0; ipart<npart; ++ipart) {
          if (ipart==jpart) continue;
          double bij=0;
          for (int k=0; k<npart; ++k) {
            bij += (*phi[islice])(k,jmoving)*(*matrix[islice])(ipart,k);
          }
          for (int k=0; k<npart; ++k) {
            (*matrix[islice])(ipart,k) -= (*matrix[islice])(jpart,k)*bij;
          }
        }
      }
    }
  }
}
