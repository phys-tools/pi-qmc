//$Id$
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
#include <cstdlib>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include "SHONodes.h"
#include "PeriodicGaussian.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "SpinModelState.h"

#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int*, const int*, double*, const int*,
                           const int*, int*);
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
extern "C" void DGETRI_F77(const int*, double*, const int*, const int*,
                           double*, const int*, int*);

SHONodes::SHONodes(const SimulationInfo &simInfo,
  const Species &species, const double omega, const double temperature,
  const int maxlevel)
  : tau(simInfo.getTau()),temperature(temperature),mass(species.mass),
    omega(omega), coshwt(cosh(omega*0.5/temperature)),
    sinhwt(sinh(omega*0.5/temperature)), c(mass*0.5*omega/sinhwt),
    npart(species.count),ifirst(species.ifirst), 
    matrix((int)(pow(2,maxlevel)+0.1)+1),
    ipiv(npart),lwork(npart*npart),work(lwork),
    notMySpecies(false), gradArray(npart), gradMatrix(npart,npart), 
    mat2(npart,npart), grad2Matrix(npart,npart), nerror(0), scale(1.) {
  for (unsigned int i=0; i<matrix.size(); ++i)  {
    matrix[i] = new Matrix(npart,npart,ColMajor());
  }
  std::cout << "SHONodes with temperature = " << temperature
            << " for species " << species.name 
            << " and omega=" << omega << std::endl;
}

SHONodes::~SHONodes() {
}

NodeModel::DetWithFlag
SHONodes::evaluate(const VArray &r1, const VArray &r2, 
                   const int islice, bool scaleMagnitude) {
  DetWithFlag result; result.err=false;
  do { // Loop if scale is wrong to avoid overflow/underflow.
    Matrix& mat(*matrix[islice]);
    mat=0;
    for(int jpart=0; jpart<npart; ++jpart) {
      Vec rj=r1(jpart+ifirst);
      double rj2=dot(rj,rj);
      for(int ipart=0; ipart<npart; ++ipart) {
        Vec ri=r2(ipart+ifirst);
        double ri2=dot(ri,ri);
        double rirj=dot(ri,rj);
        mat(ipart,jpart)= scale*exp(-c*((ri2+rj2)*coshwt-2.0*rirj));
      }
    }
    if (hasSpinModelState()) {
      SpinModelState::IArray spin=(spinModelState->getSpinState());
      for(int jpart=0; jpart<npart; ++jpart) {
        for(int ipart=0; ipart<npart; ++ipart) {
          if (spin(ipart) != spin(jpart)) mat(ipart,jpart) = 0.;
        }
      }
    }
    // Calculate determinant and inverse.
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
    result.det = det;
  }
  while (scaleMagnitude && (result.err ||
         (fabs(result.det)<1e-50 || fabs(result.det)>1e50)));
  return result;
}

void SHONodes::evaluateDistance(const VArray &r1, const VArray &r2,
                     const int islice, Array& d1, Array& d2) {
  const Matrix& mat(*matrix[islice]);
  double d2inv=0; 
  // Calculate the nodal distance from jpart particles.
  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0;
    Vec rj=r1(jpart+ifirst);
    double rj2=dot(rj,rj);
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec ri=r2(ipart+ifirst);
      double ri2=dot(ri,ri);
      double rirj=dot(ri,rj);
      //Vec grad=2*c*ri*exp(-c*(ri2*coshwt-2.0*rirj));
      Vec grad=scale*2*c*(ri-coshwt*rj)*exp(-c*((ri2+rj2)*coshwt-2*rirj));
      logGrad+=mat(jpart,ipart)*grad;
    }
    d2inv+=dot(logGrad,logGrad);
  }
  // And calculate the nodal distance from ipart particles.
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0;
    Vec ri=r2(ipart+ifirst);
    double ri2=dot(ri,ri);
    for(int jpart=0; jpart<npart; ++jpart) {
      Vec rj=r1(jpart+ifirst);
      double rj2=dot(rj,rj);
      double rirj=dot(ri,rj);
      //Vec grad=2*c*rj*exp(-c*(rj2*coshwt-2.0*rirj));
      Vec grad=scale*2*c*(rj-coshwt*ri)*exp(-c*((ri2+rj2)*coshwt-2*rirj));
      logGrad+=mat(jpart,ipart)*grad;
    }
    d2inv+=dot(logGrad,logGrad);
  }
  //return sqrt(2*mass/(d2inv*tau));
}

void SHONodes::evaluateGradLogDist(const VArray &r1, const VArray &r2,
          const int islice, VMatrix &gradd1, VMatrix& gradd2,
          const Array &dist1, const Array& dist2) {
/*  const Matrix& mat(*matrix[islice]);
  gradArray=0;
  gradMatrix=0;
  grad2Matrix=Mat(0);
  double coshwt=cosh(omega*0.5/temperature);
  double sinhwt=sinh(omega*0.5/temperature);
  for (int jpart=0; jpart<npart; ++jpart) {
    for (int ipart=0; ipart<npart; ++ipart) {
      // First calculate first derivative (loggrad) terms.
      Vec grad=-mass*omega*(r1(jpart+ifirst)*coshwt-r2(ipart+ifirst))/sinhwt;
      gradArray(jpart)+=mat(jpart,ipart)*grad;
      // Calculate more complicated second derivative terms.
      for (int kpart=0; kpart<npart; ++kpart) {
        if (kpart==jpart) {
          Mat grad2;
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              if (i==j) {
                grad2(i,j)=-mass*omega*coshwt/sinhwt;
              } else {
                grad2(i,j)=-mass*omega/sinhwt
                            *(r1(jpart+ifirst)[i]*coshwt-r2(ipart+ifirst)[i])
                            *(r1(jpart+ifirst)[j]*coshwt-r2(ipart+ifirst)[j]);
              }
            }
          }
          Vec rj=r1(jpart+ifirst);
          Vec ri=r2(ipart+ifirst);
          double ri2=dot(ri,ri);
          double rj2=dot(rj,rj);
          double rirj=dot(ri,rj);
          double delta2=dot(rj-ri,rj-ri);
          grad2*= exp( (mass*omega*(ri2+rj2)*coshwt-2*rirj)/(2.0*sinhwt)
                       - temperature*mass*delta2 );
          grad2Matrix(jpart,jpart)+=mat(jpart,ipart)*grad2;
        } else {
          gradMatrix(jpart,kpart)+=mat(kpart,ipart)*grad;
        }
      }
    }
    force(jpart)=gradArray(jpart);
  }
  for (int jpart=0; jpart<npart; ++jpart) {
    for (int kpart=0; kpart<npart; ++kpart) {
      if (kpart!=jpart) {
        for (int i=0; i<NDIM; ++i) {
          for (int j=0; j<NDIM; ++j) {
            grad2Matrix(jpart,kpart)(i,j)
              =gradArray(jpart)[i]*gradArray(kpart)[j]
              -gradMatrix(jpart,kpart)[i]*gradMatrix(kpart,jpart)[j];
          }
        }
      }
    }
  }
  for (int jpart=0; jpart<npart; ++jpart) {
    for (int kpart=0; kpart<npart; ++kpart) {
      force(jpart)-=(0.5*tau/mass)*dist*dist*grad2Matrix(jpart,kpart)*g(kpart);
    }
  }
  // Calculate terms with gradients on other beads.
  for (int jpart=0; jpart<npart; ++jpart) {
    gradArray=0;
    for (int kpart=0; kpart<npart; ++kpart) {
      Vec grad=-mass*omega*(r2(jpart+ifirst)*coshwt-r1(kpart+ifirst))/sinhwt;
      for (int ipart=0; ipart<npart; ++ipart) {
        gradArray(ipart)+=mat(kpart,ipart)*grad;
      }
    }
    for (int kpart=0; kpart<npart; ++kpart) {
      for (int i=0; i<NDIM; ++i) {
        mat2(kpart,jpart)[i]=mat(kpart,jpart)/gradArray(jpart)[i];
      }
    }
    for (int ipart=0; ipart<npart; ++ipart) {
      if (ipart!=jpart) {
        for (int kpart=0; kpart<npart; ++kpart) {
          mat2(kpart,ipart)=mat(kpart,ipart);
          for (int i=0; i<NDIM; ++i) {
            mat2(kpart,ipart)[i]-=mat2(kpart,jpart)[i]*gradArray(ipart)[i];
          }
        }
      }
    }
    // Now second derivative.
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec temp=0;
      for (int kpart=0; kpart<npart; ++kpart) {
        if (jpart==kpart) {
          Mat grad2;
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              if (i==j) {
                grad2(i,j)=-mass*omega*coshwt/sinhwt;
              } else {
                grad2(i,j)=-mass*omega/sinhwt
                            *(r1(ipart+ifirst)[i]*coshwt-r2(kpart+ifirst)[i])
                            *(r1(ipart+ifirst)[j]*coshwt-r2(kpart+ifirst)[j]);
              }
            }
          }
          Vec rj=r1(jpart+ifirst);
          Vec ri=r2(ipart+ifirst);
          double ri2=dot(ri,ri);
          double rj2=dot(rj,rj);
          double rirj=dot(ri,rj);
          double delta2=dot(rj-ri,rj-ri);
          grad2*= exp( (mass*omega*(ri2+rj2)*coshwt-2*rirj)/(2.0*sinhwt)
                       - temperature*mass*delta2 );
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              temp[i] += grad2(i,j)*mat2(ipart,kpart)[j]
                                   *gradArray(jpart)[j]*gradArray(jpart)[j];
            }
          }
        } else {
          Vec grad=-mass*omega*(r2(jpart+ifirst)*coshwt-r1(kpart+ifirst))
                   /sinhwt;
          for (int j=0; j<NDIM; ++j) {
            temp+=grad*mat2(ipart,kpart)[j]
                      *gradArray(jpart)[j]*gradArray(jpart)[j];
          }
        }
      }
      force(ipart)-=(0.5*tau/mass)*dist*dist*temp;
    }
  } */
} 

/* int main(int argc, char** argv) {
  const double SCALE_RAND=1./RAND_MAX;
  const double EPS=1e-4;
  std::cout << "Test free particle nodes" << std::endl;
  const int npart=(int)pow(2,NDIM);
  Species *species=new Species("e",npart,1.0,-1.0,2,true);  
  species->ifirst=0;
  std::vector<Species*> speciesList(1), speciesIndex(npart);
  speciesList[0]=species;
  for (int i=0; i<npart; ++i) speciesIndex[i]=species;
  const double temperature=0.01, tau=0.01;
  const int nslice = (int)(1.0/(temperature*tau));
  SuperCell *cell=new SuperCell(SuperCell::Vec(10.));
  cell->computeRecipricalVectors();
  SimulationInfo simInfo(cell,npart,speciesList,speciesIndex,temperature,
                         tau,nslice);
  SHONodes nodes(simInfo, *species, temperature, 4);
  // Put particles on a grid with random displacements.
  SHONodes::VArray r1(npart), r2(npart);
  for (int ipart=0; ipart<npart; ++ipart) {
    for (int idim=0; idim<NDIM; ++idim) {
      int i=(ipart/(int)pow(2,idim))&1;
      r1(ipart)[idim]=r2(ipart)[idim]=2.5-5.0*i;
      r1(ipart)[idim]+=1.0*(SCALE_RAND*rand()-0.5);
      r2(ipart)[idim]+=1.0*(SCALE_RAND*rand()-0.5);
    }
  }
  SHONodes::Vec delta=r1(0)-r2(0);
  nodes.evaluate(r1,r2,0);
  double d=nodes.evaluateDistance(r1,r2,0);
  std::cout << "Distance to node " << d << std::endl;
  // Now check forces.
  SHONodes::VArray f(npart);
  nodes.evaluateGradLogDist(r1,r2,0,f,d);
  for (int ipart=0; ipart<npart; ++ipart) {
    for (int idim=0; idim<NDIM; ++idim) {
      std::cout << ipart << "[" << idim << "] = " << f(ipart)[idim] << " ";
      r1(ipart)[idim]+=EPS;
      nodes.evaluate(r1,r2,0);
      double dp=nodes.evaluateDistance(r1,r2,0);
      r1(ipart)[idim]-=2*EPS;
      nodes.evaluate(r1,r2,0);
      double dm=nodes.evaluateDistance(r1,r2,0);
      r1(ipart)[idim]+=EPS;
      std::cout << "(numerical derivative: " 
                << (dp-dm)/(2*EPS*d) << ")" << std::endl;
    }
  }
  return 0;
}*/
