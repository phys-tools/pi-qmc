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
#include "FreePartNodesNoUpdate.h"
#include "sampler/SectionSamplerInterface.h"
#include "sampler/DoubleSectionChooser.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Paths.h"
#include "Beads.h"
#include "SuperCell.h"
#include <iostream>
#include <blitz/tinyvec-et.h>

#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int*, const int*, double*, const int*,
                           const int*, int*);
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
extern "C" void DGETRI_F77(const int*, double*, const int*, const int*,
                           double*, const int*, int*);

FreePartNodesNoUpdate::FreePartNodesNoUpdate(const SimulationInfo& simInfo,
  const Species& species, const double temperature, const int maxlevel)
  : tau(simInfo.getTau()),mass(species.mass),npart(species.count),
    ifirst(species.ifirst), 
    r1(npart), r2(npart), matrix((int)(pow(2,maxlevel)+0.1)+1),
    det((int)(pow(2,maxlevel)+0.1)+1),
    newDet((int)(pow(2,maxlevel)+0.1)+1),
    dist((int)(pow(2,maxlevel)+0.1)+1),
    newDist((int)(pow(2,maxlevel+1)+0.1)+1),
    ipiv(npart),lwork(npart*npart),work(lwork),
    cell(*simInfo.getSuperCell()), pg(mass*temperature,cell.a[0],500),
    notMySpecies(false), force(npart), gradArray(npart), 
    gradMatrix(npart,npart), grad2Matrix(npart,npart) {
  for (unsigned int i=0; i<matrix.size(); ++i) { 
    matrix[i]=new Matrix(npart,npart,ColMajor());
  }
  std::cout << "FreePartNodesNoUpdate with temperature = "
            << temperature << std::endl;
}

FreePartNodesNoUpdate::~FreePartNodesNoUpdate() {
  for (unsigned int i=0; i<matrix.size(); ++i) delete matrix[i];
}

void FreePartNodesNoUpdate::initialize(const DoubleSectionChooser& chooser) {
  const Beads<NDIM>& sectionBeads1=chooser.getBeads(1);
  const Beads<NDIM>& sectionBeads2=chooser.getBeads(2);
  nslice=sectionBeads1.getNSlice();
  for (int islice=0; islice<nslice; ++islice) {
    for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i+ifirst,islice);
    for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i+ifirst,islice);
    det(islice)=evaluate(r1,r2,islice);
    if (det(islice)*det(0)<=0.0)
          std::cout << "ERROR - det changed signs " << islice << std::endl;
    dist(islice)=evaluateDistance(islice);
  } 
  newDet(0)=det(0); newDist(0)=dist(0);
} 

double FreePartNodesNoUpdate::evaluateDistance(const int islice) const {
  const Matrix& mat(*matrix[islice]);
  double d2inv=0; 
  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0;
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec delta=r1(jpart)-r2(ipart);
      cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=pg.grad(fabs(delta[i]))/pg(fabs(delta[i]));
        if (delta[i]<0) grad[i]=-grad[i];
      }
      for (int i=0; i<NDIM; ++i) grad*=pg(fabs(delta[i]));
      logGrad+=mat(jpart,ipart)*grad;
    }
    d2inv+=dot(logGrad,logGrad);
  }
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0;
    for(int jpart=0; jpart<npart; ++jpart) {
      Vec delta=r1(jpart)-r2(ipart);
      cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=pg.grad(fabs(delta[i]))/pg(fabs(delta[i]));
        if (delta[i]<0) grad[i]=-grad[i];
      }
      for (int i=0; i<NDIM; ++i) grad*=pg(fabs(delta[i]));
      logGrad+=mat(jpart,ipart)*grad;
    }
    d2inv+=dot(logGrad,logGrad);
  }
  return 1./sqrt(d2inv); 
}

double FreePartNodesNoUpdate::getActionDifference(
    const SectionSamplerInterface& sampler, int level) {
  // Get ready to move paths.
  double deltaAction=0;
  const Beads<NDIM>& sectionBeads1=sampler.getSectionBeads(1);
  const Beads<NDIM>& sectionBeads2=sampler.getSectionBeads(2);
  const Beads<NDIM>& movingBeads1=sampler.getMovingBeads(1);
  const Beads<NDIM>& movingBeads2=sampler.getMovingBeads(2);
  const IArray& index1=sampler.getMovingIndex(1); 
  const IArray& index2=sampler.getMovingIndex(2); 
  if (index1(0)<ifirst||index1(0)>=ifirst+npart) {notMySpecies=true; return 0;}
  else notMySpecies=false;
  int nSlice=sectionBeads1.getNSlice();
  const int nMoving=index1.size();
  const int nStride=(int)pow(2,level+1);
  // First check for node crossing.
  for (int islice=nStride/2; islice<nSlice; islice+=nStride) {
    for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i+ifirst,islice);
    for (int i=0; i<nMoving; ++i) r1(index1(i)-ifirst)=movingBeads1(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i+ifirst,islice);
    if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i)
                                  r2(index2(i)-ifirst)=movingBeads2(i,islice);
    newDet(islice)=evaluate(r1,r2,islice);
    if (newDet(islice)*det(0)<=0) return deltaAction=2e100;
  } 
  // Calculate the nodal action if level=0;
  if (level==0) {
    const double twomovert=2*mass/tau;
    double d0=0,newd0=0,d1=0,newd1=0;
    for (int islice=0;islice<nSlice; ++islice) { 
      newd1=d1=dist(islice);
      if (islice>0) {
        for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i+ifirst,islice);
        for (int i=0; i<nMoving; ++i)
                                   r1(index1(i)-ifirst)=movingBeads1(i,islice);
        for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i+ifirst,islice);
        if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i)
                                   r2(index2(i)-ifirst)=movingBeads2(i,islice);
        newd1=newDist(islice)=evaluateDistance(islice);
        deltaAction+=log( (1-exp(-twomovert*d1*d0))
                         /(1-exp(-twomovert*newd1*newd0)) );
//        deltaAction+=-log(1-exp(-twomovert*newd1*newd0))
//                     +log(1-exp(-twomovert*d1*d0));
	//std::cout << newd1 << " " << newd0 << " " << d1 << " " << d0 << std::endl;
      }
      d0=d1; newd0=newd1; 
    }
  }
  return deltaAction;
}

double FreePartNodesNoUpdate::evaluate(const VArray& r1, const VArray& r2,
    const int islice) const {
  Matrix& mat(*matrix[islice]);
  mat=0;
  for(int jpart=0; jpart<npart; ++jpart) {
    for(int ipart=0; ipart<npart; ++ipart) {
      Vec delta(r1(jpart)-r2(ipart));
      cell.pbc(delta);
      double ear2=1;
      for (int i=0; i<NDIM; ++i) ear2*=pg(fabs(delta[i]));
      mat(ipart,jpart)=ear2;
    }
  }
  // Calculate determinant and inverse.
  int info=0;//LU decomposition
  DGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
  if (info!=0) std::cout << "BAD RETURN FROM ZGETRF!!!!" << std::endl;
  double det = 1;
  for (int i=0; i<npart; ++i) {
    det*= mat(i,i); 
    det *= (i+1==ipiv(i))?1:-1;
  }
  DGETRI_F77(&npart,mat.data(),&npart,ipiv.data(),work.data(),&lwork,&info);
  if (info!=0) std::cout << "BAD RETURN FROM ZGETRI!!!!" << std::endl;
  return det;
} 

double FreePartNodesNoUpdate::getTotalAction(const Paths&, const int level)
    const {return 0;}

void FreePartNodesNoUpdate::getBeadAction(const Paths& paths, int iPart,
    int islice, double& u, double& utau, double& ulambda, 
    Vec& fm, Vec& fp) const {
  int totNSlice=paths.getNSlice();
  fm=0; fp=0; u=utau=ulambda=0;
  // Attribute u and utau to first particle.
  // We only calculate determinants when iPart==ifirst, then store
  // the forces in the force array.
  if (iPart==ifirst) {
    // Calculate the action and the gradient of the action.
    // Calculate d_i-1
    int jslice=(islice+totNSlice/2)%totNSlice;
    for (int i=0; i<npart; ++i) r1(i)=paths(i+ifirst,islice,-1);
    for (int i=0; i<npart; ++i) r2(i)=paths(i+ifirst,jslice,-1);
    evaluate(r1, r2, 0);
    double dim1=evaluateDistance(0);
    // Calculate d_i+1
    for (int i=0; i<npart; ++i) r1(i)=paths(i+ifirst,islice,+1);
    for (int i=0; i<npart; ++i) r2(i)=paths(i+ifirst,jslice,+1);
    evaluate(r1, r2, 0);
    double dip1=evaluateDistance(0);
    // Calculate the action and the gradient of the action.
    for (int i=0; i<npart; ++i) r1(i)=paths(i+ifirst,islice);
    for (int i=0; i<npart; ++i) r2(i)=paths(i+ifirst,jslice);
    evaluate(r1, r2, 0);
    double di=evaluateDistance(0);
    // Now calculate vector to node to get derivatives of action.
    const Matrix& mat(*matrix[0]);
    gradArray=0.0;
    gradMatrix=0.0;
    for (int jpart=0; jpart<npart; ++jpart) {
      for (int ipart=0; ipart<npart; ++ipart) {
        grad2Matrix(ipart,jpart)=0.;
        // First calculate first derivative (loggrad) terms.
        Vec delta=r1(jpart)-r2(ipart);
        cell.pbc(delta);
        Vec grad;
        for (int i=0; i<NDIM; ++i) {
          grad[i]=pg.grad(fabs(delta[i]))/pg(fabs(delta[i]));
          if (delta[i]<0) grad[i]=-grad[i];
        }
        for (int i=0; i<NDIM; ++i) grad*=pg(fabs(delta[i]));
        gradArray(jpart)+=mat(jpart,ipart)*grad;
        // Calculate more complicated second derivative terms.
        for (int kpart=0; kpart<npart; ++kpart) {
          if (kpart==jpart) {
            Mat grad2;
            for (int i=0; i<NDIM; ++i) {
              for (int j=0; j<NDIM; ++j) {
                if (i==j) {
                  grad2(i,i)=pg.d2(fabs(delta[i]))/pg(fabs(delta[i]));
                } else {
                  grad2(i,j)=pg.grad(fabs(delta[i]))/pg(fabs(delta[i])) 
                            *pg.grad(fabs(delta[j]))/pg(fabs(delta[j]));
                  if (delta[i]*delta[j]<0) grad2(i,j)=-grad2(i,j);
                }
              }
            }
            for (int ii=0; ii<NDIM; ++ii) {
              for (int i=0; i<NDIM; ++i) {
                for (int j=0; j<NDIM; ++j) {
                  grad2(i,j) *= pg(fabs(delta[i]));
                }
              }
            }
            for (int i=0; i<NDIM; ++i) {
              for (int j=0; j<NDIM; ++j) {
                grad2Matrix(jpart,jpart)(i,j) += mat(jpart,ipart)*grad2(i,j);
              }
            }
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
        for (int i=0; i<NDIM; ++i) {
          for (int j=0; i<NDIM; ++i) {
            force(jpart)(i) -= di*di*grad2Matrix(jpart,kpart)(i,j)*gradArray(kpart)(j);
          }
        }
      }
    }
    // Calculate terms with gradients on other beads.
    for (int jpart=0; jpart<npart; ++jpart) {
      gradArray=0.0;
      for (int kpart=0; kpart<npart; ++kpart) {
        Vec delta=r2(jpart)-r1(kpart);
        cell.pbc(delta);
        Vec grad;
        for (int i=0; i<NDIM; ++i) {
          grad[i]=pg.grad(fabs(delta[i]))/pg(fabs(delta[i]));
          if (delta[i]<0) grad[i]=-grad[i];
        }
        for (int i=0; i<NDIM; ++i) grad*=pg(fabs(delta[i]));
        for (int ipart=0; ipart<npart; ++ipart) {
          gradArray(ipart)+=mat(kpart,ipart)*grad;
        }
      }
      blitz::Array<Vec,2> mat2(npart);
      for (int kpart=0; kpart<npart; ++kpart) {
        mat2(kpart,jpart)=mat(kpart,jpart)*gradArray(jpart)
                          /dot(gradArray(jpart),gradArray(jpart));
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
        Vec temp=0.0;
        for (int kpart=0; kpart<npart; ++kpart) {
          Vec delta=r1(ipart)-r2(kpart);
          cell.pbc(delta);
          if (jpart==kpart) {
            Mat grad2;
            for (int i=0; i<NDIM; ++i) {
              for (int j=0; j<NDIM; ++j) {
                if (i==j) {
                  grad2=-pg.d2(fabs(delta[i]))/pg(fabs(delta[i]));
                } else {
                  grad2=-pg.grad(fabs(delta[i]))/pg(fabs(delta[i]))
                        *pg.grad(fabs(delta[j]))/pg(fabs(delta[j]));
                  if (delta[i]*delta[j]<0) grad2(i,j)=-grad2(i,j);
                }
              }
            }
            for (int ii=0; ii<NDIM; ++ii) {
              for (int i=0; i<NDIM; ++i) {
                for (int j=0; j<NDIM; ++j) {
                  grad2(i,j) *= pg(fabs(delta[ii]));
                }
              }
            }
            for (int i=0; i<NDIM; ++i) {
              for (int j=0; j<NDIM; ++j) {
                temp(i) +=mat(ipart,jpart)*grad2(i,j)*gradArray(jpart)(j);
              }
            }
          } else {
            Vec grad;
            for (int i=0; i<NDIM; ++i) {
              grad[i]=pg.grad(fabs(delta[i]))/pg(fabs(delta[i]));
              if (delta[i]<0) grad[i]=-grad[i];
            }
            for (int i=0; i<NDIM; ++i) grad*=pg(fabs(delta[i]));
            temp+=grad*dot(mat2(ipart,kpart),gradArray(jpart));
          }
          force(ipart)-=di*di*temp;
        }
      }
    }
    // Now calulate factor on force.
    double xip1=2*mass*dip1*di/tau;
    double xim1=2*mass*dim1*di/tau;
    double factor=xip1*exp(-xip1)/(1-exp(-xip1))
                 +xim1*exp(-xim1)/(1-exp(-xim1));
    for (int jpart=0; jpart<npart; ++jpart) force(jpart)*=factor;
    // Calculate the nodal action.
    u = -0.5*log(1-exp(-xim1));
    utau = 0.5*xim1*exp(-xim1)/(tau*(1-exp(-xim1)));
  }
  if (iPart>=ifirst && iPart-ifirst<npart) fm=force(iPart-ifirst);
}
void FreePartNodesNoUpdate::acceptLastMove() {
  if (notMySpecies) return;
  for (int i=0; i<nslice; ++i) {
    det(i)=newDet(i);
    dist(i)=newDist(i);
  }
}
