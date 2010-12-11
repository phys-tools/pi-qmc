//$Id$
/*  Copyright (C) 2009,2010 John B. Shumway, Jr.

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
#include <blitz/tinyvec-et.h>
#include "AugmentedNodes.h"
#include "PeriodicGaussian.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Beads.h"
#include "DoubleMLSampler.h"
#include <cstdlib>

#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int*, const int*, double*, const int*,
                           const int*, int*);
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
extern "C" void DGETRI_F77(const int*, double*, const int*, const int*,
                           double*, const int*, int*);
#define ASSNDX_F77 F77_FUNC(assndx,ASSNDX)
extern "C" void ASSNDX_F77(const int *mode, double *a, const int *n, 
  const int *m, const int *ida, int *k, double *sum, int *iw, const int *idw);

AugmentedNodes::AugmentedNodes(const SimulationInfo &simInfo,
  const Species &species, const double temperature,
  const int maxlevel, const bool useUpdates, const int maxMovers,
  double density, const std::vector<const AtomicOrbitalDM*>& orbitals,
  bool useHungarian)
  : NodeModel("_"+species.name),
    tau(simInfo.getTau()),mass(species.mass),
    npart(species.count), ifirst(species.ifirst), 
    matrix((int)(pow(2,maxlevel)+0.1)+1),
    ipiv(npart),lwork(npart*npart),work(lwork),
    cell(*simInfo.getSuperCell()),
    pg(NDIM), pgp(NDIM), pgm(NDIM),
    notMySpecies(false),
    gradArray1(npart), gradArray2(npart), 
    temp1(simInfo.getNPart()), temp2(simInfo.getNPart()),
    uarray(npart,npart,ColMajor()), 
    kindex((int)(pow(2,maxlevel)+0.1)+1,npart), kwork(npart*6), nerror(0),
    useHungarian(useHungarian), scale(1.0),
    density(density), orbitals(orbitals) {
  for (unsigned int i=0; i<matrix.size(); ++i)  {
    matrix[i] = new Matrix(npart,npart,ColMajor());
  }
  if (density==0.) density = pow(mass*temperature/PI,0.5*NDIM);
  std::cout << "AugmentedNodes with temperature = "
            << temperature << std::endl;
  std::cout << "density = " << density << std::endl;
  double tempp=temperature/(1.0+EPSILON); //Larger beta (plus).
  double tempm=temperature/(1.0-EPSILON); //Smaller beta (minus).
  if (temperature!=simInfo.getTemperature()) {
    tempp=tempm=temperature;
  }
  for (int idim=0; idim<NDIM; ++idim) {
    pg[idim]=new PeriodicGaussian(mass*temperature,cell.a[idim],
                            (int)(100*cell.a[idim]*sqrt(mass*temperature)));
    pgm[idim]=new PeriodicGaussian(mass*tempm,cell.a[idim],
                             (int)(100*cell.a[idim]*sqrt(mass*tempm)));
    pgp[idim]=new PeriodicGaussian(mass*tempp,cell.a[idim],
                             (int)(100*cell.a[idim]*sqrt(mass*tempp)));
  }
  if (useUpdates) {
    updateObj = new MatrixUpdate(maxMovers,maxlevel,npart,matrix,*this);
  }
}

AugmentedNodes::~AugmentedNodes() {
  for (int idim=0; idim<NDIM; ++idim) {
    delete pg[idim]; delete pgm[idim]; delete pgp[idim];
  }
  delete updateObj;
}

NodeModel::DetWithFlag 
AugmentedNodes::evaluate(const VArray &r1, const VArray &r2, 
                         const int islice, bool scaleMagnitude) {
  DetWithFlag result; result.err=false;
  Matrix& mat(*matrix[islice]);
  do { // Loop if scale is wrong to avoid overflow/underflow.
    mat=0.;
    for(int jpart=0; jpart<npart; ++jpart) {
      for(int ipart=0; ipart<npart; ++ipart) {
        Vec delta(r1(jpart+ifirst)-r2(ipart+ifirst));
        cell.pbc(delta);
        double ear2=scale*density;
        for (int i=0; i<NDIM; ++i) ear2 *= (*pg[i])(fabs(delta[i]));
        mat(ipart,jpart) = ear2;
      }
    }
    // Add contribution from orbitals.
    for (std::vector<const AtomicOrbitalDM*>::const_iterator orb 
         = orbitals.begin(); orb != orbitals.end(); ++orb) {
      int kfirst=(*orb)->ifirst, nkpart=(*orb)->npart;
      VArray2& delta1((*orb)->getD1Array());
      for (int jpart=0; jpart<npart; ++jpart) {
        for (int kpart=0; kpart<nkpart; ++kpart) {
          Vec delta(r1(jpart+ifirst)-r1(kpart+kfirst));
          cell.pbc(delta);
          delta1(jpart,kpart)=delta;
        }
      }
      VArray2& delta2((*orb)->getD2Array());
      for(int ipart=0; ipart<npart; ++ipart) {
        for (int kpart=0; kpart<nkpart; ++kpart) {
          Vec delta(r2(ipart+ifirst)-r2(kpart+kfirst));
          cell.pbc(delta);
          delta2(ipart,kpart)=delta;
        }
      }
      (*orb)->evaluateValue(mat,scale);
    }
    for(int jpart=0; jpart<npart; ++jpart) {
      for(int ipart=0; ipart<npart; ++ipart) {
        uarray(ipart,jpart)=-log(fabs(mat(ipart,jpart))+1e-100);
      }
    }
    // Find dominant contribution to determinant (distroys uarray).
    if (useHungarian) {
      const int MODE=1; //"Case 1", sum uarray(i,j=k(i)) minimized.
      double usum=0;
      ASSNDX_F77(&MODE,uarray.data(),&npart,&npart,&npart,&kindex(islice,0),
               &usum,kwork.data(),&npart);
      for (int ipart=0; ipart<npart; ++ipart) kindex(islice,ipart)-=1;
      // Note: u(ipart,jpart=kindex(islice,ipart)) makes maximum contribution
      // or lowest total action.
    }
    // Next calculate determinant and inverse of slater matrix.
    int info=0;//LU decomposition
    DGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
    if (info!=0) {
      result.err = true;
      std::cout << "BAD RETURN FROM ZGETRF!?!!" << std::endl;
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
      std::cout << "BAD RETURN FROM ZGETRI!?!!" << std::endl;
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

void AugmentedNodes::evaluateDotDistance(const VArray &r1, const VArray &r2,
         const int islice, Array &d1, Array &d2) {
  std::vector<PeriodicGaussian*> pgSave(NDIM);
  for (int i=0; i<NDIM; ++i) pgSave[i]=pg[i];
  double tauSave=tau;
  // Calculate gradient with finite differences.
  // First calculate distance at smaller beta.
  tau = tauSave*(1.-EPSILON);
  for (int i=0; i<NDIM; ++i) pg[i]=pgm[i];
  evaluate(r1, r2, islice, false);
  evaluateDistance(r1, r2, islice, temp1, temp2);
  // Next calculate distance at larger beta.
  tau = tauSave*(1.+EPSILON);
  for (int i=0; i<NDIM; ++i) pg[i]=pgp[i];
  evaluate(r1, r2, islice, false);
  evaluateDistance(r1, r2, islice, d1, d2);
  // Now use finite differences.
  double denom=1./(2*tau*EPSILON);
  d1 = denom*(d1-temp1);
  d2 = denom*(d2-temp2);
  // Restore pg and tau to proper temperature.
  for (int i=0; i<NDIM; ++i) pg[i]=pgSave[i];
  tau=tauSave;
}

void AugmentedNodes::evaluateDistance(const VArray& r1, const VArray& r2,
                              const int islice, Array& d1, Array& d2) {
  Matrix& mat(*matrix[islice]);
  // Calculate log gradients to estimate distance.
  d1=200; d2=200; // Initialize distances to a very large value.
  // Set up distances for orbitals.
  for (std::vector<const AtomicOrbitalDM*>::const_iterator orb 
       = orbitals.begin(); orb != orbitals.end(); ++orb) {
    int kfirst=(*orb)->ifirst, nkpart=(*orb)->npart;
    VArray2& delta1((*orb)->getD1Array());
    for (int jpart=0; jpart<npart; ++jpart) {
      for (int kpart=0; kpart<nkpart; ++kpart) {
        Vec delta(r1(jpart+ifirst)-r1(kpart+kfirst));
        cell.pbc(delta);
        delta1(jpart,kpart)=delta;
      }
    }
    VArray2& delta2((*orb)->getD2Array());
    for(int ipart=0; ipart<npart; ++ipart) {
      for (int kpart=0; kpart<nkpart; ++kpart) {
        Vec delta(r2(ipart+ifirst)-r2(kpart+kfirst));
        cell.pbc(delta);
        delta2(ipart,kpart)=delta;
      }
    }
    (*orb)->evaluateValueAndGrad(0);
  }
  // Calculate gradients for slice 1.
  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0, fgrad=0.0;
    for (int ipart=0; ipart<npart; ++ipart) {
      AtomicOrbitalDM::ValAndGrad orbital;
      for (std::vector<const AtomicOrbitalDM*>::const_iterator orb 
          = orbitals.begin(); orb != orbitals.end(); ++orb) {
        AtomicOrbitalDM::ValAndGrad temp =  (**orb)(ipart,jpart);
        orbital += temp;
      }
      Vec delta=r1(jpart+ifirst)-r2(ipart+ifirst);
      cell.pbc(delta);
      Vec grad; double val=density;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=(*pg[i]).grad(fabs(delta[i]))/(*pg[i])(fabs(delta[i])+1e-300);
        val *= (*pg[i])(fabs(delta[i]));
        if (delta[i]<0) grad[i]=-grad[i];
      }
      grad *= val;
      grad += orbital.grad1;
      val += orbital.val;
      if (useHungarian && jpart==kindex(islice,ipart)) fgrad=grad/(val+1e-300);
      logGrad += mat(jpart,ipart)*grad*scale;
    }
    gradArray1(jpart) = logGrad - fgrad;
    d1(jpart+ifirst) = sqrt(2*mass/
                      ((dot(gradArray1(jpart),gradArray1(jpart))+1e-15)*tau));
  }
  // Calculate gradients for slice 2.
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0, fgrad=0.0;
    for (int jpart=0; jpart<npart; ++jpart) {
      AtomicOrbitalDM::ValAndGrad orbital;
      for (std::vector<const AtomicOrbitalDM*>::const_iterator orb 
          = orbitals.begin(); orb != orbitals.end(); ++orb) {
        AtomicOrbitalDM::ValAndGrad temp =  (**orb)(ipart,jpart);
        orbital += temp;
      }
      Vec delta=r2(ipart+ifirst)-r1(jpart+ifirst);
      cell.pbc(delta);
      Vec grad; double val=density;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=(*pg[i]).grad(fabs(delta[i]))/(*pg[i])(fabs(delta[i])+1e-300);
        val *= (*pg[i])(fabs(delta[i]));
        if (delta[i]<0) grad[i]=-grad[i];
      }
      grad *= val;
      grad += orbital.grad2;
      val += orbital.val;
      if (useHungarian && jpart==kindex(islice,ipart)) fgrad=grad/(val+1e-300);
      logGrad += mat(jpart,ipart)*grad*scale;
    }
    gradArray2(ipart) = logGrad - fgrad;
    d2(ipart+ifirst) = sqrt(2*mass/
                      ((dot(gradArray2(ipart),gradArray2(ipart))+1e-15)*tau));
  }
}

void AugmentedNodes::evaluateGradLogDist(const VArray &r1, const VArray &r2,
       const int islice, VMatrix &gradd1, VMatrix &gradd2, 
       const Array& dist1, const Array& dist2) {
  gradd1=0.; gradd2=0.; 
  const Matrix& mat(*matrix[islice]);
  // ipart is the index of the particle in the gradient.
  // jpart is the index of the particle in the distance.
  for (int ipart=0; ipart<npart; ++ipart) {
    for (int jpart=0; jpart<npart; ++jpart) {
      // First treat d1 terms.
      if (ipart==jpart) {
        Mat d2ii;
        d2ii = 0.0;
        for(int kpart=0; kpart<npart; ++kpart) {
          Vec delta=r1(ipart+ifirst)-r2(kpart+ifirst);
          cell.pbc(delta);
          Mat d2iik;
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              if (i==j) {
                d2iik(i,j)=(*pg[i]).d2(fabs(delta[i]))
                          /(*pg[i])(fabs(delta[i]));
              } else {
                d2iik(i,j)=(*pg[i]).grad(fabs(delta[i]))
                          /(*pg[i])(fabs(delta[i])) 
                          *(*pg[j]).grad(fabs(delta[j]))
                          /(*pg[j])(fabs(delta[j]));
                if (delta[i]*delta[j]<0) d2iik(i,j)=-d2iik(i,j);
              }
            }
          }
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              d2ii(i,j) += mat(ipart,kpart)*d2iik(i,j);
            }
          }
        }
        for (int i=0; i<NDIM; ++i) {
          for (int j=0; j<NDIM; ++j) {
            gradd1(ipart+ifirst,jpart+ifirst)(i)=gradArray1(jpart)(j)*d2ii(j,i);
          }
        }
      } else { // ipart!=jpart
        // Two derivatives act on different columns, so you get a 2x2 det.
        Vec g00(0.0),g01(0.0),g10(0.0),g11(0.0);
        for(int kpart=0; kpart<npart; ++kpart) {
          // First derivative acts on particle ipart.
          Vec delta=r1(ipart+ifirst)-r2(kpart+ifirst);
          cell.pbc(delta);
          Vec grad;
          for (int i=0; i<NDIM; ++i) {
            grad[i]=(*pg[i]).grad(fabs(delta[i]))/(*pg[i])(fabs(delta[i]));
            if (delta[i]<0) grad[i]=-grad[i];
          }
          for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]));
          g00 += mat(ipart,kpart)*grad;
          g10 += mat(jpart,kpart)*grad;
          // Next derivative acts on particle jpart.
          delta=r1(jpart+ifirst)-r2(kpart+ifirst);
          cell.pbc(delta);
          for (int i=0; i<NDIM; ++i) {
            grad[i]=(*pg[i]).grad(fabs(delta[i]))/(*pg[i])(fabs(delta[i]));
            if (delta[i]<0) grad[i]=-grad[i];
          }
          for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]));
          g01 += mat(ipart,kpart)*grad;
          g11 += mat(jpart,kpart)*grad;
        }
        for (int i=0; i<NDIM; ++i) {
          for (int j=0; j<NDIM; ++j) {
            gradd1(ipart+ifirst,jpart+ifirst)(i)=gradArray1(jpart)(j)
                                 *(g00(i)*g11(j)-g10(i)*g01(j));
          }
        }
      }
      // Then treat d2 terms.
      // Derivative on row is a sum of derivatives on each column.
      for(int lpart=0; lpart<npart; ++lpart) {
        if (ipart==lpart) {
          // Only one term is non-zero. NEED TO CHECK THESE DERIVATIVES
          Vec delta=r1(ipart+ifirst)-r2(jpart+ifirst);
          cell.pbc(delta);
          Mat d2ii;
          d2ii = 0.0;
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              if (i==j) {
                d2ii(i,j)=(*pg[i]).d2(fabs(delta[i]))
                         /(*pg[i])(fabs(delta[i]));
              } else {
                d2ii(i,j)=(*pg[i]).grad(fabs(delta[i]))
                         /(*pg[i])(fabs(delta[i])) 
                         *(*pg[j]).grad(fabs(delta[j]))
                         /(*pg[j])(fabs(delta[j]));
                if (delta[i]*delta[j]>=0) d2ii(i,j)=-d2ii(i,j);
              }
            }
          }
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
               d2ii(i,j) *= mat(ipart,jpart);
            }
          }
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              gradd2(ipart+ifirst,jpart+ifirst)(i)=gradArray2(jpart)(j)*d2ii(j,i);
            }
          }
        } else { // ipart!=lpart
          // Two derivatives act on different columns, so you get a 2x2 det.
          Vec g00(0.0),g01(0.0),g10(0.0),g11(0.0);
          for(int kpart=0; kpart<npart; ++kpart) {
            // First derivative acts on particle ipart.
            Vec delta=r1(ipart+ifirst)-r2(kpart+ifirst);
            cell.pbc(delta);
            Vec grad;
            for (int i=0; i<NDIM; ++i) {
              grad[i]=(*pg[i]).grad(fabs(delta[i]))/(*pg[i])(fabs(delta[i]));
              if (delta[i]<0) grad[i]=-grad[i];
            }
            for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]));
            g00 += mat(ipart,kpart)*grad;
            g10 += mat(lpart,kpart)*grad;
          }
          // Next derivative acts on particle jpart.
          Vec delta=r2(jpart+ifirst)-r1(lpart+ifirst);
          cell.pbc(delta);
          Vec grad(0.0);
          for (int i=0; i<NDIM; ++i) {
            grad[i]=(*pg[i]).grad(fabs(delta[i]))/(*pg[i])(fabs(delta[i]));
            if (delta[i]<0) grad[i]=-grad[i];
          }
          for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]));
          g01 = mat(ipart,jpart)*grad;
          g11 = mat(lpart,jpart)*grad;
          for (int i=0; i<NDIM; ++i) {
            for (int j=0; j<NDIM; ++j) {
              gradd2(ipart+ifirst,jpart+ifirst)(i)=gradArray2(jpart)(j)
                                   *(g00(i)*g11(j)-g10(i)*g01(j));
            }
          }
        }
      }
      // Put in prefactors.
      gradd1(ipart+ifirst,jpart+ifirst)
        *= -(tau/mass)*dist1(jpart+ifirst)*dist1(jpart+ifirst);
      gradd2(ipart+ifirst,jpart+ifirst)
        *= -(tau/mass)*dist2(jpart+ifirst)*dist2(jpart+ifirst);
      // Then add linear (j independent) contribution.
      gradd1(ipart+ifirst,jpart+ifirst)+=gradArray1(ipart);
      gradd2(ipart+ifirst,jpart+ifirst)+=gradArray1(ipart);
    }
  }
}

const double AugmentedNodes::EPSILON=1e-6;

AugmentedNodes::MatrixUpdate::MatrixUpdate(int maxMovers, int maxlevel, 
    int npart, std::vector<Matrix*> &matrix, const AugmentedNodes &fpNodes)
  : fpNodes(fpNodes), maxMovers(maxMovers), npart(npart),
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

double AugmentedNodes::MatrixUpdate::evaluateChange(
    const DoubleMLSampler &sampler, int islice) {
  isNewMatrixUpdated=false;
  const Beads<NDIM>& sectionBeads2=sampler.getSectionBeads(2);
  const Beads<NDIM>& movingBeads1=sampler.getMovingBeads(1);
  // Get info on moving paths of this species.
  const IArray &movingIndex(sampler.getMovingIndex(1));
  nMoving=0;
  for (int i=0; i<movingIndex.size(); ++i) {
    if (movingIndex(i) >= fpNodes.ifirst &&
        movingIndex(i) < fpNodes.npart+fpNodes.ifirst) {
      index1(nMoving) = movingIndex(i);
      mindex1(nMoving) = i;
      ++nMoving;
    }
  }
  // Compute new slater matrix elements for moving particles.
  for (int jmoving=0; jmoving<nMoving; ++jmoving) {
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec delta(movingBeads1(mindex1(jmoving),islice)
               -sectionBeads2(ipart+fpNodes.ifirst,islice));
      fpNodes.cell.pbc(delta);
      double ear2=1;
      for (int i=0; i<NDIM; ++i) ear2*=(*fpNodes.pg[i])(fabs(delta[i]));
      (*phi[islice])(ipart,jmoving)=ear2;
    }
    for (int imoving=0; imoving<nMoving; ++imoving) {
      int ipart=index1(imoving)-fpNodes.ifirst;
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
      int ipart=index1(imoving)-fpNodes.ifirst;
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

void AugmentedNodes::MatrixUpdate::evaluateNewInverse(const int islice) {
  *newMatrix[islice]=*matrix[islice];
  for (int jmoving=0; jmoving<nMoving; ++jmoving) {
    int jpart=index1(jmoving)-fpNodes.ifirst;
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

void AugmentedNodes::MatrixUpdate::evaluateNewDistance(
    const VArray& r1, const VArray& r2,
    const int islice, Array &d1, Array &d2) {
  // Calculate log gradients to estimate distance.
  d1=200; d2=200; // Initialize distances to a very large value.
  const Matrix& mat(*newMatrix[islice]);
  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0, fgrad=0.0;
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec delta=r1(jpart+fpNodes.ifirst)-r2(ipart+fpNodes.ifirst);
      fpNodes.cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=(*fpNodes.pg[i]).grad(fabs(delta[i]))/(*fpNodes.pg[i])(fabs(delta[i]));
        if (delta[i]<0) grad[i]=-grad[i];
      }
      if (ipart==fpNodes.kindex(islice,jpart)) fgrad=grad;
      for (int i=0; i<NDIM; ++i) grad*=(*fpNodes.pg[i])(fabs(delta[i]));
      logGrad+=mat(jpart,ipart)*grad;
    }
    fpNodes.gradArray1(jpart)=logGrad-fgrad;
    d1(jpart+fpNodes.ifirst)=sqrt(2*fpNodes.mass
      /((dot(logGrad,logGrad)+1e-15)*fpNodes.tau));
  }
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0, fgrad=0.0;
    for(int jpart=0; jpart<npart; ++jpart) {
      Vec delta=r2(ipart+fpNodes.ifirst)-r1(jpart+fpNodes.ifirst);
      fpNodes.cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=(*fpNodes.pg[i]).grad(fabs(delta[i]))/(*fpNodes.pg[i])(fabs(delta[i]));
        if (delta[i]<0) grad[i]=-grad[i];
      }
      if (ipart==fpNodes.kindex(islice,jpart)) fgrad=grad;
      for (int i=0; i<NDIM; ++i) grad*=(*fpNodes.pg[i])(fabs(delta[i]));
      logGrad+=mat(jpart,ipart)*grad;
    }
    fpNodes.gradArray2(ipart)=logGrad-fgrad;
    d2(ipart+fpNodes.ifirst)=sqrt(2*fpNodes.mass
      /((dot(logGrad,logGrad)+1e-15)*fpNodes.tau));
  }
}

void AugmentedNodes::MatrixUpdate::acceptLastMove(int nslice) {
  if (isNewMatrixUpdated) {
    for (int islice=1; islice<nslice-1; ++islice) {
      Matrix* ptr=matrix[islice];    
      matrix[islice]=newMatrix[islice];
      newMatrix[islice]=ptr;
    }
  } else {
    for (int islice=1; islice<nslice-1; ++islice) {
      for (int jmoving=0; jmoving<nMoving; ++jmoving) {
        int jpart=index1(jmoving)-fpNodes.ifirst;
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

const double AugmentedNodes::PI = acos(-1.0);

AugmentedNodes::AtomicOrbitalDM::AtomicOrbitalDM(
    int ifirst, int npart, int nfermion, double weight)
  : weight(weight), ifirst(ifirst), npart(npart),
    dist1(nfermion,npart), dist2(nfermion,npart) {
}


AugmentedNodes::Atomic1sDM::Atomic1sDM(
    double Z, int ifirst, int npart, int nfermion, double weight)
  : AtomicOrbitalDM(ifirst, npart, nfermion, weight),
    npart(npart), nfermion(nfermion), coef(sqrt(weight*Z*Z*Z/PI)),
    work1(nfermion,npart), work2(nfermion,npart), Z(Z) {
  std::cout << "1s: Z="<<Z<<", "<<weight<<", "<<ifirst << std::endl;
}


AugmentedNodes::Atomic2sDM::Atomic2sDM(
    double Z, int ifirst, int npart, double weight)
  : AtomicOrbitalDM(ifirst, npart, npart, weight), Z(Z) {
  std::cout << "2s: Z="<<Z<<", "<<weight<<", "<<ifirst << std::endl;
}

AugmentedNodes::Atomic2pDM::Atomic2pDM(
    double Z, int ifirst, int npart, double weight)
  : AtomicOrbitalDM(ifirst, npart, npart, weight), Z(Z) {
  std::cout << "2p: Z="<<Z<<", "<<weight<<", "<<ifirst << std::endl;
}

AugmentedNodes::Atomic2spDM::Atomic2spDM(
    double Z, int ifirst, int npart, double weight, double alpha)
  : AtomicOrbitalDM(ifirst, npart, npart, weight), Z(Z), alpha(alpha) {
  std::cout << "2sp: Z=" << Z << ", " << alpha 
            << weight << ", " << ifirst << std::endl;
}

AugmentedNodes::AtomicOrbitalDM::ValueAndGradient
AugmentedNodes::Atomic1sDM::operator()(double r1, double r2,
    double costheta) const {
  ValueAndGradient result;
  result.value = weight*exp(-Z*(r1+r2))*Z*Z*Z/PI;
  result.gradr1 = -Z*result.value;
  result.gradr2 = -Z*result.value;
  result.gradcostheta = 0;
  return result;
}

AugmentedNodes::AtomicOrbitalDM::ValueAndGradient
AugmentedNodes::Atomic2sDM::operator()(double r1, double r2,
    double costheta) const {
  ValueAndGradient result;
  double temp  = weight*exp(-0.5*Z*(r1+r2))*Z*Z*Z/(32*PI);
  result.value = temp*(2.-Z*r1)*(2.-Z*r2);
  result.gradr1 = 0.5*Z*(Z*r1-4.)*(2.-Z*r2);
  result.gradr2 = 0.5*Z*(Z*r2-4.)*(2.-Z*r1);
  result.gradcostheta = 0;
  return result;
}

AugmentedNodes::AtomicOrbitalDM::ValueAndGradient
AugmentedNodes::Atomic2pDM::operator()(double r1, double r2,
    double costheta) const {
  ValueAndGradient result;
  double temp = weight*exp(-0.5*Z*(r1+r2))*Z*Z*Z*Z*Z/(32*PI);
  result.value = temp*r1*r2*(costheta+1.);
  result.gradr1 = (1.-0.5*Z*r1)*r2*temp*(costheta+1.);
  result.gradr2 = (1.-0.5*Z*r2)*r1*temp*(costheta+1.);
  result.gradcostheta = temp*r1*r2;
  return result;
}

AugmentedNodes::AtomicOrbitalDM::ValueAndGradient
AugmentedNodes::Atomic2spDM::operator()(double r1, double r2,
    double costheta) const {
  ValueAndGradient result;
  double temp  = weight*exp(-0.5*Z*(r1+r2)-alpha*Z*Z*(r1*r1+r2*r2
                             -2*r1*r2*costheta))*Z*Z*Z/PI;
  result.value = temp*(1+4*alpha*Z*Z*r1*r2);
  result.gradr1 = (2*alpha*Z*Z*(r2*costheta-r1)-0.5*Z)*result.value
                  +4*Z*Z*alpha*r2*temp;
  result.gradr2 = (2*alpha*Z*Z*(r1*costheta-r2)-0.5*Z)*result.value
                  +4*Z*Z*alpha*r1*temp;
  result.gradcostheta = 2*alpha*Z*Z*r1*r2*result.value;
  return result;
}

void
AugmentedNodes::Atomic1sDM::evaluateValue(Matrix &mat, double scale) const {
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      work1(i,k) = sqrt(dot(dist1(i,k),dist1(i,k)));
      work2(i,k) = sqrt(dot(dist2(i,k),dist2(i,k)));
    }
  }
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      work1(i,k) = coef*exp(-Z*work1(i,k));
      work2(i,k) = coef*exp(-Z*work2(i,k));
    }
  }
  for (int i=0; i<nfermion; ++i) {
    for (int j=0; j<nfermion; ++j) {
      for (int k=0; k<npart; ++k) {
        mat(i,j) += scale*work1(j,k)*work2(i,k);
      }
    }
  }
}

void
AugmentedNodes::Atomic1sDM::evaluateValueAndGrad(double scale) const {
  // First set the work arrays to the scalar distances.
  for (int i=0; i<nfermion; ++i) {
    for (int j=0; j<npart; ++j) {
      work1(i,j) = sqrt(dot(dist1(i,j),dist1(i,j)));
      work2(i,j) = sqrt(dot(dist2(i,j),dist2(i,j)));
    }
  }
  // Then rescale the distance arrays to their unit vectors.
  for (int i=0; i<nfermion; ++i) {
    for (int j=0; j<npart; ++j) {
      dist1(i,j) /= (work1(i,j)+1e-300);
      dist2(i,j) /= (work2(i,j)+1e-300);
    }
  }
  // Then evaluate the orbitals.
  for (int i=0; i<nfermion; ++i) {
    for (int j=0; j<npart; ++j) {
      work1(i,j) = coef*exp(-Z*work1(i,j));
      work2(i,j) = coef*exp(-Z*work2(i,j));
    }
  }
  // Finally, set the distance arrays to the gradients.
  for (int i=0; i<nfermion; ++i) {
    for (int j=0; j<npart; ++j) {
      dist1(i,j) *= -Z*work1(i,j);
      dist2(i,j) *= -Z*work2(i,j);
    }
  }
}

AugmentedNodes::AtomicOrbitalDM::ValAndGrad
AugmentedNodes::Atomic1sDM::operator()(int i, int j) const {
  ValAndGrad result;
  result.val = 0.;
  result.grad1 = 0.;
  result.grad2 = 0.;
  for (int k=0; k<npart; ++k) {
    result.val += work1(j,k)*work2(i,k);
    result.grad1 += dist1(j,k)*work2(i,k);
    result.grad2 += work1(j,k)*dist2(i,k);
  }
  return result;
}
