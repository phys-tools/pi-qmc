 //$Id$
/*  Copyright (C) 2004-2009 John B. Shumway, Jr. and Saad A. Khairallah

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
#include "FreeParticleNodes.h"
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

FreeParticleNodes::FreeParticleNodes(const SimulationInfo &simInfo,
  const Species &species, const double temperature, const int maxlevel,
				     const bool useUpdates, const int maxMovers, 
				     const bool useHungarian, const int useIterations,
				     const double nodalFactor)
  : NodeModel("_"+species.name), maxlevel(maxlevel),
    tau(simInfo.getTau()),mass(species.mass),npart(species.count),
    ifirst(species.ifirst), 
    matrix((int)(pow(2,maxlevel)+0.1)+1), 
    romatrix((int)(pow(2,maxlevel)+0.1)+1),
    ipiv(npart),lwork(npart*npart),work(lwork),
    cell(*simInfo.getSuperCell()), pg(NDIM), pgp(NDIM), pgm(NDIM),
    notMySpecies(false),
    gradArray1(npart), gradArray2(npart), 
    temp1(simInfo.getNPart()), temp2(simInfo.getNPart()),
    uarray(npart,npart,ColMajor()), 
    kindex((int)(pow(2,maxlevel)+0.1)+1,npart), kwork(npart*6), nerror(0),
    useHungarian(useHungarian), scale(1.0), useIterations(useIterations),
    nodalFactor(nodalFactor){
  for (unsigned int i=0; i<matrix.size(); ++i)  {
    matrix[i] = new Matrix(npart,npart,ColMajor()); 
    romatrix[i] = new Matrix(npart,npart,ColMajor());
  }
  std::cout << "FreeParticleNodes with temperature = "
            << temperature << std::endl;
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

FreeParticleNodes::~FreeParticleNodes() {
  for (int idim=0; idim<NDIM; ++idim) {
    delete pg[idim]; delete pgm[idim]; delete pgp[idim];
  }
  delete updateObj;
}

NodeModel::DetWithFlag
FreeParticleNodes::evaluate(const VArray &r1, const VArray &r2, 
                            const int islice, bool scaleMagnitude) {
  DetWithFlag result; result.err=false;
  do { // Loop if scale is wrong to avoid overflow/underflow.
    Matrix& mat(*matrix[islice]);
    mat=0;
    for(int jpart=0; jpart<npart; ++jpart) {
      for(int ipart=0; ipart<npart; ++ipart) {
        Vec delta(r1(jpart+ifirst)-r2(ipart+ifirst));
        cell.pbc(delta);
        double ear2=scale;
        for (int i=0; i<NDIM; ++i) ear2*=(*pg[i])(fabs(delta[i]));
        mat(ipart,jpart)=ear2;
	(*romatrix[islice])(ipart,jpart)=ear2;
      }
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
      // Note: mat(ipart,jpart=kindex(islice,ipart)) makes maximum contribution
      // or lowest total action.
    }
    // Next calculate determinant and inverse of slater matrix.
    int info=0;//LU decomposition
    DGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
    if (info!=0) {
      result.err = true;
      std::cout << "BAD RETURN FROM ZGETRF!!!! " << islice << std::endl;
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
      std::cout << "BAD RETURN FROM ZGETRI!!!! " << islice << std::endl;
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
    result.det = det; //std :: cout <<"0- det "<<det << " for islice "<<islice<<std ::endl;
  } 
  while (scaleMagnitude && (result.err ||
         (fabs(result.det)<1e-50 || fabs(result.det)>1e50)));
  return result;
}

void FreeParticleNodes::evaluateDotDistance(const VArray &r1, const VArray &r2,
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
  double denom=1./(2*tauSave*EPSILON);
  d1 = denom*(d1-temp1);
  d2 = denom*(d2-temp2);
  // Restore pg and tau to proper temperature.
  for (int i=0; i<NDIM; ++i) pg[i]=pgSave[i];
  tau=tauSave;
}
////////////////////
void FreeParticleNodes::evaluateDistance(const VArray& r1, const VArray& r2,
                              const int islice, Array& d1, Array& d2) {
  Matrix& mat(*matrix[islice]);
  // Calculate log gradients to estimate distance.
  d1=200.; d2=200.; // Initialize distances to a very large value.

 //Array maxDist2(npart);
  //maxDist2=2000; getMaxDist2( r1, maxDist2);

  for (int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0, fgrad=0.0;
    for (int ipart=0; ipart<npart; ++ipart) {
      Vec delta=r1(jpart+ifirst)-r2(ipart+ifirst);
      cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
	grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
	if (delta[i]<0) grad[i]=-grad[i];
      }
      if (useHungarian && jpart==kindex(islice,ipart)) fgrad=grad;
      for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
      logGrad+=mat(jpart,ipart)*grad*scale;
    }
    gradArray1(jpart)=logGrad-fgrad;//

//cell.pbc( gradArray1(jpart));

   //double maxDist=sqrt(maxDist2(jpart)*.5);
    double nodalDist = sqrt(1/((dot(gradArray1(jpart),gradArray1(jpart))+1e-15)));
    //if (nodalDist> maxDist) { nodalDist=maxDist;
    // }
    // d1(jpart+ifirst)=sqrt(2*mass/
    //		  ((dot(gradArray1(jpart),gradArray1(jpart))+1e-15)*tau));
    d1(jpart+ifirst)=sqrt(2*mass/tau)*nodalDist;

  }
  //maxDist2=2000; 
  //getMaxDist2( r2, maxDist2); 
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec logGrad=0.0, fgrad=0.0;
    for(int jpart=0; jpart<npart; ++jpart) {
      Vec delta=r2(ipart+ifirst)-r1(jpart+ifirst);
      cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
	grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
	if (delta[i]<0) grad[i]=-grad[i];
      }
      if (useHungarian && jpart==kindex(islice,ipart)) fgrad = grad;
      for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
      logGrad+=mat(jpart,ipart)*grad*scale;
    }
    gradArray2(ipart)=logGrad-fgrad; //
    //double maxDist=sqrt(maxDist2(ipart+ifirst)*.5);
    double nodalDist = sqrt(1/((dot(gradArray2(ipart),gradArray2(ipart))+1e-15)));
    //if (nodalDist> maxDist) nodalDist=maxDist;

//cell.pbc( gradArray2(ipart));

    // d2(ipart+ifirst)=sqrt(2*mass/
    //		  ((dot(gradArray2(ipart),gradArray2(ipart))+1e-15)*tau));
    d2(ipart+ifirst)=sqrt(2*mass/tau)*nodalDist;
  }
  //newtonRaphson(r1, r2, islice, d1, 1);
  //newtonRaphson(r2, r1, islice, d2, 2);
}

void FreeParticleNodes::newtonRaphson(const VArray& r1, const VArray& r2, const int  islice, Array& d, int  section) {
  d=200;
  Matrix& mat(*matrix[islice]);	
  Matrix invMat(npart,npart);
  Matrix matNew(npart,npart); 
  Matrix matSav(npart,npart); 
  double det; 
  double initialDet;
 
  VArray gradArray((section==1)?gradArray1:gradArray2);
  double radius2Convergence =nodalFactor*sqrt(tau*0.5/mass);
  double dfactor=sqrt(2*mass/tau);
  int info=0; 
  // 
  int iter; double prevNodalDist;
  Vec r1jnew, r1j, logGrad, fgrad, gradf;
  Array maxDist2(npart);
  maxDist2=2000;
  getMaxDist2(r1, maxDist2);
  double nodalDist; 
  if (section==2) mat.transposeSelf(1,0);
  matSav=*romatrix[islice];
  
  IArray2 localKindex((int)(pow(2,maxlevel)+0.1)+1,npart);
  localKindex=kindex;
  Vec normal;

  if (section==2) matSav.transposeSelf(1,0);
  matNew=matSav;
  getDet(matNew, initialDet);//// redundant: gives same answer for each section. Somethingalready computed under evaluate. Make initialDet[nslice]. saves time.
  matNew=matSav;
  invMat=mat;

  for (int jpart=0; jpart<npart; ++jpart) {
        
      r1j=r1(jpart+ifirst);
      if (useIterations>0) {
	matSav=*romatrix[islice];
	if (section==2) matSav.transposeSelf(1,0);
	matNew=matSav;
	invMat=mat;
      }
      prevNodalDist=0.0;
      iter=0;
      while(iter<=useIterations){
	logGrad=0.0;
	fgrad=0.0;
        gradf=0.0;
	if (iter==0){
	  for (int ipart=0; ipart<npart; ++ipart) {
	    Vec delta=r1j-r2(ipart+ifirst);
	    cell.pbc(delta);
	    Vec grad; 
	    for (int i=0; i<NDIM; ++i) {
	      grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
	      if (delta[i]<0) grad[i]=-grad[i];
	    }
	    if (useHungarian) if( section==1 && jpart==localKindex(islice,ipart) || 
				  section==2 && ipart==localKindex(islice,jpart) ) fgrad=grad;
	    for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
	    logGrad+=invMat(jpart,ipart)*grad*scale;
	  }
	  gradArray(jpart)=logGrad-fgrad;
	  gradf=gradArray(jpart);
	} else {
	  r1j=r1(jpart+ifirst);
	  matNew=matSav;
	  invMat=matSav;
	  gradf=dot(gradArray(jpart),normal)*normal;
	}
	//////////
	
	getRnew(r1j,gradf, r1jnew);
	
	Vec dr = r1jnew-r1(jpart+ifirst);//no pbc
	nodalDist=dot(dr,dr);
	nodalDist =sqrt(nodalDist);
	
	if (nodalDist<radius2Convergence){

	  double detAtr1jnew;
	  
	  getDetAtrnew(r1jnew, jpart, r2, matNew, detAtr1jnew);  
	  matNew=matSav;
	  int signprev=(initialDet>=0)?1:-1;// can be moved out
	  int signnew=(detAtr1jnew>=0)?1:-1;//To avoid underflow.
	  if (signprev!=signnew) {
	    double dd=rootBisectionSearch(jpart, r1j, r1jnew, 
					  r2, matNew, matSav, detAtr1jnew);
	    Vec dr = r1(jpart+ifirst)-r1jnew;
	    nodalDist =sqrt(dot(dr,dr));
	    
	    iter++;
	    if (fabs(prevNodalDist-nodalDist)<10e-8) {
	      iter=useIterations+300;
	    } else {
	      if (iter <=useIterations){
		prevNodalDist=nodalDist;
		findNormAtr1jnew(r1jnew, jpart, r2, matNew, normal, localKindex, islice, info, iter, section);
		matNew=matSav;
	      }
	    }
	    
	  } else {// (signprev!=signnew) no crossing
	    iter++; 
	    if (fabs(prevNodalDist-nodalDist)<10e-8) {
	      iter=useIterations+200;
	    } else {
	      if (iter <=useIterations){
		prevNodalDist=nodalDist;
		findNormAtr1jnew(r1jnew, jpart, r2, matNew, normal, localKindex, islice, info, iter, section);
		matNew=matSav;
	      }
	    }
	  }//(signprev!=signnew)
	  /*if (nodalDist< maxDist2(jpart) && nodalDist < radius2Convergence ) {std :: cout <<"iter :: (nodalDist maxDist) "<<iter<<"  "<<nodalDist<<" "<<maxDist2(jpart"<<radius2Convergence<<std::endl;
    	plotNRPoints(r1j, r1jnew, gradf,normal,iter);
	plotRho2D(jpart, r1, r2, matSav,iter);}*/
	} else {
	  iter=useIterations+100;
	}//(nodalDist<radius2Convergence)
      }//(iter<=useIterations)
      
      //db
      /*
	if (nodalDist> maxDist2(jpart) &&nodalDist < radius2Convergence ) std :: cout <<"iter :: (nodalDist maxDist) "<<iter<<"  "<<nodalDist<<" "<<maxDist2(jpart)<<".   "<<radius2Convergence<<std::endl;
	if (sqrt(dot(gradf,gradf))<5*sqrt(tau*0.5/mass)  && nodalDist<radius2Convergence &&  radius2Convergence<maxDist2(jpart)  ){
	std :: cout <<"ultra rare case/bug "<<iter<<"  "<<nodalDist<<" "<<maxDist2(jpart)<<".   "<<radius2Convergence<<std::endl; 
	}
      */
      
      if (nodalDist> maxDist2(jpart))	nodalDist = maxDist2(jpart);
        
      d(jpart+ifirst)=dfactor*nodalDist;

  }// for loop jpart

  if (section==2) mat.transposeSelf(1,0);
}


void  FreeParticleNodes::plotNRPoints(const Vec & xprev, const Vec &Xnext, const Vec& gradf, const Vec& normal, const int &iter){

  Vec xnext= Xnext;
  cell.pbc(xnext);
  
  static int counter=0; 
  std::ofstream file; 
  // file = new std::ofstream("newtonRaphsonPts");  
  file.open("newtonRaphsonPts", std::fstream::out | std::fstream::app);
  file.precision(15);
      
  
  file<<"%number    iteration    Rprev   GradAtRprev  normalAtRnext Rnext\n";
  file<<counter<<" "<<iter<<"  ";
  for (int i=0;i<NDIM;i++)
    file<<xprev[i]<<"  ";

  double normGradLogf=dot(gradf,gradf); 
  for (int i=0;i<NDIM;i++) 
    file<<gradf[i]/(normGradLogf+1e-50)<<"  ";
 
 double norm=dot(normal,normal); 
  for (int i=0;i<NDIM;i++) 
    file<<normal[i]/(norm+1e-50)<<"  ";

  for (int i=0;i<NDIM;i++)
    file<<xnext[i]<<"  ";
  file<<std::endl;
  

  file.flush();
  file.close();

  if (counter>50)  {
    counter=0;
    std :: cout <<"Deleting newtonRaphsonPts\n"; 
    file.open("newtonRaphsonPts", std::fstream::out | std::fstream::trunc);    
    file.close();
  }
  counter++;
  
}

double FreeParticleNodes::rootBisectionSearch(const int& jpart,  const Vec & xold, Vec &midpt, const VArray & r2, 
					      Matrix &matNew, const Matrix &matInitial, double &detAtmidpt ){
  //Bisection. midpt is the root
  int signa;
  int signmidpt;
  Vec a=xold;
  Vec b=midpt; 

  double detAta;
  double dist2=dot(b-a,b-a);
  while( dist2>10e-20){
   
    midpt=a+0.5*(b-a); 

    getDetAtrnew(midpt, jpart, r2, matNew, detAtmidpt);   
    matNew=matInitial; 
    getDetAtrnew(a, jpart, r2, matNew, detAta);   
    matNew=matInitial;  
 
    signmidpt=(detAtmidpt>=0)?1:-1;
    signa=(detAta>=0)?1:-1;
    if (signa == signmidpt){
      a=midpt;
    }
    else {
      b=midpt;
    }
    dist2=dot(b-a,b-a);
  }
  
  Vec delta = midpt-xold;
  //cell.pbc(delta);/////////////sak
  return sqrt( dot(delta,delta) );
}

void FreeParticleNodes::getRnew(const Vec &rold,  const Vec &gradArray, Vec &Rnew){
  // calculates next Newton Raphson point rnew and return the distance2
  double normGradLogf=dot(gradArray,gradArray); 
  Rnew = rold-gradArray/(normGradLogf+1e-50);
  //cell.pbc(Rnew);
}

void  FreeParticleNodes::getDetAtrnew( const Vec &r1j, const int& jpart, const VArray& r2, Matrix &matNew, double &detAtr1jnew){
  //calculate determinant detAtr1jnew given the change due to  r1jnew
  Vec r1jnew=r1j; 
  cell.pbc(r1jnew);
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec delta= r1jnew-r2(ipart+ifirst);
    cell.pbc(delta);
    double ea2=scale;
    for (int i=0; i<NDIM; ++i) ea2*=((*pg[i])(fabs(delta[i])+1e-300));
    matNew(ipart,jpart)=ea2;
  }
  getDet(matNew, detAtr1jnew);
}
void  FreeParticleNodes:: findNormAtr1jnew(const Vec &r1j, const int& jpart, const VArray& r2, Matrix &matNew, Vec & normal, IArray2 &localKindex, const int &islice, int &info, int &iter, const int &section){
  Vec r1jnew=r1j; 
  cell.pbc(r1jnew);
  double det;
  //calculate determinant detAtr1jnew given the change due to  r1jnew
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec delta= r1jnew-r2(ipart+ifirst);
    cell.pbc(delta);
    double ea2=scale;
    for (int i=0; i<NDIM; ++i) ea2*=((*pg[i])(fabs(delta[i])+1e-300));
    matNew(ipart,jpart)=ea2;
  }
  getDetInvMat(matNew, det, localKindex, islice, info, iter);

  Vec fgrad;  fgrad=0;
  Vec logGrad; logGrad=0;
  for (int ipart=0; ipart<npart; ++ipart) {
    Vec delta=r1jnew-r2(ipart+ifirst);
    cell.pbc(delta);
    Vec grad; 
    for (int i=0; i<NDIM; ++i) {
      grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
      if (delta[i]<0) grad[i]=-grad[i];
    }
    if (useHungarian) if( section==1 && jpart==localKindex(islice,ipart) || 
			  section==2 && ipart==localKindex(islice,jpart) ) fgrad=grad;
    for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
    logGrad+=matNew(jpart,ipart)*grad*scale;
  }
  normal=logGrad-fgrad;
  normal=normal/sqrt(1e-200+dot(normal,normal));
}

void FreeParticleNodes:: getDet( Matrix &romat, double &det){
  // calculate determinant and inverse of slater matrix.
  int info=0;
  IArray ipiv(npart);
  DGETRF_F77(&npart,&npart,romat.data(),&npart,ipiv.data(),&info);
  
  if (info!=0) {std ::cout <<"WARNING:: info!=0 in  getDet, romat="<<romat<<std :: endl;return;}
  det = 1.0;
  for (int i=0; i<npart; ++i) {
    det*= romat(i,i);
    det *= (i+1==ipiv(i))?1:-1;
  }
}

void FreeParticleNodes:: getDetInvMat( Matrix &invMat, double &det, IArray2 &localKindex, const int &islice, int &info, int &iter){

  //Check hungarian
  if (useHungarian && iter>0) {
    Matrix uarray(npart,npart);
    for(int jpart=0; jpart<npart; ++jpart) {
      for(int ipart=0; ipart<npart; ++ipart) {
        uarray(ipart,jpart)=-log(fabs(invMat(ipart,jpart))+1e-100);
      }
    }
    IArray kwork(npart*6);
    const int MODE=1; //"Case 1", sum uarray(i,j=k(i)) minimized.
    double usum=0;
    ASSNDX_F77(&MODE,uarray.data(),&npart,&npart,&npart,&localKindex(islice,0),
	       &usum,kwork.data(),&npart);
    for (int ipart=0; ipart<npart; ++ipart) localKindex(islice,ipart)-=1;
  }
  
  // calculate determinant and inverse of slater matrix.
  info=0;//LU decomposition   
  IArray ipiv(npart);
  DGETRF_F77(&npart,&npart,invMat.data(),&npart,ipiv.data(),&info);
  if (info!=0){std ::cout <<"WARNING:: info!=0 in getDetInvMat, invMat="<<invMat<<std :: endl;return;}
  det = 1.0;
  for (int i=0; i<npart; ++i) {
    det*= invMat(i,i);
    det *= (i+1==ipiv(i))?1:-1;
  }
  
  DGETRI_F77(&npart,invMat.data(),&npart,ipiv.data(),work.data(),&lwork,&info);
  
}  
void  FreeParticleNodes::plotRho2D(const int jpart, const VArray& r1, const VArray& r2, const Matrix &matInitial, const int &iter){
  static int countFiles=0; 

  //plot for 2D the density matrix.
  Matrix matNew(npart,npart); 
  matNew=matInitial;//case0 
  std::ofstream file;     
  file.open("rhoin2D", std::fstream::out | std::ios::app);
  
  Vec r1jnew; double detAtr1jnew; 
  int grid=100; 
  double invgridx=cell.a[0]/grid;
  double invgridy=cell.a[1]/grid;
 
  //case0- contour plots
  for (int x=0; x<=grid; x++){
    r1jnew[0]=-cell.a[0]*.5 + x*invgridx;
    file <<countFiles<<"  "<<iter<<"  ";       
    for (int y=0; y<=grid; y++){
      r1jnew[1]=-cell.a[1]*.5 + y*invgridy; 
      getDetAtrnew(r1jnew, jpart, r2, matNew, detAtr1jnew);   
      matNew=matInitial; 
      file <<"     "<<detAtr1jnew;//this goes with ...
    }
    file <<std::endl;// ... this line.
  }
  
  file.flush();
  file.close();
  if (countFiles>20) {file.open("rhoin2D", std::fstream::out | std::ios::trunc);  file.close();}

  //case1- 3D plot 
  file.open("rhoin2Dv2", std::fstream::out | std::ios::app);
  for (int x=0; x<=grid; x++){
    r1jnew[0]=-cell.a[0]*.5 + x*invgridx;
    for (int y=0; y<=grid; y++){
      r1jnew[1]=-cell.a[1]*.5 + y*invgridy; 
      getDetAtrnew(r1jnew, jpart, r2, matNew, detAtr1jnew);   
      matNew=matInitial;  
      file <<countFiles<<"  "<<iter<<"  ";   
      file <<r1jnew[0]<<"    "<<r1jnew[1]<<"     "<<detAtr1jnew<<std::endl;
    }
  }
  file.flush();
  file.close();
  if (countFiles>20) {file.open("rhoin2Dv2",  std::fstream::out | std::ios::trunc);  file.close();}

  file.open("coord2Dr1",  std::fstream::out |std::ios::app);
  for (int ipart=0; ipart<npart; ++ipart) {
    file <<countFiles<<"  "<<"  "<<iter<<"  "<<ipart<<"  ";
    for(int k=0;k<NDIM; k++) 
      file<<r1(ipart+ifirst)[k] <<"   ";
    file<<std::endl; 
  }
  file.flush();
  file.close();
  if (countFiles>20) {file.open("coord2Dr1", std::fstream::out |  std::ios::trunc);  file.close();}

  file.open("coord2Dr2",  std::fstream::out |std::ios::app);
  for (int ipart=0; ipart<npart; ++ipart) {
    file <<countFiles<<"  "<<"  "<<iter<<"  "<<ipart<<"  ";
    for(int k=0;k<NDIM; k++)
      file<<r2(ipart+ifirst)[k] <<"   ";
    file<<std::endl; 
  }
  file.flush();
  file.close();  
  if (countFiles>50) {countFiles=0; file.open("coord2Dr2", std::fstream::out |  std::ios::trunc);  file.close(); std :: cout <<"Deleting all debug files. \n";}
  countFiles++;
}

void FreeParticleNodes::getMaxDist2(const VArray& r, Array& maxDist2){
  double tmpDist; 
  Vec dist;
  for (int i=0; i<npart; i++){
    for (int j=0; j<npart; j++){
      if (i!=j) {
	dist=r(i+ifirst)-r(j+ifirst);
	cell.pbc(dist);
	tmpDist=dot(dist,dist);
	if (tmpDist < maxDist2(i)){
	  maxDist2(i)=tmpDist;
	}
      }
    }
    maxDist2(i)= sqrt(maxDist2(i)*.5);
  }
}

///////////////
void FreeParticleNodes::evaluateGradLogDist(const VArray &r1, const VArray &r2,
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
            grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
            if (delta[i]<0) grad[i]=-grad[i];
          }
          for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
          g00 += mat(ipart,kpart)*grad;
          g10 += mat(jpart,kpart)*grad;
          // Next derivative acts on particle jpart.
          delta=r1(jpart+ifirst)-r2(kpart+ifirst);
          cell.pbc(delta);
          for (int i=0; i<NDIM; ++i) {
            grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
            if (delta[i]<0) grad[i]=-grad[i];
          }
          for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
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
                         /((*pg[i])(fabs(delta[i]))+1e-300);
              } else {
                d2ii(i,j)=(*pg[i]).grad(fabs(delta[i]))
                         /((*pg[i])(fabs(delta[i]))+1e-300) 
                         *(*pg[j]).grad(fabs(delta[j]))
                         /((*pg[j])(fabs(delta[j]))+1e-300);
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
              grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
              if (delta[i]<0) grad[i]=-grad[i];
            }
            for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
            g00 += mat(ipart,kpart)*grad;
            g10 += mat(lpart,kpart)*grad;
          }
          // Next derivative acts on particle jpart.
          Vec delta=r2(jpart+ifirst)-r1(lpart+ifirst);
          cell.pbc(delta);
          Vec grad(0.0);
          for (int i=0; i<NDIM; ++i) {
            grad[i]=(*pg[i]).grad(fabs(delta[i]))/((*pg[i])(fabs(delta[i]))+1e-300);
            if (delta[i]<0) grad[i]=-grad[i];
          }
          for (int i=0; i<NDIM; ++i) grad*=(*pg[i])(fabs(delta[i]))+1e-300;
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

/*int main(int argc, char** argv) {
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
  FreeParticleNodes nodes(simInfo, *species, temperature, 4);
  // Put particles on a grid with random displacements.
  FreeParticleNodes::VArray r1(npart), r2(npart);
  for (int ipart=0; ipart<npart; ++ipart) {
    for (int idim=0; idim<NDIM; ++idim) {
      int i=(ipart/(int)pow(2,idim))&1;
      r1(ipart)[idim]=r2(ipart)[idim]=2.5-5.0*i;
      r1(ipart)[idim]+=1.0*(SCALE_RAND*rand()-0.5);
      r2(ipart)[idim]+=1.0*(SCALE_RAND*rand()-0.5);
    }
  }
  FreeParticleNodes::Vec delta=r1(0.0)-r2(0.0);
  nodes.evaluate(r1,r2,0);
  double d=nodes.evaluateDistance(r1,r2,0);
  std::cout << "Distance to node " << d << std::endl;
  // Now check forces.
  FreeParticleNodes::VArray f(npart);
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

const double FreeParticleNodes::EPSILON=1e-6;//6

FreeParticleNodes::MatrixUpdate::MatrixUpdate(int maxMovers, int maxlevel, 
    int npart, std::vector<Matrix*> &matrix, const FreeParticleNodes &fpNodes)
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

double FreeParticleNodes::MatrixUpdate::evaluateChange(
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
      double ear2=fpNodes.scale;
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

void FreeParticleNodes::MatrixUpdate::evaluateNewInverse(const int islice) {
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

void FreeParticleNodes::MatrixUpdate::evaluateNewDistance(
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
        grad[i]=(*fpNodes.pg[i]).grad(fabs(delta[i]))/((*fpNodes.pg[i])(fabs(delta[i]))+1e-300);
        if (delta[i]<0) grad[i]=-grad[i];
      }
      if (fpNodes.useHungarian 
          && jpart==fpNodes.kindex(islice,ipart)) fgrad=grad;
      for (int i=0; i<NDIM; ++i) grad*=(*fpNodes.pg[i])(fabs(delta[i]))+1e-300;
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
        grad[i]=(*fpNodes.pg[i]).grad(fabs(delta[i]))/((*fpNodes.pg[i])(fabs(delta[i]))+1e-300);
        if (delta[i]<0) grad[i]=-grad[i];
      }
      if (fpNodes.useHungarian 
          && jpart==fpNodes.kindex(islice,ipart)) fgrad=grad;
      for (int i=0; i<NDIM; ++i) grad*=(*fpNodes.pg[i])(fabs(delta[i]))+1e-300;
      logGrad+=mat(jpart,ipart)*grad;
    }
    fpNodes.gradArray2(ipart)=logGrad-fgrad;
    d2(ipart+fpNodes.ifirst)=sqrt(2*fpNodes.mass
      /((dot(logGrad,logGrad)+1e-15)*fpNodes.tau));
  }
}

void FreeParticleNodes::MatrixUpdate::acceptLastMove(int nslice) {
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
