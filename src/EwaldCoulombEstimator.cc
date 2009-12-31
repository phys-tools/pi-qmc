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

// constructor trad ewald
EwaldCoulombEstimator::EwaldCoulombEstimator(
  const SimulationInfo& simInfo, const Action* action, const double epsilon,
  const double rcut, const double kcut,
  MPIManager *mpi, const std::string& unitName, double scale, double shift, 
  const double kappa, const int nImages, const bool testEwald)
  : ScalarEstimator("coulomb_energy",unitName,scale,shift),
    testEwald(testEwald), kappa(kappa), kcut(kcut),
    ewaldSum(*new TradEwaldSum(*simInfo.getSuperCell(), 
                                simInfo.getNPart(),rcut,kcut,kappa)),
    cell(*simInfo.getSuperCell()),
    energy(0), etot(0), enorm(0), vgrid(1001), nradial(1001),
    rcut(rcut), dr(rcut/1000), drinv(1./dr),
    action(action), epsilon(epsilon),q(simInfo.getNPart()),
    r(simInfo.getNPart()), nImages(nImages), sphereR(0.),
    mpi(mpi) {
  
  for (int i=0; i<q.size(); ++i) q(i)=simInfo.getPartSpecies(i).charge;
  ewaldSum.getQArray() = q;
  ewaldSum.evalSelfEnergy();
  vgrid(0) = -ewaldSum.evalFR0()/epsilon;
  for (int i=1; i<nradial; ++i) {
    vgrid(i) = -ewaldSum.evalFR(i*dr)/epsilon;
  }
  
  // set up for summing over images for trad ewald sum
  for (int i=0; i< NDIM; i++){
    sphereR = (cell[i]>sphereR)?cell[i]:sphereR;
  }
  sphereR *=nImages;
  boxImageVecs.resize(0);
  findBoxImageVectors(cell);
  std::cout << "In Ewald CEstimator Using Trad Ewald with sphereR: " << sphereR << std::endl;

}

// constructor opt ewald
EwaldCoulombEstimator::EwaldCoulombEstimator(const SimulationInfo& simInfo,
  const Action* action, const double epsilon, const double rcut, 
  const double kcut, MPIManager *mpi, const std::string& unitName, 
  double scale, double shift)
  : ScalarEstimator("coulomb_energy",unitName,scale,shift),
    testEwald(false),
    ewaldSum(*new OptEwaldSum(*simInfo.getSuperCell(),
                               simInfo.getNPart(),rcut,kcut,4*kcut,8)),
    cell(*simInfo.getSuperCell()),
    energy(0), etot(0), enorm(0), vgrid(1001), nradial(1001),
    rcut(rcut), dr(rcut/1000), drinv(1./dr),
    action(action), epsilon(epsilon),q(simInfo.getNPart()),
    r(simInfo.getNPart()), nImages(0), mpi(mpi) {
   
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

  if (nImages >1){
    for (int jpart=0; jpart<ipart; ++jpart) {
      for (unsigned int img=0; img<boxImageVecs.size(); img++){
	Vec boxImage;
	for (int l=0; l<NDIM; l++) boxImage[l]=boxImageVecs[img][l];// eventually change data strucure to use tinyvecs.
		
	Vec delta=end-paths(jpart,islice);
	cell.pbc(delta);
	double r=sqrt(dot(delta+boxImage,delta+boxImage));
	energy+=q(ipart)*q(jpart)*(1./(r*epsilon) -ewaldSum.evalFR(r)/epsilon);/// could use the vgrid later after testing
	  }
    }
  }else{
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
  }
 
  // Add long range contribution.
  if (ipart==0) {
    paths.getSlice(islice,r);
    energy += ewaldSum.evalLongRange(r)/epsilon;
   }
 
  //check if testing is required.
  if (mpi && mpi->isMain() && testEwald) testEwaldTotalCharge(paths);
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

void EwaldCoulombEstimator::testEwaldTotalCharge( const Paths& paths){
  //test kcut and nImages by calculating total charge in the box. 
  int islice=1; 
  int n=500;  
  double h=cell.a[0]/n;
  int totk=ewaldSum.gettotk(); 
  Vec dk=ewaldSum.getDeltak(); 
  std::complex<double> I(0,1);
  double intqr=0;
  std::complex<double> INTqofk=0;
  
  for (int jpart=0; jpart<paths.getNPart(); ++jpart) {
    for (unsigned int img=0; img<boxImageVecs.size(); img++){
      Vec boxImage;
      for (int l=0; l<NDIM; l++) boxImage[l]=boxImageVecs[img][l];
         // eventually change data strucure to use tinyvecs.

      Vec rj=paths(jpart,islice);
      
      Vec rim =rj+boxImage;
      double tmpx;
      tmpx=  exp(-kappa*kappa*(cell.a[0]/2-rim[0])*(cell.a[0]/2-rim[0]) );
      tmpx+=exp(-kappa*kappa*(-cell.a[0]/2-rim[0])*(-cell.a[0]/2-rim[0]) );
      for (int i=2;i<=n/2;i++)
        tmpx += 2* exp(-kappa*kappa*(h*(2*i-2)-rim[0]-cell.a[0]/2  )
                         *(h*(2*i-2)-rim[0]-cell.a[0]/2  ));
      for (int i=1;i<=n/2;i++)
        tmpx += 4*exp(-kappa*kappa*(h*(2*i-1)-rim[0]-cell.a[0]/2 )
                         *(h*(2*i-1)-rim[0]-cell.a[0]/2 ));
      double tmpy;
      tmpy=  exp(-kappa*kappa*(cell.a[1]/2-rim[1])*(cell.a[1]/2-rim[1]) );
      tmpy+=exp(-kappa*kappa*(-cell.a[1]/2-rim[1])*(-cell.a[1]/2-rim[1]) );
      for (int i=2;i<=n/2;i++)
        tmpy += 2*exp(-kappa*kappa*(h*(2*i-2)-rim[1]-cell.a[1]/2  )
                         *(h*(2*i-2)-rim[1]-cell.a[1]/2  ));
      for (int i=1;i<=n/2;i++)
        tmpy += 4*exp(-kappa*kappa*(h*(2*i-1)-rim[1]-cell.a[1]/2 )
                         *(h*(2*i-1)-rim[1]-cell.a[1]/2 ));
      double tmpz;
      tmpz=  exp(-kappa*kappa*(cell.a[2]/2-rim[2])*(cell.a[2]/2-rim[2]) );
      tmpz+= exp(-kappa*kappa*(-cell.a[2]/2-rim[2])*(-cell.a[2]/2-rim[2]) );
      for (int i=2;i<=n/2;i++)
        tmpz += 2*exp(-kappa*kappa*(h*(2*i-2)-rim[2]-cell.a[2]/2  )
                         *(h*(2*i-2)-rim[2]-cell.a[2]/2  ));
      for (int i=1;i<=n/2;i++)
        tmpz += 4*exp(-kappa*kappa*(h*(2*i-1)-rim[2]-cell.a[2]/2 )
                         *(h*(2*i-1)-rim[2]-cell.a[2]/2 ));
      intqr+=q(jpart)*pow(kappa,3)/pow(3.141592653,1.5)
                     *tmpx*tmpy*tmpz*h*h*h/27.;
    }
  }
  std::cout << " Total Charge from space integral over space " 
            << "charge density ::  "<< intqr << std::endl;

#ifdef _OPENMP
#pragma omp parallel            // reduction(+:INTqofk)
#endif
  {
    for (int jpart=0; jpart<paths.getNPart(); ++jpart) {
#ifdef _OPENMP
#pragma omp parallel for private(l,boxImage,kk)  reduction(+:INTqofk)
#endif
      for (unsigned int img=0; img<boxImageVecs.size(); img++){
	Vec boxImage;
	for (int l=0; l<NDIM; l++) boxImage[l]=boxImageVecs[img][l];
            // eventually change data strucure to use tinyvecs.
	
	Vec rj=paths(jpart,islice);
	
	Vec rim =rj+boxImage;

	for (int kk=0; kk<totk; kk++){
	  double h=cell.a[0]/n; 
	  Vec k=ewaldSum.getkvec(kk)*dk;
	  std::complex<double>tmpx;
	  tmpx=  exp(-I*k[0]*cell.a[0]*0.5)
            *exp(-kappa*kappa*(cell.a[0]*0.5-rim[0])*(cell.a[0]*0.5-rim[0]) );
	  tmpx+=exp(I*k[0]*cell.a[0]*0.5)
            *exp(-kappa*kappa*(-cell.a[0]*0.5-rim[0])*(-cell.a[0]*0.5-rim[0]) );
	  for (int i=2;i<=n*0.5;i++)
            tmpx+=2.0* exp( -I*k[0]*( h*(2*i-2)-cell.a[0]*0.5 ))  
               *exp(-kappa*kappa*(h*(2*i-2)-rim[0]-cell.a[0]*0.5  )
                      *(h*(2*i-2)-rim[0]-cell.a[0]*0.5  ));
	  for (int i=1;i<=n*0.5;i++)
            tmpx+=4.0*exp( -I*k[0]*( h*(2*i-1)-cell.a[0]*0.5))
               *exp(-kappa*kappa*(h*(2*i-1)-rim[0]-cell.a[0]*0.5 )
                      *(h*(2*i-1)-rim[0]-cell.a[0]*0.5 ));
	  std::complex<double> tmpy;
	  tmpy=  exp(-I*k[1]*cell.a[1]*0.5)
            *exp(-kappa*kappa*(cell.a[1]*0.5-rim[1])*(cell.a[1]*0.5-rim[1]) );
	  tmpy+=exp(I*k[1]*cell.a[1]*0.5)
            *exp(-kappa*kappa*(-cell.a[1]*0.5-rim[1])*(-cell.a[1]*0.5-rim[1]) );
	  for (int i=2;i<=n*0.5;i++)
            tmpy+=2.0*exp( -I*k[1]*( h*(2*i-2)-cell.a[1]*0.5 ))  
               *exp(-kappa*kappa*(h*(2*i-2)-rim[1]-cell.a[1]*0.5  )
                       *(h*(2*i-2)-rim[1]-cell.a[1]*0.5  ));
	  for (int i=1;i<=n*0.5;i++)  
            tmpy+=4.0*exp( -I*k[1]*( h*(2*i-1)-cell.a[1]*0.5))
                     *exp(-kappa*kappa*(h*(2*i-1)-rim[1]-cell.a[1]*0.5 )
                             *(h*(2*i-1)-rim[1]-cell.a[1]*0.5 ));
	  std::complex<double> tmpz;
	  tmpz=  exp(-I*k[2]*cell.a[2]*0.5)
            *exp(-kappa*kappa*(cell.a[2]*0.5-rim[2])*(cell.a[2]*0.5-rim[2]) );
	  tmpz+=exp(I*k[2]*cell.a[2]*0.5)
            *exp(-kappa*kappa*(-cell.a[2]*0.5-rim[2])*(-cell.a[2]*0.5-rim[2]) );
	  for (int i=2;i<=n*0.5;i++)  
             tmpz+=2.0* exp( -I*k[2]*( h*(2*i-2)-cell.a[2]*0.5 ))  
               *exp(-kappa*kappa*(h*(2*i-2)-rim[2]-cell.a[2]*0.5  )
                       *(h*(2*i-2)-rim[2]-cell.a[2]*0.5  ));
	  for (int i=1;i<=n*0.5;i++)  
             tmpz+=4.0* exp( -I*k[2]*( h*(2*i-1)-cell.a[2]*0.5))  
               *exp(-kappa*kappa*(h*(2*i-1)-rim[2]-cell.a[2]*0.5 )
                       *(h*(2*i-1)-rim[2]-cell.a[2]*0.5 ));
	  INTqofk+=q(jpart)*pow(kappa,3)/pow(3.141592653,1.5)
                                        /27.0*tmpx*tmpy*tmpz*h*h*h;
	}
      }
    } //end of omp parallel section
    std::cout << " Total Charge from k-space integral over k-space" 
              << " charge density ::  " << INTqofk <<  std::endl;
  }
  
}

void EwaldCoulombEstimator::findBoxImageVectors(const SuperCell &a) { 
 // 3D case first... to be extended to 2D and 1D soon...after testing 3D case
  std :: vector<std :: vector<double> > vertices(8); 
  for (unsigned int i=0; i< vertices.size(); i++){
    vertices[i].resize(NDIM);
  }
  vertices[0][0]=-a[0]/2;  vertices[0][1]=-a[1]/2;  vertices[0][2]=-a[2]/2;
  vertices[1][0]=a[0]/2;  vertices[1][1]=-a[1]/2;  vertices[1][2]=-a[2]/2;
  vertices[2][0]=a[0]/2;  vertices[2][1]=a[1]/2;  vertices[2][2]=-a[2]/2;
  vertices[3][0]=-a[0]/2;  vertices[3][1]=a[1]/2;  vertices[3][2]=-a[2]/2;
  vertices[4][0]=-a[0]/2;  vertices[4][1]=-a[1]/2;  vertices[4][2]=a[2]/2;
  vertices[5][0]=a[0]/2;  vertices[5][1]=-a[1]/2;  vertices[5][2]=a[2]/2;
  vertices[6][0]=a[0]/2;  vertices[6][1]=a[1]/2;  vertices[6][2]=a[2]/2;
  vertices[7][0]=-a[0]/2;  vertices[7][1]=a[1]/2;  vertices[7][2]=a[2]/2;

  // Check that all the box images are inside the sphere of radius 
  // nImages*max(cell side) and save the location of these boxes 
  // (center of box).
  std :: vector<double> L(NDIM);
  for (int nx=-nImages; nx<=nImages; nx++){
    for (int ny=-nImages; ny<=nImages; ny++){
      for (int nz=-nImages; nz<=nImages; nz++){
	
	L[0]=nx*a[0]; L[1]=ny*a[1]; L[2]=nz*a[2];
	int flag=1;
	for (int v=0;v<8;v++){
	  double tmpR=0;
	  for (int i=0; i< NDIM;i++) {
	    tmpR+=(L[i]+vertices[v][i])*(L[i]+vertices[v][i]);
	  }
	  if (sqrt(tmpR) > sphereR) {
	    flag=0;
	    break;
	  }
	}
	if (flag==1) {
	  boxImageVecs.push_back(L);
	}
	
      }
    }
  }  
}
