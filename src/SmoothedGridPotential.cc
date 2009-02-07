// $Id$
/*  Copyright (C) 2004-2008 John B. Shumway, Jr.

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
#include "SmoothedGridPotential.h"
#include "Beads.h"

SmoothedGridPotential::SmoothedGridPotential(const SimulationInfo& simInfo, const int maxLevel,
                                             const std::string& filename)
  : tau(simInfo.getTau()), npart(simInfo.getNPart()), nlevel(maxLevel), vindex(npart) {
  std::cout << "Smoothed Grid Potential nlevel = " << nlevel << std::endl;
  // Read the bandoffsets from grid.h5.
  hid_t fileID = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t groupID = H5Gopen2(fileID, (filename=="emagrids.h5"?"boffset":"/"),
                           H5P_DEFAULT);
#else
  hid_t groupID = H5Gopen(fileID, (filename=="emagrids.h5"?"boffset":"/"));
#endif
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dataSetID = H5Dopen2(groupID, "vh", H5P_DEFAULT);
#else
  hid_t dataSetID = H5Dopen(groupID, "vh");
#endif
  hid_t dataSpaceID = H5Dget_space(dataSetID);
  hsize_t dims[3];
  H5Sget_simple_extent_dims(dataSpaceID, dims, NULL);
  H5Sclose(dataSpaceID);
  n(0)=dims[0]; n(1)=dims[1]; n(2)=dims[2];
  vhtemp.resize(n);
  H5Dread(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          vhtemp.data()); 
  vetemp.resize(n);
  H5Dclose(dataSetID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dataSetID = H5Dopen2(groupID, "ve", H5P_DEFAULT );
#else
  dataSetID = H5Dopen(groupID, "ve");
#endif
  H5Dread(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          vetemp.data()); 
  H5Dclose(dataSetID);
  // Read the grid spacing.
  double a=0;
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dataSetID = H5Dopen2(groupID, "a", H5P_DEFAULT);
#else
  dataSetID = H5Dopen(groupID, "a");
#endif
  H5Dread(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&a);
  H5Dclose(dataSetID);
  b=1.0/a;
  H5Gclose(groupID);
  H5Fclose(fileID);
  // Convert grids from eV to Ha.
  if(filename=="boffset.h5") {
    vhtemp*=EVTOHA;
    vetemp*=EVTOHA;
  }
  // Flip sign of hole grid if read from emagrid.h5
  if(filename=="emagrids.h5")
    vhtemp*=-1;

  CArray3 kvegrid, kvhgrid, vesmooth, vhsmooth;
  kvegrid.resize(n);
  kvhgrid.resize(n);
  vesmooth.resize(n);
  vhsmooth.resize(n);
  vegrid.resize(nlevel+1);
  vhgrid.resize(nlevel+1);
  //what do I do if only one e or h?
  /*Vec e_m, h_m;
  for(int i=0; i<simInfo.getNPart(); i++) {
    std::string name=simInfo.getPartSpecies(i).name;
    if(name.substr(0,1)=="h"){
      h_m=((simInfo.getPartSpecies(i).anMass!=0)?(*simInfo.getPartSpecies(i).anMass)
         :(Vec(simInfo.getPartSpecies(i).mass,simInfo.getPartSpecies(i).mass,simInfo.getPartSpecies(i).mass)));
    }
    else if(name.substr(0,1)=="e"){
      e_m=((simInfo.getPartSpecies(i).anMass!=0)?(*simInfo.getPartSpecies(i).anMass)
         :(Vec(simInfo.getPartSpecies(i).mass,simInfo.getPartSpecies(i).mass,simInfo.getPartSpecies(i).mass)));
    }
  }*/
  Vec e_m(0.067, 0.067, 0.067), h_m(0.08, 0.08, 0.45); //need to read these in
  std::cout<<"e_m = ["<<e_m(0)<<", "<<e_m(1)<<", "<<e_m(2)<<"]"<<std::endl;
  std::cout<<"h_m = ["<<h_m(0)<<", "<<h_m(1)<<", "<<h_m(2)<<"]"<<std::endl;
  //FFT
  fftw_plan phin = fftw_plan_dft(3, n.data(), (fftw_complex*)kvhgrid.data(), (fftw_complex*)kvhgrid.data(), FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan pein = fftw_plan_dft(3, n.data(), (fftw_complex*)kvegrid.data(), (fftw_complex*)kvegrid.data(), FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan phout = fftw_plan_dft(3, n.data(), (fftw_complex*)vhsmooth.data(), (fftw_complex*)vhsmooth.data(), FFTW_BACKWARD, FFTW_MEASURE);
  fftw_plan peout = fftw_plan_dft(3, n.data(), (fftw_complex*)vesmooth.data(), (fftw_complex*)vesmooth.data(), FFTW_BACKWARD, FFTW_MEASURE);
  for(int i=0; i<n(0); i++){
    for(int j=0; j<n(1); j++){
      for(int k=0; k<n(2); k++){
        kvhgrid(i,j,k)=std::complex<double>(vhtemp(i,j,k),0.);
        kvegrid(i,j,k)=std::complex<double>(vetemp(i,j,k),0.);
      }
    }
  }
  fftw_execute(phin);
  fftw_execute(pein);
  std::string name;

  //loop over levels
  std::cout<<"Smoothing level " << std::flush;
  for(int ilevel=0; ilevel<=nlevel; ++ilevel){
    std::cout << ilevel << " " << std::flush;
    vhgrid(ilevel).resize(n);
    vegrid(ilevel).resize(n);
    vhsmooth=0., vesmooth=0.;
    //smooth
    Vec h_alpha=0., e_alpha=0.;
    Vec h_norm=0., e_norm=0.;
    double sigma=pow(2,ilevel)*tau/24.;
    for(int i=0; i<3; i++){
      h_alpha(i)=sigma/h_m(i);
      e_alpha(i)=sigma/e_m(i);
    }
    double TWObPI=b*TWOPI;
    //initialize periodic guassians and normalize
    for(int i=0; i<3; i++){
      pgh(i)=new PeriodicGaussian(h_alpha(i),TWObPI,n(i));
      pge(i)=new PeriodicGaussian(e_alpha(i),TWObPI,n(i));
      h_norm(i)=1./(*pgh(i))(0.);
      e_norm(i)=1./(*pge(i))(0.);
    }
    double xhfactor=0., yhfactor=0., xefactor=0., yefactor=0.;
    double ki, kj, kk;
    for(int i=0; i<n(0); i++){
      ki=((i<=n(0)/2)?i*TWObPI/n(0):(n(0)-i)*TWObPI/n(0));
      xhfactor=(*pgh(0))(ki);
      xefactor=(*pge(0))(ki);
      for(int j=0; j<n(1); j++){
        kj=((j<=n(1)/2)?j*TWObPI/n(1):(n(1)-j)*TWObPI/n(1));
        yhfactor=(*pgh(1))(kj);
        yefactor=(*pge(1))(kj);
        for(int k=0; k<n(2); k++){
          kk=((k<=n(2)/2)?k*TWObPI/n(2):(n(2)-k)*TWObPI/n(2));
          vhsmooth(i,j,k)=kvhgrid(i,j,k)*blitz::product(h_norm)*xhfactor*yhfactor*(*pgh(2))(kk);
          vesmooth(i,j,k)=kvegrid(i,j,k)*blitz::product(e_norm)*xefactor*yefactor*(*pge(2))(kk);
        }
      }
    }
    //iFFT
    fftw_execute(peout);
    fftw_execute(phout);
    for(int i=0; i<n(0); i++){
      for(int j=0; j<n(1); j++){
        for(int k=0; k<n(2); k++){
          vegrid(ilevel)(i,j,k)=vesmooth(i,j,k).real();
          vhgrid(ilevel)(i,j,k)=vhsmooth(i,j,k).real();
        }
      }
    }
    double norm=1./(double)blitz::product(n);
    vhgrid(ilevel)*=norm;
    vegrid(ilevel)*=norm;
    /*std::cout << "\tmin (eV)\tmean (eV)\tmax (eV)" <<std::endl;
    std::cout << "h\t" << HATOEV*min(vhgrid(ilevel))
              << "\t" << HATOEV*mean(vhgrid(ilevel))
              << "\t" << HATOEV*max(vhgrid(ilevel)) << std::endl;
    std::cout << "e\t" << HATOEV*min(vegrid(ilevel))
              << "\t" << HATOEV*mean(vegrid(ilevel))
              << "\t" << HATOEV*max(vegrid(ilevel)) << std::endl;*/

  }
  std::cout << std::endl << std::endl;
  //clean up
  fftw_destroy_plan(phin);
  fftw_destroy_plan(pein);
  fftw_destroy_plan(phout);
  fftw_destroy_plan(peout);
  fftw_cleanup();

  // Setup the vindex.
  for(int i=0; i<simInfo.getNPart(); i++) {
    std::string name=simInfo.getPartSpecies(i).name;
    if(name.substr(0,1)=="h")
      vindex(i)=0;
    else if(name.substr(0,1)=="e")
      vindex(i)=1;
  }
  //writeout
  /*{
    std::cout << "Writing the smoothed band offsets to emagrids.smoothed.h5" << std::endl;
    system("cp emagrids.h5 emagrids.smoothed.h5");
    std::string name;
    H5::H5File file("emagrids.smoothed.h5", H5F_ACC_RDWR);
    H5::Group Group(file.createGroup("boffset/smoothed"));
    hsize_t dims[] = {n(0),n(1),n(2)};
    hsize_t dims1[] = {1};
    hsize_t dims3[] = {3};
    H5::DataSpace dataSpace1(1,dims1);
    H5::DataSpace dataSpace3(1,dims3);
    H5::DataSpace dataSpace(3,dims);
    double a=1./b;
    H5::DataSet(Group.createDataSet("a", H5::PredType::NATIVE_DOUBLE, dataSpace1)).write(&a,H5::PredType::NATIVE_DOUBLE);
    H5::DataSet(Group.createDataSet("tau", H5::PredType::NATIVE_DOUBLE, dataSpace1)).write(&tau,H5::PredType::NATIVE_DOUBLE);
    H5::DataSet(Group.createDataSet("nlevel", H5::PredType::NATIVE_INT, dataSpace1)).write(&nlevel,H5::PredType::NATIVE_INT);
    H5::DataSet(Group.createDataSet("m_e", H5::PredType::NATIVE_DOUBLE, dataSpace3)).write(e_m.data(),H5::PredType::NATIVE_DOUBLE);
    H5::DataSet(Group.createDataSet("m_h", H5::PredType::NATIVE_DOUBLE, dataSpace3)).write(h_m.data(),H5::PredType::NATIVE_DOUBLE);
    std::cout << "Writing level " << std::flush;
    //loop over levels
    for(int ilevel=0; ilevel<=nlevel; ++ilevel){
      std::cout << ilevel << " " << std::flush;
      std::ostringstream temp;
      temp << ilevel;
      name=temp.str();
      H5::Group GroupS(Group.createGroup(name.c_str()));
      H5::DataSet(GroupS.createDataSet("vh", H5::PredType::NATIVE_DOUBLE, dataSpace)).write(vhgrid(ilevel).data(),
                                      H5::PredType::NATIVE_DOUBLE);
      H5::DataSet(GroupS.createDataSet("ve", H5::PredType::NATIVE_DOUBLE, dataSpace)).write(vegrid(ilevel).data(),
                                      H5::PredType::NATIVE_DOUBLE);
      H5::DataSet(GroupS.createDataSet("level", H5::PredType::NATIVE_INT, dataSpace1)).write(&ilevel,
                                      H5::PredType::NATIVE_INT);
    }
    std::cout << std::endl;
  }*/
}

SmoothedGridPotential::~SmoothedGridPotential(){
  for(int i=0; i<3; i++){
    delete pgh(i);
    delete pge(i);
  }
}

double SmoothedGridPotential::getActionDifference(
    const MultiLevelSampler& sampler, const int ilevel) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,ilevel);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex();
  const int nMoving=index.size();
  double deltaAction=0;
  // Read values off of the grid.
  for (int islice=0; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Add action for moving beads. (Evaluate v at midpoint)
      Vec r=movingBeads(iMoving,islice);
      cell.pbc(r);
      Vec delta=movingBeads(iMoving,islice+nStride);
      delta-=r; cell.pbc(delta); delta*=0.5; r+=delta;
      cell.pbc(r);
      deltaAction+=v(r,i,ilevel)*tau*nStride;
      // Subtract action for old beads.
      r=sectionBeads(i,islice);
      cell.pbc(r);
      delta=sectionBeads(i,islice+nStride);
      delta-=r; cell.pbc(delta); delta*=0.5; r+=delta;
      cell.pbc(r);
      deltaAction-=v(r,i,ilevel)*tau*nStride;
    }
  }
  return deltaAction;
}

double SmoothedGridPotential::v(Vec r, const int i, int ilevel) const {
  r*=b;
  int i1=(int)floor(r[0]), i2=(int)floor(r[1]), i3=(int)floor(r[2]);
  double x=r[0]-i1, y=r[1]-i2, z=r[2]-i3;
  i1+=n(0)/2;
  i2+=n(1)/2;
  i3+=n(2)/2;
  if (i1<0) {i1=0; x=0;};
  if (i2<0) {i2=0; y=0;};
  if (i3<0) {i3=0; z=0;};
  if (i1>n(0)-2) {i1=n(0)-2; x=1;};
  if (i2>n(1)-2) {i2=n(1)-2; y=1;};
  if (i3>n(2)-2) {i3=n(2)-2; z=1;};
  //if need a level beyond the levels smoothed, use the last level smoothed
  if(ilevel>nlevel) ilevel=nlevel;
  double V=0.;
  const Array3& U( (vindex(i)==0) ? vhgrid(ilevel) : vegrid(ilevel) );
  V =(1-z)*( (1-y)*( (1-x)*U(i1,i2,i3)    +x*U(i1+1,i2,i3) )
              + y *( (1-x)*U(i1,i2+1,i3)  +x*U(i1+1,i2+1,i3) ) )
      + z *( (1-y)*( (1-x)*U(i1,i2,i3+1)  +x*U(i1+1,i2,i3+1) )
              + y *( (1-x)*U(i1,i2+1,i3+1)+x*U(i1+1,i2+1,i3+1)));
  return V;
}

double SmoothedGridPotential::getTotalAction(const Paths& paths, int ilevel) const {
  return 0;
}


void SmoothedGridPotential::getBeadAction(const Paths& paths, int ipart,
    int islice, double& u, double& utau, double& ulambda, 
    Vec& fm, Vec& fp) const {
  Vec r=paths(ipart,islice);
  r*=b;
  int i1=(int)floor(r[0]), i2=(int)floor(r[1]), i3=(int)floor(r[2]);
  double x=r[0]-i1, y=r[1]-i2, z=r[2]-i3;
  i1+=n(0)/2; i2+=n(1)/2; i3+=n(2)/2;
  if (i1<0) {i1=0; x=0;} 
  if (i2<0) {i2=0; y=0;} 
  if (i3<0) {i3=0; z=0;}
  if (i1>n(0)-2) {i1=n(0)-2; x=1;}
  if (i2>n(1)-2) {i2=n(1)-2; y=1;}
  if (i3>n(2)-2) {i3=n(2)-2; z=1;}
  const Array3& U( (vindex(ipart)==0) ? vhgrid(0) : vegrid(0) );
  utau =(1-z)*( (1-y)*( (1-x)*U(i1,i2,i3)    +x*U(i1+1,i2,i3) )
                 + y *( (1-x)*U(i1,i2+1,i3)  +x*U(i1+1,i2+1,i3) ) )
         + z *( (1-y)*( (1-x)*U(i1,i2,i3+1)  +x*U(i1+1,i2,i3+1) )
                 + y *( (1-x)*U(i1,i2+1,i3+1)+x*U(i1+1,i2+1,i3+1)));
  u = utau*tau;
  ulambda=0.;
  fm[0] =(1-z)*( (1-y)*( (1-x)*U(i1,i2,i3)    +x*U(i1+1,i2,i3) )
                  + y *( (1-x)*U(i1,i2+1,i3)  +x*U(i1+1,i2+1,i3) ) )
          + z *( (1-y)*( (1-x)*U(i1,i2,i3+1)  +x*U(i1+1,i2,i3+1) )
                  + y *( (1-x)*U(i1,i2+1,i3+1)+x*U(i1+1,i2+1,i3+1)));
  fm[1] =(1-z)*( (1-y)*( (1-x)*U(i1,i2,i3)    +x*U(i1+1,i2,i3) )
                  + y *( (1-x)*U(i1,i2+1,i3)  +x*U(i1+1,i2+1,i3) ) )
          + z *( (1-y)*( (1-x)*U(i1,i2,i3+1)  +x*U(i1+1,i2,i3+1) )
                  + y *( (1-x)*U(i1,i2+1,i3+1)+x*U(i1+1,i2+1,i3+1)));
  fm[2] =(1-z)*( (1-y)*( (1-x)*U(i1,i2,i3)    +x*U(i1+1,i2,i3) )
                  + y *( (1-x)*U(i1,i2+1,i3)  +x*U(i1+1,i2+1,i3) ) )
          + z *( (1-y)*( (1-x)*U(i1,i2,i3+1)  +x*U(i1+1,i2,i3+1) )
                  + y *( (1-x)*U(i1,i2+1,i3+1)+x*U(i1+1,i2+1,i3+1)));
  fm*=0.5;
  fp=fm;
}
