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
#include "H5ReportBuilder.h"
#include "ScalarEstimator.h"
#include "EstimatorManager.h"
#include "ArrayBlockedEstimator.h"

H5ReportBuilder::H5ReportBuilder(const std::string& filename,
  const EstimatorManager::SimInfoWriter *simInfoWriter) 
  : filename(filename), simInfoWriter(simInfoWriter),
    fileID(0), writingGroupID(0), stepAttrID(0), dataset(0) {
  fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  simInfoWriter->writeH5(fileID);
}

H5ReportBuilder::~H5ReportBuilder() {
  H5Fclose(fileID);
}

void H5ReportBuilder::startWritingGroup(EstimatorManager& manager) {
  nstep=manager.nstep;
  istep=0;
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  writingGroupID = H5Gcreate2(fileID,"estimators",H5P_DEFAULT,
                             H5P_DEFAULT,H5P_DEFAULT);
#else
  writingGroupID = H5Gcreate(fileID,"estimators",0);
#endif
  hsize_t dims = 1;
  hid_t dataspaceID = H5Screate_simple(1, &dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  stepAttrID = H5Acreate2(writingGroupID, "nstep", H5T_NATIVE_INT, 
                         dataspaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
  stepAttrID = H5Acreate(writingGroupID, "nstep", H5T_NATIVE_INT, 
                         dataspaceID, H5P_DEFAULT);
#endif
  H5Sclose(dataspaceID);
  H5Awrite(stepAttrID, H5T_NATIVE_INT, &istep);
  dataset.resize(0);
  for (EstimatorManager::EstimatorIter est=manager.estimator.begin();
       est!=manager.estimator.end(); ++est) {
    (*est)->startReport(*this);
  }
}

void H5ReportBuilder::startScalarReport(const ScalarEstimator& est) {
  hsize_t dims=nstep;
  hid_t dataSpaceID = H5Screate_simple(1, &dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dataSetID = H5Dcreate2(writingGroupID, est.getName().c_str(),
                               H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT, 
                               H5P_DEFAULT, H5P_DEFAULT);
#else
  hid_t dataSetID = H5Dcreate(writingGroupID, est.getName().c_str(),
                              H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT);
#endif
  H5Sclose(dataSpaceID);
  const std::string& unitName(est.getUnitName());
  if (unitName!="") {
    hsize_t dims=1;
    hid_t dataSpaceID = H5Screate_simple(1, &dims, NULL);
    hid_t strType = H5Tcopy(H5T_C_S1);
    H5Tset_size(strType, unitName.length());
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t attrID = H5Acreate2(dataSetID, "unit", strType, dataSpaceID,
                              H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t attrID = H5Acreate(dataSetID, "unit", strType, dataSpaceID,
                             H5P_DEFAULT);
#endif
    H5Awrite(attrID, strType, unitName.c_str());
    H5Sclose(dataSpaceID);
    H5Aclose(attrID);
  }
  dataset.push_back(dataSetID);
}

void H5ReportBuilder::startAccRejReport(const AccRejEstimator& est) {

}

void H5ReportBuilder::reportAccRejStep(const AccRejEstimator& est) {
}

void H5ReportBuilder::writeStep(EstimatorManager& manager) {
  dset=dataset.begin();
  for (EstimatorManager::EstimatorIter est=manager.estimator.begin();
       est!=manager.estimator.end(); ++est) {
    (*est)->reportStep(*this);
  }
  H5Awrite(stepAttrID, H5T_NATIVE_INT, &(++istep));
  H5Fflush(fileID, H5F_SCOPE_LOCAL);
  if (istep == nstep) {
    H5Aclose(stepAttrID);
    for (DataSetIter d=dataset.begin(); d!=dataset.end(); ++d) {
      H5Dclose(*d);
    }
    H5Gclose(writingGroupID);
    dataset.resize(0);
  }
}

void H5ReportBuilder::reportScalarStep(const ScalarEstimator& est) {
  hid_t space = H5Dget_space(*dset); 
  H5Sselect_none(space);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR==6)&&(H5_VERS_RELEASE>=7))
  hsize_t element = istep;
  H5Sselect_elements(space, H5S_SELECT_SET, 1, &element);
#else
  hsize_t element[][1] = {{istep}};
  H5Sselect_elements(space, H5S_SELECT_SET, 1, (const hsize_t **)element);
#endif

  float a=est.getValue();
  hsize_t dims=1, maxdims=1;
  hid_t mspace = H5Screate_simple(1, &dims, &maxdims);
  H5Dwrite(*dset, H5T_NATIVE_FLOAT, mspace, space, H5P_DEFAULT, &a);
  dset++;
}
 
void H5ReportBuilder::reportArrayBlockedStep(const ArrayBlockedEstimator& est) {
  est.normalize();
  H5Dwrite(*dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           est.getData());
  dset++;
  if (est.hasError()) {
    H5Dwrite(*dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             est.getError());
    dset++;
  }
  est.unnormalize();
}

void H5ReportBuilder::startArrayBlockedReport(const ArrayBlockedEstimator& est) {
  hsize_t dims[est.getNDim()];
  for (int i=0; i<est.getNDim(); ++i) dims[i]=est.getExtent(i);

  hid_t dataSpaceID = H5Screate_simple(est.getNDim(), dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dataSetID = H5Dcreate2(writingGroupID, est.getName().c_str(),
                               H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT,
                               H5P_DEFAULT, H5P_DEFAULT);
#else
  hid_t dataSetID = H5Dcreate(writingGroupID, est.getName().c_str(),
                              H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT);
#endif

  if (est.hasScale()) {
    hsize_t dim=est.getNDim();
    hid_t attrSpaceID = H5Screate_simple(1, &dim, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t attrID = H5Acreate2(dataSetID, "scale", H5T_NATIVE_DOUBLE,
                              attrSpaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t attrID = H5Acreate(dataSetID, "scale",
                             H5T_NATIVE_DOUBLE, attrSpaceID, H5P_DEFAULT);
#endif
    H5Sclose(attrSpaceID);
    H5Awrite(attrID, H5T_NATIVE_DOUBLE, est.getScale()); 
    H5Aclose(attrID);
  }

  if (est.hasOrigin()) {
    hsize_t dim=est.getNDim();
    hid_t attrSpaceID = H5Screate_simple(1, &dim, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t attrID = H5Acreate2(dataSetID, "origin", H5T_NATIVE_DOUBLE,
                              attrSpaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t attrID = H5Acreate(dataSetID, "origin",
                             H5T_NATIVE_DOUBLE, attrSpaceID, H5P_DEFAULT);
#endif
    H5Sclose(attrSpaceID);
    H5Awrite(attrID, H5T_NATIVE_DOUBLE, est.getOrigin()); 
    H5Aclose(attrID);
  }

  H5Sclose(dataSpaceID);
  dataset.push_back(dataSetID);

  if (est.hasError()) {
    hid_t dataSpaceID = H5Screate_simple(est.getNDim(), dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dataSetID = H5Dcreate2(writingGroupID, (est.getName()+"_err").c_str(),
                                H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t dataSetID = H5Dcreate(writingGroupID, (est.getName()+"_err").c_str(),
                                H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT);
#endif
    H5Sclose(dataSpaceID);
    dataset.push_back(dataSetID);
  }
}
