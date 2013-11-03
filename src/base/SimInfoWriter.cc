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

#include "SimInfoWriter.h"
#include "SimulationInfo.h"
#include "util/SuperCell.h"
#include <string>
#include <iostream>
#include <vector>

SimInfoWriter::SimInfoWriter(const SimulationInfo &simInfo)
  : simInfo(simInfo) {
}

void SimInfoWriter::writeH5(hid_t fileID) const {
  std::cout << "Writing simulation info to pimc.h5." << std::endl;
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t groupID = H5Gcreate2(fileID, "simInfo", H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
#else
  hid_t groupID = H5Gcreate(fileID, "simInfo", 0);
#endif
  hsize_t dims = 1;
  hid_t dataSpaceID = H5Screate_simple(1, &dims, NULL);

#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dataSetID = H5Dcreate2(groupID, "temperature",
                         H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
#else
  hid_t dataSetID = H5Dcreate(groupID, "temperature",
                              H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT);
#endif
  H5Dwrite(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &simInfo.temperature);
  H5Dclose(dataSetID);

#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dataSetID = H5Dcreate2(groupID, "deltaTau",
                         H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
#else
  dataSetID = H5Dcreate(groupID, "deltaTau",
                              H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT);
#endif
  H5Dwrite(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &simInfo.tau);
  H5Dclose(dataSetID);

#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dataSetID = H5Dcreate2(groupID, "nslice",
                         H5T_NATIVE_INT, dataSpaceID, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
#else
  dataSetID = H5Dcreate(groupID, "nslice",
                              H5T_NATIVE_INT, dataSpaceID, H5P_DEFAULT);
#endif
  H5Dwrite(dataSetID, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &simInfo.nslice);
  H5Dclose(dataSetID);

#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dataSetID = H5Dcreate2(groupID, "npart",
                         H5T_NATIVE_INT, dataSpaceID, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
#else
  dataSetID = H5Dcreate(groupID, "npart",
                              H5T_NATIVE_INT, dataSpaceID, H5P_DEFAULT);
#endif
  H5Dwrite(dataSetID, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &simInfo.npart);
  H5Dclose(dataSetID);
  H5Sclose(dataSpaceID);

  dims = NDIM;
  dataSpaceID = H5Screate_simple(1, &dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dataSetID = H5Dcreate2(groupID, "superCell",
                         H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
#else
  dataSetID = H5Dcreate(groupID, "superCell",
                              H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT);
#endif
  H5Dwrite(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &simInfo.superCell->a);
  H5Dclose(dataSetID);
  H5Sclose(dataSpaceID);
  H5Gclose(groupID);

  H5Fflush(fileID, H5F_SCOPE_LOCAL);
}
