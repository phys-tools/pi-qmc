// $Id$
/*  Copyright (C) 2007 John B. Shumway, Jr.

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
#include "RandomNumGenerator.h"
#include "stats/MPIManager.h"
#include "StructReader.h"
#include "Paths.h"
#include "Species.h"
#include <blitz/tinyvec-et.h>
#include "SuperCell.h"
#include <hdf5.h>

StructReader::StructReader(Paths& paths, const std::string& filename,
                           MPIManager *mpi)
   : paths(paths), filename(filename), mpi(mpi) {
}


void StructReader::run() {
  int nslice=paths.getNSlice();
  int ifirstSlice=paths.getLowestOwnedSlice(false);
  int ilastSlice=paths.getHighestOwnedSlice(false);
  int npart=0;
  // Read the particle positions from the struct.h5 file.
  if (!mpi || mpi->isCloneMain()) {
    std::cout << "Reading positions from file " << filename << std::endl;
    hid_t fileID = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dataSetID = H5Dopen2(fileID, "coords/coord", H5P_DEFAULT);
#else
    hid_t dataSetID = H5Dopen(fileID, "coords/coord");
#endif
    hid_t dataSpaceID = H5Dget_space(dataSetID);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataSpaceID, dims, NULL);
    H5Sclose(dataSpaceID);
    npart=dims[0];
    VArray coords(npart);
    H5Dread(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            coords.data());
    H5Dclose(dataSetID);
    H5Fclose(fileID);
    for (int ipart=0; ipart<npart; ++ipart) {
       paths(ipart,ifirstSlice)=coords(ipart);
    }
  }
  // Copy first slice to workers.
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Bcast(&paths(0,ifirstSlice), npart*NDIM,
                               MPI::DOUBLE,0);
  } 
#endif
  // Copy to other slices.
  for (int islice=ifirstSlice; islice<=ilastSlice; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      paths(ipart,islice)=paths(ipart,ifirstSlice);
    } 
  }
  if (paths.isDouble()) {
    for (int islice=ifirstSlice; islice<=ilastSlice; ++islice) {
      for (int ipart=0; ipart<npart; ++ipart) {
        paths(ipart,islice+nslice/2)=paths(ipart,islice);
      } 
    }
  } 
  paths.setBuffers();
}
