#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "SpinSetter.h"
#include "base/Paths.h"
#include "base/Species.h"
#include "stats/MPIManager.h"
#include "util/RandomNumGenerator.h"
#include "util/SuperCell.h"
#include <cstdlib>
#include <blitz/tinyvec-et.h>

SpinSetter::SpinSetter(Paths& paths,MPIManager *mpi)
   : paths(paths), npart(paths.getNPart()), mpi(mpi) {
}

void SpinSetter::run() {
std::cout << "setting spin" << std::endl;
  int nslice=paths.getNSlice();
  int ifirstSlice=paths.getLowestOwnedSlice(false);
  int ilastSlice=paths.getHighestOwnedSlice(false);
  SuperCell cell=paths.getSuperCell();
  if (!mpi || mpi->isCloneMain()) {
    for (int ipart=0; ipart<npart; ++ipart) {
      (*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,ifirstSlice,1)))
      =SVec(-0.485332,-0.433568,0.476574,0.375401);
    }
  }
  // Copy first slice to workers.
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Bcast(paths.getAuxBead(0,ifirstSlice,1), npart*4,
                               MPI::DOUBLE,0);
  } 
#endif
  // Copy to other slices.
  for (int islice=ifirstSlice+1; islice<ilastSlice; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      (*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,islice,1)))
        =(*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,ifirstSlice,1)));
    } 
  }
  if (paths.isDouble()) {
  for (int islice=ifirstSlice; islice<ilastSlice; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      (*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,islice+nslice/2,1)))
        =(*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,islice,1)));
    } 
  }
  }
}
