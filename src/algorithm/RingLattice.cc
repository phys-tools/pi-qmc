#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "RingLattice.h"
#include "base/Paths.h"
#include "base/Species.h"
#include "stats/MPIManager.h"
#include <blitz/tinyvec-et.h>
#include "util/SuperCell.h"

RingLattice::RingLattice(Paths& paths, double radius, double angle0, 
  double anglef, double anglex, MPIManager *mpi) 
  : radius(radius), anglex(anglex), angle0(angle0), anglef(anglef),
    ifirst(0),npart(paths.getNPart()), paths(paths), mpi(mpi) {
}

RingLattice::RingLattice(Paths& paths, double radius, double angle0, 
  double anglef, double anglex, const Species& species, MPIManager *mpi)
  : radius(radius), anglex(anglex), angle0(angle0), anglef(anglef), 
    ifirst(species.ifirst), npart(species.count), paths(paths), mpi(mpi) {
}

void RingLattice::run() {
#if NDIM==2 || NDIM==3
  int nslice = paths.getNSlice();
  // Place the particles on the ring.
  int ifirstSlice = paths.getLowestOwnedSlice(false);
  int ilastSlice = paths.getHighestOwnedSlice(false);
  SuperCell cell=paths.getSuperCell();
  if (!mpi || mpi->isCloneMain()) {
    int ipart = ifirst;
    double angleDist = (anglef - angle0) / (npart-1);
    for (int i=0; i<npart; ++i) {
      paths(ipart, ifirstSlice) =
#if NDIM==2 
                                  Vec(radius * cos(angle0 + angleDist * i + anglex),
                                      radius * sin(angle0 + angleDist * i + anglex));
#endif
#if NDIM==3
                                  Vec(radius * cos(angle0 + angleDist * i + anglex),
                                      radius * sin(angle0 + angleDist * i + anglex),
                                      0);
#endif
      cell.pbc(paths(ipart, ifirstSlice));
      ++ipart;
    }
  }
  // Copy first slice to workers.
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Bcast(&paths(ifirst,ifirstSlice), npart*NDIM,
                               MPI::DOUBLE,0);
  }
#endif
  // Copy to other slices.
  for (int islice=ifirstSlice+1; islice<=ilastSlice; ++islice) {
    for (int ipart=ifirst; ipart<ifirst+npart; ++ipart) {
      paths(ipart,islice)=paths(ipart,ifirstSlice);
    }
  }
  if (paths.isDouble()) {
    for (int islice=ifirstSlice; islice<=ilastSlice; ++islice) {
      for (int ipart=ifirst; ipart<ifirst+npart; ++ipart) {
        paths(ipart,islice+nslice/2)=paths(ipart,islice);
      }
    }
  }
  paths.setBuffers();
#endif
}
