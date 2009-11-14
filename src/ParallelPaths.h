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
#ifndef __ParallelPaths_h_
#define __ParallelPaths_h_
template <int TDIM> class Beads;
class Permutation;
class MPIManager;
class BeadFactory;
class SuperCell;
#include "Paths.h"
#include <blitz/array.h>

/// Storage for paths in a parallel simulation.
/// This is parallelization in imaginary time.
/// The parallelization strategy is like parallelzing a lattice model in 1D.
/// The simplest idea is to break time up into nworker domains, each 
/// each approximately equal in size.
/// The size of the domain stored on the worker is nProcSlice.
/// For example, if we have 1024 slices and four workers, 
/// nslice=256.
/// The boundaries between slices should be frozen. 
/// We accomplish this by freezing slice i=0 and slice i=nProcSlice.
/// Note that slice i=0 is slice i=nprocSlice on the previous worker,
/// and slice i=nprocSlice is slice i=0 on the next worker.
/// For single-slice or multi-level sampling, each worker 
/// should only move beads on slices i=1 through i=nprocSlice-1. 
/// However, for measurements and displace moves, it is important that
/// all slices, including the buffers, be owned by exactly one worker.
/// That is other workers may have copies of buffer slice data, but only
/// the owner (owning worker) should evalauate estimators or perform
/// displacement movements on the beads.
/// We adopt the convention that each worker owns i=nprocSlice.
/// To compute estimators or make displace moves, it is necessary that
/// each worker also have a read-only copy of slice i=nprocSlice+1.
///
/// To summarize, each worker has nproc+2 slices stored on it.
/// Slice i=0 is read only and is owned by the previous worker.
/// Slice i=nprocSlice is owned by the worker, but is only moved in
/// displace moves.
/// Slice i=nprocSlice+1 is read only and is owned by the next worker.
/// 
/// To communicate these ideas with the rest of the code, which may
/// deal with serial paths, parallel paths, or double parallel paths,
/// we have three informative methods: getLowestOwnedSlice,
/// getHighestOwnedSlice, and getHighestSampledSlice.
/// These should return 1, nprocSlice, and nprocSlice-nsampling.
/// 
/// Finally, we have shiftWorkers and setBuffers method to communicate
/// slices between workers.
/// @version $Revision$
/// @author John Shumway
class ParallelPaths : public Paths {
public:
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  ParallelPaths(int npart, int nslice, double tau, const SuperCell& cell,
                MPIManager &mpi, const BeadFactory&);
  /// Destructor.
  virtual ~ParallelPaths();
  /// Loop over links, calling a LinkSummable object.
  virtual void sumOverLinks(LinkSummable&) const;
  /// Get a reference to a bead.
  virtual Vec& operator()(int ipart, int islice);
  /// Get a const reference to a bead.
  virtual const Vec& operator()(int ipart, int islice) const;
  /// Get a reference to a bead by offset.
  virtual Vec& operator()(int ipart, int islice, int istep);
  /// Get a const reference to a bead by offset.
  virtual const Vec&
    operator()(int ipart, int islice, int istep) const;
  /// Get a relative displacement a bead by offset.
  virtual Vec delta(int ipart, int islice, int istep) const;
  /// Get beads.
  virtual void getBeads(int ifirstSlice, Beads<NDIM>& ) const;
  /// Get a slice.
  virtual void getSlice(int islice, VArray& ) const;
  /// Get the number of slices on this processor.
  virtual int getNProcSlice(){return nprocSlice;}
  /// Get the number of slices unique to this processor.
  virtual int getNUniqueSlice(){return nprocSlice-1;}
  /// Get auxialiary bead.
  virtual const void* getAuxBead(int ipart, int islice, int iaux) const;
  /// Get auxialiary bead.
  virtual void* getAuxBead(int ipart, int islice, int iaux);
  /// Put beads.
  virtual void putBeads(int ifirstSlice,
                        const Beads<NDIM>&, const Permutation&) const;
  virtual void putDoubleBeads(
                 int ifirstSlice1,Beads<NDIM>&, Permutation&,
                 int ifirstSlice2,Beads<NDIM>&, Permutation&) const;
  /// Get the global permuation.
  /// @bug Placeholder, not correct for MPI.
    virtual const Permutation& getGlobalPermutation() const;
    virtual   const Permutation& getPermutation() const {return permutation;}
  virtual int getLowestOwnedSlice(bool d) const {return ifirst+1;}
  virtual int getHighestOwnedSlice(bool d) const {return ifirst+nprocSlice;}
  virtual int getHighestSampledSlice(int n, bool d) const {
    return ifirst+nprocSlice-n;}
  virtual bool isOwnedSlice(int islice) const {
    return islice > ifirst && islice <= ifirst+nprocSlice;}
  virtual void shift(int ishift);
  virtual void setBuffers();
  virtual bool is() const {return true;}
  virtual void clearPermutation();
private:
  /// Worker nubmer.
  const int iworker;
  /// Number of workers.
  const int nworker;
  /// The first slice on this processor.
  int ifirst;
  /// The number of slices on this processor.
  int nprocSlice; 
  /// Storage for this process's beads.
  Beads<NDIM> &beads;
  /// Buffers.
  Beads<NDIM> &buffer;
  /// Storage for this process's permutation.
  Permutation &permutation;  
  Permutation &globalPermutation;
  /// Storage for the inverse permutation.
  Permutation &inversePermutation;
  /// MPIManager;
  MPIManager &mpi;
  /// Number of slices on this processor.
  IArray npSlice;
};
#endif
