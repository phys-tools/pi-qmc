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
#ifndef __DoubleParallelPaths_h_
#define __DoubleParallelPaths_h_
template <int TDIM> class Beads;
class Permutation;
class MPIManager;
class BeadFactory;
class SuperCell;
#include "Paths.h"
#include <blitz/array.h>

/// Storage for paths in a serial process.
/// @version $Revision$
/// @author John Shumway
class DoubleParallelPaths : public Paths {
public:
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  DoubleParallelPaths(int npart, int nslice, double tau, const SuperCell& cell,
                MPIManager &mpi, const BeadFactory&);
  /// Destructor.
  virtual ~DoubleParallelPaths();
  /// Loop over links, calling a LinkSummable object.
  virtual void sumOverLinks(LinkSummable&) const;
  /// Get a reference to a bead.
  virtual Vec& operator()(const int ipart, const int islice);
  /// Get a const reference to a bead.
  virtual const Vec& operator()(const int ipart, const int islice) const;
  /// Get a reference to a bead by offset.
  virtual Vec& operator()(const int ipart, const int islice, const int istep);
  /// Get a const reference to a bead by offset.
  virtual const Vec&
    operator()(const int ipart, const int islice, const int istep) const;
  /// Get a relative displacement a bead by offset.
  virtual Vec delta(const int ipart, const int islice, const int istep) const;
  /// Get beads.
  virtual void getBeads(int ifirstSlice, Beads<NDIM>& ) const;

  /// Get a slice.
  virtual void getSlice(int islice, VArray& ) const;
  /// Get the number of slices on this processor.
  virtual int getNProcSlice(){return nprocSlice;}
  /// Get the number of slices unique to this processor.
  virtual int getNUniqueSlice(){return nprocSlice-1;}
  /// Get auxialiary bead.
  virtual const void* getAuxBead(const int ipart, const int islice, 
                                 const int iaux) const;
  /// Get auxialiary bead.
  virtual void* getAuxBead(const int ipart, const int islice, 
                                 const int iaux);
  /// Put beads.
  virtual void putBeads(int ifirstSlice,
                        const Beads<NDIM>&, const Permutation&) const;
  virtual void putDoubleBeads(
                 int ifirstSlice1,Beads<NDIM>&, Permutation&,
                 int ifirstSlice2,Beads<NDIM>&, Permutation&) const;
  /// Get the global permuation.
  /// @bug Placeholder, not correct for MPI.
  virtual   const Permutation& getGlobalPermutation() const;

  virtual const Permutation& getPermutation() const {return permutation1;}
  virtual int getLowestOwnedSlice(bool d) const {return ifirst+1;}
  virtual int getHighestOwnedSlice(bool d) const { return ifirst+nprocSlice;}
  virtual int getHighestSampledSlice(int n, bool d) const {
    return ifirst+nprocSlice-n;}
  virtual bool isOwnedSlice(int islice) const {
    const int jslice = (islice-ifirst+1)%(nslice/2)+ifirst-1;
    return jslice>ifirst && jslice <= ifirst+nprocSlice;}
  virtual void shift(const int ishift);
  virtual void setBuffers();
  virtual bool isDouble() const {return true;}
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
  Beads<NDIM> &beads1, &beads2;
  /// Buffers.
  Beads<NDIM> &buffer1, &buffer2;
  /// Storage for this process's permutation.
  Permutation &permutation1, &permutation2;
  Permutation &globalPermutation;
  /// Storage for the inverse permutation.
  Permutation &inversePermutation1, &inversePermutation2;
  /// MPIManager;
  MPIManager &mpi;
  /// Number of slices on each of the processors.
  IArray npSlice;
};
#endif
