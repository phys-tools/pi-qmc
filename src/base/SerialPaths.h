// $Id$
/*  Copyright (C) 2004-2006,2009 John B. Shumway, Jr.

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
#ifndef __SerialPaths_h_
#define __SerialPaths_h_
template <int TDIM> class Beads;
class Permutation;
class SuperCell;
class BeadFactory;
#include "Paths.h"

/// Storage for paths in a serial process.
/// We own slices 0 through nslice-1.
/// Access methods accept slices -nslice+1 through 2*nslice-1.
/// @version $Revision$
/// @author John Shumway
class SerialPaths : public Paths {
public:
  /// Constructor.
  SerialPaths(int npart, int nslice, double tau, const SuperCell& cell,
              const BeadFactory&);
  /// Destructor.
  virtual ~SerialPaths();
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
  /// Get the global permuatation.
  virtual const Permutation& getPermutation() const {return permutation;}
  virtual const Permutation& getGlobalPermutation() const {return permutation;}
  virtual int getLowestOwnedSlice(bool d) const {return 0;}
  virtual int getHighestOwnedSlice(bool d) const {
    return d?nslice/2-1:nslice - 1;}
  virtual int getHighestSampledSlice(int n, bool d) const {
    return (d?nslice/2-1:nslice - 1);}
  virtual bool isOwnedSlice(int islice) const {
    return islice>=0 && islice<nslice;}
  virtual void shift(const int ishift);
  virtual void clearPermutation();
private:
  void putBeads(int ifirstSlice, const Beads<NDIM>&, const Permutation&, 
                int ifirst, int nbslice) const;
  /// Storage for the beads.
  Beads<NDIM>& beads;
  /// Storage shifting uffers.
  Beads<NDIM> &buffer1,&buffer2;
  /// Storage for the permutation.
  Permutation& permutation;
  /// Storage for the inverse permutation.
  Permutation& inversePermutation;
};
#endif
