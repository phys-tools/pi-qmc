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
#ifndef __Beads_h_
#define __Beads_h_
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <iostream>
#include <vector>
class Spin;
class SuperCell;
class BeadFactory;
#include "util/Permutation.h"

class BeadsBase {
public:
  BeadsBase(int npart, int nslice) 
  : nAuxBeads(1), auxBeads(nAuxBeads), npart(npart), nslice(nslice) {}
  typedef std::vector<BeadsBase*> BeadsArray;
  typedef blitz::Array<int,1> IArray;
  /// Clone method.
  virtual BeadsBase* clone(const int npart, const int nslice)=0;
  BeadsBase* cloneBase(const int npart, const int nslice) {
    BeadsBase* newBeads = auxBeads[0]->clone(npart,nslice);
    newBeads->nAuxBeads=auxBeads[0]->nAuxBeads;
    newBeads->auxBeads.resize(newBeads->nAuxBeads);
    newBeads->auxBeads[0]=newBeads;
    for (int i=1; i<nAuxBeads; ++i) {
      newBeads->auxBeads[i]=auxBeads[i]->clone(npart,nslice);
      newBeads->auxBeads[i]->auxBeads.resize(nAuxBeads);
    }
    for (int i=1; i<nAuxBeads; ++i) {
      for (int j=0; j<nAuxBeads; ++j) {
        newBeads->auxBeads[i]->auxBeads[j]=newBeads->auxBeads[j];
      }
    }
    return newBeads;
  }
  /// Virtual destructor.
  virtual ~BeadsBase() {}
  /// Copy operatation.
  BeadsBase& operator=(const BeadsBase& beadsIn) {
    for (int i=0; i<nAuxBeads; ++i) beadsIn.auxBeads[i]->vcopy(auxBeads[i]);
    return *this;
  }
  /// Copy a slice.
  void copySlice(int isliceIn, BeadsBase& beadsOut, int isliceOut) const {
    for (int i=0; i<nAuxBeads; ++i) {
      auxBeads[i]->vcopySlice(isliceIn, beadsOut.auxBeads[i], isliceOut);
    }
  }
  void copySlice(const Permutation& p, int isliceIn,
                 BeadsBase& beadsOut, int isliceOut) const {
    for (int i=0; i<nAuxBeads; ++i) {
      auxBeads[i]->vcopySlice(p, isliceIn, beadsOut.auxBeads[i], isliceOut);
    }
  }
  void copySlice(int isliceIn, BeadsBase& beadsOut, 
                 const Permutation& p, int isliceOut) const {
    for (int i=0; i<nAuxBeads; ++i) {
      auxBeads[i]->vcopySlice(isliceIn, beadsOut.auxBeads[i], p, isliceOut);
    }
  }
  void copySlice(const IArray& indexIn, int isliceIn, BeadsBase& beadsOut, 
                 const IArray& indexOut, int isliceOut) const {
    for (int i=0; i<nAuxBeads; ++i) {
      auxBeads[i]->vcopySlice(indexIn, isliceIn, beadsOut.auxBeads[i], 
                              indexOut, isliceOut);
    }
  }
  /// Get the number of auxiliary beads.
  int getAuxBeadCount() const {return nAuxBeads;}
  /// Get a pointer to an auxiliary bead.
  virtual const void* getAuxBead(const int ipart, const int islice, 
      const int iaux) const {
    return auxBeads[iaux]->vgetAuxBead(ipart,islice);
  }
  /// Get a pointer to an auxiliary bead.
  virtual void* getAuxBead(const int ipart, const int islice, 
      const int iaux) {
    return auxBeads[iaux]->vgetAuxBead(ipart,islice);
  }
  /// Support for printing.
  virtual void print(std::ostream& out) const=0;
  /// Return a pointer to the auxiliary beads.
  const BeadsBase* getAuxBeads(const int i=0) const {return auxBeads[i];}
  BeadsBase* getAuxBeads(const int i=0) {return auxBeads[i];}
  /// BeadFactory is a friend so that it can add to auxBeads.
  friend class BeadFactory;
protected:
  /// Number of auxBeads.
  int nAuxBeads;
  /// Auxiliary beads.
  BeadsArray auxBeads;  
  /// Number of particles.
  int npart;
  /// Number of slices.
  int nslice;
  virtual void vcopySlice(int isliceIn, 
                          BeadsBase* beadsOut, int isliceOut)const=0;
  virtual void vcopySlice(const Permutation& p, int isliceIn, 
                          BeadsBase* beadsOut, int isliceOut)const=0;
  virtual void vcopySlice(int isliceIn, BeadsBase* beadsOut,
                          const Permutation& p, int isliceOut)const=0;
  virtual void vcopySlice(const IArray& indexIn, int isliceIn, 
     BeadsBase* beadsOut, const IArray& indexOut, int isliceOut)const=0;
  virtual const void* vgetAuxBead(const int ipart, const int islice) const=0;
  virtual void* vgetAuxBead(const int ipart, const int islice)=0;
  virtual void vcopy(BeadsBase* beadsOut)const=0;
};

/// Storage for beads on the path.
/// @version $Revision$
/// @author John Shumway 
/// @bug Permute doesn't move auxilliary beads.
template <int TDIM>
class Beads : public BeadsBase{
public:
  /// Typedefs and constants.
  typedef blitz::TinyVector<double,TDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  typedef blitz::Array<Vec,2> VArray2;
  /// Constructor
  Beads(const int npart, const int nslice);
  /// Copy constructor.
  Beads(const Beads&);
  /// Clone method.
  virtual BeadsBase* clone(const int npart, const int nslice) {
    return new Beads(npart,nslice);
  }
  /// Virtual destructor
  ~Beads() {
    // Destroy the other beads if this is the main set of beads.
    if (nAuxBeads>0 && auxBeads[0]==this) {
      for (int i=1; i<nAuxBeads; ++i) {
        delete(auxBeads[i]); auxBeads[i]=0;
      }
    }
    if (nAuxBeads>0) auxBeads[0]=0;
  };
  /// Support for printing.
  virtual void print(std::ostream& out) const;
  /// Const access
  const Vec operator() (const int ipart, const int islice) const {
    return coord(ipart,islice);
  }
  /// Non const access
  Vec& operator() (const int ipart, const int islice) {
    return coord(ipart,islice);
  }
  /// Get a relative displacement.
  const Vec delta(const int ipart, const int islice, const int istep) const {
    return coord(ipart,islice)-coord(ipart,islice+istep);
  }
  /// Return the number of particles.
  int getNPart() const {return npart;}
  /// Return the number of slices.
  int getNSlice() const {return nslice;}
  /// Set the supercell.
  void setSuperCell(const SuperCell* c) {cell = c;}
  /// Return a pointer to the supercell.
  const SuperCell* getSuperCell() const {return cell;}
  /// Permute the beads.
  void permute(const Permutation&);
  /// Return a reference to the coordinate array (for MPI calls).
  VArray2& getCoordArray() {return coord;}
protected:
  /// Virtual method to copy a slice.
  virtual void vcopySlice(int isliceIn, 
                          BeadsBase* beadsOut, int isliceOut) const {
    Beads& beads = *dynamic_cast<Beads*>(beadsOut);
    for (int i=0; i<npart; ++i) beads.coord(i,isliceOut)=coord(i,isliceIn);
  }
  virtual void vcopySlice(const Permutation& p, int isliceIn, 
                          BeadsBase* beadsOut, int isliceOut) const {
    Beads& beads = *dynamic_cast<Beads*>(beadsOut);
    for (int i=0; i<npart; ++i) beads.coord(i,isliceOut)=coord(p[i],isliceIn);
  }
  virtual void vcopySlice(int isliceIn, BeadsBase* beadsOut, 
                          const Permutation& p, int isliceOut) const {
    Beads& beads = *dynamic_cast<Beads*>(beadsOut);
    for (int i=0; i<npart; ++i) beads.coord(p[i],isliceOut)=coord(i,isliceIn);
  }
  virtual void vcopySlice(const IArray& indexIn, int isliceIn, 
      BeadsBase* beadsOut, const IArray& indexOut, int isliceOut)const {
    Beads& beads = *dynamic_cast<Beads*>(beadsOut);
    for (int i=0; i<indexIn.size(); ++i) {
      beads.coord(indexOut(i),isliceOut)=coord(indexIn(i),isliceIn);
    }
  }
  virtual const void* vgetAuxBead(const int ipart, const int islice) const {
    return &coord(ipart,islice);}
  virtual void* vgetAuxBead(const int ipart, const int islice) {
    return &coord(ipart,islice);}
  virtual void vcopy(BeadsBase* beadsOut) const {
    Beads& beads = *dynamic_cast<Beads*>(beadsOut);
    beads.coord=coord;
  }
private:
  /// Storage for the bead coordinates.
  VArray2 coord;
  /// The super cell.
  const SuperCell* cell;
};



template <int TDIM>
Beads<TDIM>::Beads(const int npart, const int nslice) 
 : BeadsBase(npart, nslice),
   coord(npart,nslice,blitz::ColumnMajorArray<2>()), cell(0) {
 auxBeads[0]=this;
 coord=0.0;
}

template <int TDIM>
Beads<TDIM>::Beads(const Beads& b) 
 : BeadsBase(b.npart,b.nslice), coord(b.coord.copy()), cell(b.cell) {
}

template <int TDIM>
std::ostream& operator<<(std::ostream& out, const Beads<TDIM>& p) {
  p.print(out);
  return out;
}

template <int TDIM>
void Beads<TDIM>::print(std::ostream& out) const {
  out << "Beads (npart=" << npart << ", nslice=" << nslice << ")" << std::endl;
  for (int islice=0; islice<nslice; ++islice) {
    out << "  Slice " << islice << ":" << std::endl;
    for (int ipart=0; ipart<npart; ++ipart)  {
      out << "    Particle " << ipart << ": "
          << coord(ipart,islice) << std::endl;
    }
  }
}

template <int TDIM>
void Beads<TDIM>::permute(const Permutation& p) {
  VArray buffer(npart);
  for (int islice=0; islice<nslice; ++islice) {
    // Swap paths.
    for (int i=0; i<npart; ++i) buffer(i)=coord(p[i],islice);
    for (int i=0; i<npart; ++i) coord(i,islice)=buffer(i);
  }
}
#endif
