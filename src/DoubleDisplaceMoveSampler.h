// $Id: DoubleDisplaceMoveSampler.h 22  2009-05-18 Saad Khairallah $
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
#ifndef __DoubleDisplaceMoveSampler_h_
#define __DoubleDisplaceMoveSampler_h_
template <int TDIM> class Beads;
class SuperCell;
class Action;
class DoubleAction;
class Mover;
class Paths;
class PathsChooser;
class ParticleChooser;
//class PermutationChooser;
class Permutation;
class AccRejEstimator;
class BeadFactory;
class MPIManager;

#include "Algorithm.h"
#include <vector>
#include <blitz/array.h>
#include <iostream>
#include "DisplaceMoveSampler.h"

class DoubleDisplaceMoveSampler : public DisplaceMoveSampler {
public:
  typedef blitz::Array<int,1> IArray;
  DoubleDisplaceMoveSampler(const int nmoving, Paths&,  const double dist, 
			    const double freq, ParticleChooser&, Mover&, Action*, 
			    const int nrepeat, const BeadFactory&, const MPIManager* mpi, 
			    DoubleAction* doubleAction);
  virtual ~DoubleDisplaceMoveSampler();
  virtual void run();
  Beads<NDIM>& getMovingBeads(const int i)  {
    return i==1?*movingBeads1:*movingBeads2;}
  const Beads<NDIM>& getMovingBeads(const int i) const  {
    return i==1?*movingBeads1:*movingBeads2;}
  Beads<NDIM>& getPathsBeads(const int i) {
    return i==1?*pathsBeads1:*pathsBeads2;}
  const Beads<NDIM>& getPathsBeads(const int i) const {
    return i==1?*pathsBeads1:*pathsBeads2;}
  
  IArray&  getMovingIndex() {return *movingIndex;}
  const IArray&  getMovingIndex() const {return *movingIndex;}
  void setAction(Action*, const int level=0); /// not needed
  const SuperCell& getSuperCell() const {return cell;}
  double getDist() {return dist;}
  int getNSlice() {return nslice;}
  int getNSlice() const {return nslice;}
  void activateSection(const int i);
  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);

protected:
  bool tryMove(int imovingNonPerm);
  
  Beads<NDIM> *pathsBeads1;
  Beads<NDIM> *pathsBeads2;
  Beads<NDIM> *movingBeads1;
  Beads<NDIM>  *movingBeads2;
  DoubleAction *doubleAction;
  AccRejEstimator* accRejEst;
  IArray *movingIndex2;
  int iFirstSlice, iFirstSlice2;
};
#endif
