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
#ifndef __BlitzArrayBlkdEst_h_
#define __BlitzArrayBlkdEst_h_
#include "ArrayBlockedEstimator.h"
#include "EstimatorReportBuilder.h"
#include "MPIManager.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
class Paths;
/// Base class for scalar estimators.
/// @version $Revision$
/// @author John Shumway
template <int N>
class BlitzArrayBlkdEst : public ArrayBlockedEstimator {
public:
  typedef blitz::Array<float,N> ArrayN;
  typedef blitz::TinyVector<int,N> IVecN;
  typedef blitz::TinyVector<double,N> VecN;
  /// Construct by giving the name and size.
  BlitzArrayBlkdEst(const std::string &name, const std::string &typeString,
                    const IVecN &n, bool hasError);
  /// Virtual destructor.
  virtual ~BlitzArrayBlkdEst() {}
  /// Clear the value.
  virtual void reset()=0;
  /// Average over clones.
  virtual void averageOverClones(const MPIManager* mpi);
  /// Get the number of dimensions in the array.
  virtual int getNDim() const {return rank;}
  /// Get the extent of the array in dimensions idim.
  virtual int getExtent(const int idim) const {return n(idim);}
  /// Get a pointer to the data.
  virtual const void* getData() const {return accumvalue.data();}
  virtual const void* getError() const {
    return hasErrorFlag?accumvalue2.data():(float*)0;
  }
  virtual void normalize() const {
    accumvalue/=accumnorm; //Temporarily normalize the data to write it out.
    if (hasErrorFlag) { //Temperorily convert accumvalue2 into rms error.
      accumvalue2 /= accumnorm;
      accumvalue2 -= accumvalue*accumvalue;
      accumvalue2 = sqrt(fabs(accumvalue2/iblock));
    }
  }
  virtual void unnormalize() const {
    if (hasErrorFlag) {
      accumvalue2 *= accumvalue2*iblock;
      accumvalue2 += accumvalue*accumvalue;
      accumvalue2 *= accumnorm;
    }
    accumvalue*=accumnorm;
  }
  virtual bool hasScale() const {return scale!=0;}
  virtual bool hasOrigin() const {return origin!=0;}
  virtual const void* getScale() const {return scale;}
  virtual const void* getOrigin() const {return origin;}
  virtual bool hasMin() const {return min!=0;}
  virtual bool hasMax() const {return max!=0;}
  virtual const void* getMin() const {return min;}
  virtual const void* getMax() const {return max;}
protected:
  static const int rank=N;
  const IVecN n;
  mutable ArrayN value,accumvalue,accumvalue2;
  double norm, accumnorm;
  bool hasError;
  int iblock;
  VecN* scale;
  VecN* origin;
  VecN* min;
  VecN* max;
};



template <int N>
BlitzArrayBlkdEst<N>::BlitzArrayBlkdEst(const std::string& name, 
  const std::string& typeString, const IVecN& n, bool hasError)
 : ArrayBlockedEstimator(name,typeString,hasError),
   n(n), value(n), accumvalue(n),
   accumvalue2(hasError?n:n*0), norm(0), accumnorm(0),
   hasError(hasError), iblock(0), scale(0), origin(0), min(0), max(0) {
  value=0; accumvalue=0; accumvalue2=0;
}

template<int N>
void BlitzArrayBlkdEst<N>::averageOverClones(const MPIManager* mpi) {
#ifdef ENABLE_MPI
  if (mpi && mpi->isCloneMain()) {
    int rank = mpi->getCloneComm().Get_rank();
    int size = mpi->getCloneComm().Get_size();
    if (size>1) {
      reset();
      if (rank==0) {
#if MPI_VERSION==2
        mpi->getCloneComm().Reduce(MPI::IN_PLACE,&norm,1,MPI::DOUBLE,
                                   MPI::SUM,0);
        mpi->getCloneComm().Reduce(MPI::IN_PLACE,value.data(),value.size(),
                                   MPI::FLOAT,MPI::SUM,0);
#else
        double nbuff;
        ArrayN vbuff(n);
        mpi->getCloneComm().Reduce(&norm,&nbuff,1,MPI::DOUBLE,
                                   MPI::SUM,0);
        mpi->getCloneComm().Reduce(value.data(),vbuff.data(),value.size(),
                                   MPI::FLOAT,MPI::SUM,0);
        norm=nbuff;
        value=vbuff;
#endif
      } else {
        mpi->getCloneComm().Reduce(&norm,NULL,1,MPI::DOUBLE,MPI::SUM,0);
        mpi->getCloneComm().Reduce(value.data(),NULL,value.size(),
                                   MPI::FLOAT,MPI::SUM,0);
      }
    }
  }
#endif
  // Next add value to accumvalue and accumvalue2.
  accumvalue += value/norm;
  if (hasErrorFlag) accumvalue2 += (value*value)/(norm*norm);
  accumnorm+=1.;
  value=0.; norm=0;
  ++iblock;
}

#endif
