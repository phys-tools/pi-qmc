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
#ifndef __RandomNumGenerator_h_
#define __RandomNumGenerator_h_


#ifdef ENABLE_MPI
#include <mpi.h>
#define USE_MPI
#endif

#ifdef ENABLE_SPRNG
#define SIMPLE_SPRNG
#include <sprng.h>
#endif

#include <cstdlib>
#include <cmath>
#include <blitz/array.h>

/// The random number generator. 
/// @version $Revision$
/// @author John Shumway
class RandomNumGenerator{
public:
  ///
  static void seed(const int iseed) {
    int rank=0;
#ifdef ENABLE_MPI
    rank=MPI::COMM_WORLD.Get_rank();
#endif
#ifdef ENABLE_SPRNG
    init_sprng(DEFAULT_RNG_TYPE,iseed,SPRNG_DEFAULT);
#else
    //add 149*rank to the process number for parallel seeding
    srand(iseed+149*rank);
#endif
  }
  /// Get a random number.
#ifdef ENABLE_SPRNG
  static double getRand() {return sprng();}
#else
  static double getRand() {return rand()*SCALE_RAND;}
#endif

  /** Fill a blitz array with uniform random numbers. */
  static void makeRand(blitz::Array<double,1>& a) {
    for (int i=0; i<a.size(); ++i) a(i)=getRand();
  }

  /** Fill a blitz array with Gaussian distributed random
      numbers using Box-Mueller algorithm. */
  static void makeGaussRand(blitz::Array<double,1>& a) {
    makeGaussRand((double*)&a(0),a.size());
  }

  /** Fill a blitz TinyVector Array` with Gaussian distributed random
      numbers using Box-Mueller algorithm. */
  template <int TDIM>
  static void makeGaussRand(blitz::Array<blitz::TinyVector<double,TDIM>,1>& a) {
    makeGaussRand((double*)&a(0),TDIM*a.size());
  }

private:
  static const double SCALE_RAND;
  /** Fill a c-style array with Gaussian distributed random
      numbers using Box-Mueller algorithm. */
  static void makeGaussRand(double* a, const int n) {
    for (int i=0; i+1<n; i+=2) {
      double temp1=1-0.9999999999*getRand(), temp2=getRand();
      a[i]  =sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
      a[i+1]=sqrt(-2.0*log(temp1))*sin(6.283185306*temp2);
    }
    if (n%2==1) {
      double temp1=1-0.9999999999*getRand(), temp2=getRand();
      a[n-1]=sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
    }
  }

};
#endif
