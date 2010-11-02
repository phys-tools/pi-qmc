// $Id$
/*  Copyright (C) 2004-2006, 2009 John B. Shumway, Jr.

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
#ifndef __WireEwald
#define __WireEwald
class SuperCell;
#include <blitz/tinyvec.h>
#include <blitz/array.h>

/** Class for evaluating long range Coulomb interaction at one dimension.
For a nano-wire with a screening layer, the long range Coulomb interaction 
is a well-behaved smooth function that can be approximated by a polynomial 
containing only even order terms, which enable us to significantly improve the 
speed of the calculation.
*/

class WireEwald {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::Array<double,NDIM> DArray;
  typedef blitz::Array<double,1> Array1;
  typedef blitz::Array<int,1> IArray;
  /// Constructor calculate the polynomial by fitting.
  WireEwald(const SuperCell&, const double d, const double rmax, 
            const Vec fitWidth, const IVec gridSize);
  /// Virtual destructor.
  virtual ~WireEwald();
  /// Return the value of the Coulomb interaction.
  double operator()(const Vec &);
  /// Get the maximum difference between the polynomail and the real interaction.
  double getMaxDiff();
private:
  /// Fit the Coulomb interaction by a polynomial, using Levenberg-Marquardt method.
  void curveFit();
  /// Polynomial function
#if NDIM==2
  inline double polynomial(int, int, double, double, const Array1&);
#endif
#if NDIM==3
  inline double polynomial(int, int, int, double, double, double, const Array1&);
#endif 
  /// The SuperCell.
  const SuperCell cell;
  /// The screening distance.
  const double d;
  /// The max range of Coulomb interaction.
  const double rmax;
  /// The min separation between particles.
  const double dx, dy;
#if NDIM==3
  const double dz;
#endif
  /// The size of grid for the Coulomb interaction.
  const IVec gridSize;
  /// Highest order x^NNy^MM.
  int NN,MM;
#if NDIM==3
  int LL;
#endif
  /// The coefficients of the polynomial.
  Array1 polyCoeff;
  /// The Coulomb potential and the corresponding coordinates.
  Array1 v,x,y;
#if NDIM==3
  Array1 z;
#endif
};
#endif
