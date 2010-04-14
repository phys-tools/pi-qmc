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
#ifndef __PeriodicGaussian_h_
#define __PeriodicGaussian_h_
#include <blitz/array.h>

/** Function object for a periodic gaussian.
  * @version $Revision$
  * @author John Shumway */
class PeriodicGaussian {
public:
  typedef blitz::TinyVector<double,4> Vec4;
  typedef blitz::Array<Vec4,1> V4Array;
  /// Constructor.
  PeriodicGaussian(const double a, const double d, const int n);
  /// Get the Gaussian value.
  inline double operator()(const double x) const {
    const int i=(int)(x*dxInv); 
    const double r=x-i*dx;
    const Vec4 &f=grid(i);
    return f[0]+r*(f[1]+r*(f[2]+r*f[3]));
  }
  /// Get the gradient value.
  double grad(const double x) const {
    const int i=(int)(x*dxInv);
    const double r=x-i*dx;
    const Vec4 &f=grid(i);
    return f[1]+r*(2.0*f[2]+r*3.0*f[3]);
  }
  /// Get the second derivative value.
  double d2(const double x) const {
    const int i=(int)(x*dxInv); 
    const double r=x-i*dx;
    const Vec4 &f=grid(i);
    return 2.0*f[2]+r*6.0*f[3];
  }
 
private:
  /// The coefficient.
  const double a;
  /// The period.
  const double d;
  /// The grid.
  V4Array grid;
  /// The grid spacing.
  double dx, dxInv;
};
#endif
