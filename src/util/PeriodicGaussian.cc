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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PeriodicGaussian.h"
#include <iostream>

PeriodicGaussian::PeriodicGaussian(const double a, const double d, const int n)
  : a(a), d(d), grid(n+1), dx(d/(2*n)), dxInv(1.0/dx) {
  grid=0.0;
  // First set grid to function value and first derivative.
  //const int nImage=0;
  const int nImage=(int)(10.0/(sqrt(a)*d));
  for (int iImage=-nImage; iImage<=nImage; ++ iImage) {
    for (int i=0; i<=n; ++i) {
      double x=i*dx;
      double expax2 = exp(-a*pow((x-iImage*d),2));
      grid(i)[0]+= expax2;
      grid(i)[1]+= -2*a*(x-iImage*d)*expax2;
    }
  }
  // Then compute cubic polynomial coefficients by matching function
  // values and first derivatives.
  grid(0)[1]=0; grid(n)[1]=0;
  for (int i=0; i<n; ++i) {
    Vec4 &f0(grid(i)), &f1(grid(i+1));
    f0[3]=((f1[1]+f0[1])*dx-2*(f1[0]-f0[0]))*dxInv*dxInv*dxInv;
    f0[2]=(f1[0]-f0[0]-dx*(f0[1]+f0[3]*dx*dx))*dxInv*dxInv;
  }
}


/*int main(int argc, char** argv) {
  std::cout << "#Periodic gausian test" << std::endl;
  PeriodicGaussian pg(0.01,10,3);
  for (int i=0; i<=100; ++i) {
    double x=i*0.05; std::cout << x << " " << pg(x) << std::endl;
  }
  return 0;
}*/
