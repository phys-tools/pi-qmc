// $Id$
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#ifndef __PairDistance_h_
#define __PairDistance_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

typedef blitz::TinyVector<double,NDIM> Vec;
/// Base class for distance functions.
class PairDistance {public: 
  virtual double operator()(const Vec &r1, const Vec &r2, 
                            const SuperCell &cell)const=0;

};

  /// Distance taken from radial separation.
  class PairRadial : public PairDistance { public:
    PairRadial(int idir=-1) : mask(1.0) {if (idir!=-1) mask(idir)=0;}
    virtual double operator()(const Vec &r1, const Vec &r2, 
                              const SuperCell &cell) const {
      Vec delta=r1-r2; cell.pbc(delta);
      double radius2=0;
      for (int i=0; i<NDIM; ++i) radius2 += delta(i)*delta(i)*mask(i);
      return sqrt(radius2);
    }
    Vec mask;
  };
  /// Distance taken from cartesian position of particle 1.
  class PairCart1 : public PairDistance { public:
    PairCart1(int idim) : idim(idim){};
    int idim;
    virtual double operator()(const Vec &r1, const Vec &r2, 
                              const SuperCell &cell) const {return r1[idim];};
  };
  /// Distance taken from cartesian position of particle 2.
  class PairCart2 : public PairDistance { public:
    PairCart2(int idim) : idim(idim){};
    int idim;
    virtual double operator()(const Vec &r1, const Vec &r2, 
                              const SuperCell &cell)const {return r2[idim];};
  };
  /// Distance taken from cartesian separation.
  class PairCart : public PairDistance { public:
    PairCart(int idim) : idim(idim){};
    int idim;
    virtual double operator()(const Vec &r1, const Vec &r2, 
                              const SuperCell &cell)const {
      Vec delta=r1-r2; cell.pbc(delta);
      return delta[idim];
    }
  };
#endif
