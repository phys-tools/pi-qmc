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
#ifndef __LinkSummable_h_
#define __LinkSummable_h_
#include "stats/Estimator.h"
#include <blitz/tinyvec.h>
class Paths;
/// Interface class for ojects that can be summed over links.
/// @version $Revision$
/// @author John Shumway
class LinkSummable {
public:
typedef blitz::TinyVector<double,NDIM> Vec;
  /// Virtual destructor.
  virtual ~LinkSummable() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice) {};
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          int ipart, int islice, const Paths&)=0;
  /// Finalize the calculation.
  virtual void endCalc(int nslice) {};
};
#endif
