//$Id$
/*  Copyright (C) 2010 John B. Shumway, Jr.

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

#ifndef __EMARateEstimator_h_
#define __EMARateEstimator_h_
#include "stats/ScalarEstimator.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>
class Paths;
class SimulationInfo;
class SuperCell;
/** 
 *  @author John Shumway  */
class EMARateEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  EMARateEstimator(const SimulationInfo& simInfo, double C);
  virtual ~EMARateEstimator() {}
  virtual void initCalc(const int nslice, const int firstSlice);
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  virtual void endCalc(const int nslice);
  virtual double calcValue() {return sum/norm;}
  virtual void reset() {sum=norm=0.;}
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  const double dtau;
  const double masse;
  const double massh;
  const double C;
  const SuperCell& cell;
  double actionDifference;
  double sum, norm;
};

#endif
