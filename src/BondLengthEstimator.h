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
#ifndef __BondLengthEstimator_h_
#define __BondLengthEstimator_h_
#include "stats/ScalarEstimator.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>
#include <string>
#include <iostream>
class Paths;
class Action;
class SimulationInfo;
class Species;
class SuperCell;
/** Bond Length estimator.
 *  @todo add species to header (eg bond_length_a_b, etc)
 *  @version $Revision $ 
 *  @author Matthew Harowitz and John Shumway */
class BondLengthEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<std::string,1> SArray;
  /// Constructor.
  BondLengthEstimator(const SimulationInfo& simInfo, 
    const Species& species1, const Species& species2, 
    const std::string& unitName);
  /// Virtual destructor.
  virtual ~BondLengthEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int nfirst);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of bond-length estimator.
  virtual double calcValue() {return value=avlength/norm2;}
  /// Clear value of bond-length estimator.
  virtual void reset() {avlength=norm2=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// Bond length.
  double length;
  /// The normalization of the bond length (for multiple bonds).
  double norm1;
  /// The normalization of the bond length.
  double norm2;
  /// The Average length.
  double avlength;
  /// The names of the particles.
  SArray names;
  ///The species to measure the bond length of.
  std::string spec1;
  std::string spec2;
  /// The SuperCell, for wrapping periodic boundary conditions.
  const SuperCell &cell;
};

#endif
