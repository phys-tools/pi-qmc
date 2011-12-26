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
#ifndef __MainSpinParser_h_
#define __MainSpinParser_h_

#include "../src/XMLParser.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
class SimulationInfo;

/** Class for parsing xml data into a SimulationInfo object.
 *  @version $Revision$
 * @author John Shumway */
class MainSpinParser : XMLParser {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Constructor.
  MainSpinParser(SimulationInfo &simInfo) : simInfo(simInfo) {};
  /// Virtual destructor.
  virtual ~MainSpinParser();
  /// Parse some xml, storing the results in simulation info.
  virtual void parse(const xmlXPathContextPtr& ctxt);
  /// Get a reference to the simulation info, that remains valid
  /// until the parse method is called again.
  //SimulationInfo& getMainSpin() {return *simInfo;} 
private:
  /// The simulation info.
  SimulationInfo& simInfo;
};
#endif
