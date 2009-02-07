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
#ifndef __SimInfoWriter_h_
#define __SimInfoWriter_h_

#include <string>
#include <iostream>
#include <vector>
#include "stats/EstimatorManager.h"
class SimulationInfo;

/// Class for writing out the SimulationInfo.
/// used in parsing and setup. 
/// @version $Revision$
/// @author John Shumway
class SimInfoWriter : public EstimatorManager::SimInfoWriter {
public: 
  /// Constructor.
  SimInfoWriter(const SimulationInfo &simInfo);
  /// Write to pimc.h5.
  virtual void writeH5(hid_t) const;
  //virtual void writeStdout(std::ostream&) const {};
  //virtual void writeAscii(std::ostream&) const {};
private:
  const SimulationInfo &simInfo;
};
#endif
