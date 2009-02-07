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
#ifndef __Units_h_
#define __Units_h_
#include <string>
#include <map>
class Paths;
/// Class for unit translations.
/// @version $Revision$
/// @author John Shumway
class Units {
public:
  /// Construct by giving the name.
  Units(const std::string& eunit, const std::string& lunit);
  /// Destructor.
  ~Units();
  /// Get conversion factor for length input with unit.
  double getLengthScaleIn(const std::string& unit) const;
  /// Get conversion factor for energy input with unit.
  double getEnergyScaleIn(const std::string& unit) const;
  /// Get conversion factor for mass input with unit.
  double getMassScaleIn(const std::string& unit) const;
protected:
  typedef std::map<std::string,double> unitMap;
  /// The name of the internal energy unit.
  const std::string internalEnergyUnit;
  /// The name of the internal length unit.
  const std::string internalLengthUnit;
  /// The conversion factors to different length units.
  unitMap lengthOut;
  /// The conversion factors to different energy units.
  unitMap energyOut;
  /// The conversion factors to different mass units.
  unitMap massOut;
};
#endif
