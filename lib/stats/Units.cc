//$Id$
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
#include "Units.h"
#include <iostream>

Units::Units(const std::string& eunit, const std::string& lunit)
  : internalEnergyUnit(eunit), internalLengthUnit(lunit) {
  /// To set these values, express the atomic unit in terms of named unit.
  ///
  /// Set conversion from bohr radius to common length units.
  lengthOut["nm"]=0.05291772086;
  lengthOut["A"] =0.5291772086;
  lengthOut["pm"]=52.91772086;
  lengthOut["um"]=0.05291772086e-3;
  lengthOut["a0"]=1.0;
  /// Set conversion from Hartree to common energy units.
  energyOut["eV"]=27.2113845;
  energyOut["meV"]=27211.3845;
  energyOut["ueV"]=27211384.5;
  energyOut["K"]=3.1577465e5;
  energyOut["mK"]=3.1577465e8;
  energyOut["uK"]=3.1577465e11;
  energyOut["nK"]=3.1577465e14;
  energyOut["Ha"]=1.0;
  energyOut["Ry"]=2.0;
  energyOut["cm-1"]=219474.6314;
  energyOut["kcal/mol"]=627.5095;
  energyOut["Hz"]=6.579683920e15;
  energyOut["kHz"]=6.579683920e12;
  energyOut["MHz"]=6.579683920e9;
  energyOut["GHz"]=6.579683920e6;
  energyOut["THz"]=6.579683920e3;
  energyOut["PHz"]=6.579683920;
  /// Set conversion from electron mass to common mass units.
  massOut["amu"]=9.1093826e-31/1.66053886e-27;
  massOut["m_p"]=9.1093826e-31/1.67262171e-27;
  massOut["m_e"]=1.0;
}

double Units::getLengthScaleIn(const std::string& unit) const {
  unitMap::const_iterator i=lengthOut.find(unit);
  if (i==lengthOut.end()) {
    return 1; //Unit not found, may want to throw error.
  } else {
    return 1.0/i->second;
  }
}

double Units::getEnergyScaleIn(const std::string& unit) const {
  unitMap::const_iterator i=energyOut.find(unit);
  if (i==energyOut.end()) {
    return 1; //Unit not found, may want to throw error.
  } else {
    return 1.0/i->second;
  }
}

double Units::getMassScaleIn(const std::string& unit) const {
  unitMap::const_iterator i=massOut.find(unit);
  if (i==massOut.end()) {
    return 1; //Unit not found, may want to throw error.
  } else {
    return 1.0/i->second;
  }
}
