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
#ifndef __SimulationInfo_h_
#define __SimulationInfo_h_

#include <string>
#include <iostream>
#include <vector>
class SuperCell;
class Units;
class BeadFactory;
#include "Species.h"

/// Class for storing the description of the simulation,
/// used in parsing and setup. 
/// @version $Revision$
/// @author John Shumway
class SimulationInfo {
public:
  /// Constructor
  SimulationInfo();
  /// Constructor
  SimulationInfo(SuperCell*, int npart, std::vector<Species*>&,
    std::vector<Species*>&, double temperature, double tau, int nslice);
  /// Destructor.
  virtual ~SimulationInfo();
  /// Get a pointer to the supercell.
  SuperCell* getSuperCell() const {return superCell;}
  /// Get the number of species.
  int getNSpecies() const {return speciesList.size();}
  /// Get species by species index.
  const Species& getSpecies(const int i) const {return *speciesList[i];}
  /// Get species by particle index.
  const Species& getPartSpecies(const int i) const {return *speciesIndex[i];}
  /// Get species by name.
  const Species& getSpecies(const std::string&) const;
  /// Get the number of particles.
  int getNPart() const {return npart;}
  /// Get the temperature value.
  double getTemperature() const {return temperature;}
  /// Get the time step.
  double getTau() const {return tau;}
  void setTau(double value) {tau = value;}
  /// Get the number of slices.
  int getNSlice() const {return nslice;}
  /// Get information about units.
  const Units* getUnits() const {return units;}
  /// Get a reference to the BeadFactory.
  BeadFactory& getBeadFactory() const {return beadFactory;}
  /// The energy scale for spinors.
  double spinOmega;
private:
  /// The supercell.
  SuperCell* superCell;
  /// Unit information.
  Units* units;
  /// The number of particles.
  int npart;
  /// The species list.
  std::vector<Species*> speciesList;
  /// Index from each particle to it's species.
  std::vector<Species*> speciesIndex;
  /// The temperature.
  double temperature;
  /// The timestep.
  double tau;
  /// The number of slices.
  int nslice;
  /// The bead factory.
  BeadFactory& beadFactory;
friend class SimInfoParser;
friend class SimInfoWriter;
friend std::ostream& operator<<(std::ostream&, const SimulationInfo&);
};
#endif
