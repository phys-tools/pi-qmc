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
#include "SpinDemo.h"
#include <fstream>

void SpinDemo::generate() const {
  std::string filename("spin.xml");
  std::cout << "Writing spin demo file: " << filename << "\n" << std::endl;

  std::string xmlstring="\
<?xml version=\"1.0\"?>\n\
<Simulation>\n\
  <SuperCell a=\"5\" x=\"1\" y=\"1\" z=\"1\"/>\n\
  <Species name=\"e\" count=\"1\" mass=\"1.0\" charge=\"-1\" \
isStatic=\"true\"/>\n\
  <Temperature value=\"1.0\" tau=\"0.0004\"/>\n\
  <Spin omega=\"50\"/>\n\
  <Action>\n\
    <SpringAction/>\n\
    <SpinAction bx=\"0\" by=\"0\" bz=\"0.50\"/>\n\
    <SpinFixedPhaseAction model=\"SpinPhase\" bz=\"0.50\"/>\n\
  </Action>\n\
  <Estimators>\n\
    <ThermalEnergyEstimator/>\n\
    <SpinEstimator/>\n\
  </Estimators>\n\
  <PIMC>\n\
    <RandomGenerator iseed=\"211\"/>\n\
    <SetSpin/>\n\
    <!-- Thermalize -->\n\
    <Loop nrepeat=\"1000\">\n\
      <ChooseSection nlevel=\"7\">\n\
         <Sample npart=\"1\" nrepeat=\"10\" mover=\"Spin\"/>\n\
      </ChooseSection>\n\
    </Loop>\n\
    <!-- Sample data -->\n\
    <Loop nrepeat=\"100\">\n\
      <Loop nrepeat=\"50\">\n\
        <Loop nrepeat=\"31\">\n\
          <ChooseSection nlevel=\"5\">\n\
            <Sample npart=\"1\" mover=\"Spin\"/>\n\
            <Sample npart=\"1\" mover=\"Spin\" both=\"true\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Loop nrepeat=\"7\">\n\
          <ChooseSection nlevel=\"7\">\n\
            <Sample npart=\"1\" nrepeat=\"20\" mover=\"Spin\"/>\n\
            <Sample npart=\"1\" nrepeat=\"20\" mover=\"Spin\" both=\"true\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Loop nrepeat=\"4\">\n\
          <ChooseSection nlevel=\"8\">\n\
            <Sample npart=\"1\" nrepeat=\"50\" mover=\"Spin\"/>\n\
            <Sample npart=\"1\" nrepeat=\"50\" mover=\"Spin\" both=\"true\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Measure estimator=\"all\"/>\n\
      </Loop>\n\
      <Collect estimator=\"all\"/>\n\
      <WritePaths file=\"paths.out\"/>\n\
    </Loop>\n\
  </PIMC>\n\
</Simulation>";

  std::ofstream file(filename.c_str());
  file << xmlstring << std::endl;
}
