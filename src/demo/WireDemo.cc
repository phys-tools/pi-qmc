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
#include "WireDemo.h"
#include <fstream>
#include <config.h>
void WireDemo::generate() const {
  std::cout << "Writing quantum wire input file: wire.xml\n" << std::endl;

  std::string xmlstring="\
<?xml version=\"1.0\"?>\n\
<Simulation>\n\
  <SuperCell a=\"300 nm\" x=\"1\" y=\"1\"/>\n\
  <Species name=\"edn\" count=\"10\" mass=\"0.0667 m_e\" charge=\"-1\"/>\n\
  <Temperature value=\"12 K\" nslice=\"200\"/>\n\
  <Action>\n\
    <SpringAction/>\n";

  if (NDIM>1) {
    xmlstring+="\
    <SHOAction omega=\"30 meV\" ndim=\"1\"/>\n";
  }
  xmlstring+="\
    <FixedNodeAction species=\"edn\" model=\"WireNodes\" omega=\"30 meV\"\n\
                     noNodalAction=\"true\"/>\n\
    <!--CoulombAction norder=\"0\" rmax=\"800 nm\" useEwald=\"true\" ewaldNDim=\"1\"\n\
               kcut=\"30\" screenDist=\"100 nm\" epsilon=\"12\" dumpFiles=\"true\"/-->\n\
  </Action>\n\
  <Estimators>\n\
    <ThermalEnergyEstimator/>\n\
    <ConductivityEstimator nbin=\"50\" ndbin=\"25\"/>\n\
    <!--PairCFEstimator species1=\"edn\" species2=\"edn\" name=\"g_edn_edn\"/-->\n\
  </Estimators>\n\
  <PIMC nclone=\"4\" nworker=\"1\">\n\
    <RandomGenerator/>\n\
    <SetCubicLattice nx=\"10\" ny=\"1\" a=\"30 nm\"/>\n\
    <!--ReadPaths file=\"paths.in\"/-->\n\
    <!-- Thermalize -->\n\
    <Loop nrepeat=\"500\">\n\
      <ShiftWorkers maxShift=\"99\">\n\
        <ChooseSection nlevel=\"7\">\n\
          <Sample npart=\"1\" nrepeat=\"10\" mover=\"Free\" species=\"edn\"/>\n\
        </ChooseSection>\n\
      </ShiftWorkers>\n\
    </Loop>\n\
    <!-- Sample data -->\n\
    <Loop nrepeat=\"1000\">\n\
      <Loop nrepeat=\"20\">\n\
        <ShiftWorkers maxShift=\"99\">\n\
          <ChooseSection nlevel=\"7\">\n\
            <Sample npart=\"1\" nrepeat=\"40\" mover=\"Free\"/>\n\
            <Sample npart=\"2\" nrepeat=\"1000\" mover=\"Free\"/>\n\
            <Sample npart=\"3\" nrepeat=\"5000\" mover=\"Free\"/>\n\
          </ChooseSection>\n\
          <Loop nrepeat=\"4\">\n\
            <ChooseSection nlevel=\"5\">\n\
              <Sample npart=\"1\" nrepeat=\"10\" mover=\"Free\"/>\n\
              <Sample npart=\"2\" nrepeat=\"20\" mover=\"Free\"/>\n\
              <Sample npart=\"3\" nrepeat=\"50\" mover=\"Free\"/>\n\
              <Sample npart=\"4\" nrepeat=\"100\" mover=\"Free\"/>\n\
            </ChooseSection>\n\
          </Loop>\n\
        </ShiftWorkers>\n\
        <Measure estimator=\"all\"/>\n\
      </Loop>\n\
      <Collect estimator=\"all\"/>\n\
      <WritePaths file=\"paths.out\"/>\n\
    </Loop>\n\
  </PIMC>\n\
</Simulation>";

  std::ofstream file("wire.xml");
  file << xmlstring << std::endl;
}
