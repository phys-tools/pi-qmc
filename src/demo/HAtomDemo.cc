#include "HAtomDemo.h"
#include <string>
#include <fstream>

void HAtomDemo::generate() const {
  std::cout << "Writing hydrogen atom input file: hatom.xml\n" << std::endl;

  std::string xmlstring="\
<?xml version=\"1.0\"?>\n\
<Simulation>\n\
  <SuperCell a=\"25 A\" x=\"1\" y=\"1\" z=\"1\"/>\n\
  <Species name=\"eup\" count=\"1\" mass=\"1 m_e\" charge=\"-1\"/>\n\
  <Species name=\"H\" count=\"1\" mass=\"1.006780 amu\" charge=\"+1\"/>\n\
  <Temperature value=\"2400 K\" nslice=\"2048\"/>\n\
  <Action>\n\
    <SpringAction/>\n\
    <CoulombAction norder=\"3\" rmax=\"8.0\" ngridPoints=\"1000\"\n\
                   dumpFiles=\"false\"/>\n\
  </Action>\n\
  <Estimators>\n\
    <ThermalEnergyEstimator unit=\"eV\"/>\n\
    <VirialEnergyEstimator unit=\"eV\" nwindow=\"512\"/>\n\
  </Estimators>\n\
  <PIMC>\n\
    <RandomGenerator/>\n\
    <SetCubicLattice nx=\"2\" ny=\"2\" nz=\"2\" a=\"1.\"/>\n\
    <!-- Thermalize -->\n\
    <Loop nrepeat=\"10000\">\n\
      <ChooseSection nlevel=\"7\">\n\
         <Sample npart=\"1\" mover=\"Free\" species=\"H\"/>\n\
         <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
      </ChooseSection>\n\
    </Loop>\n\
    <!-- Sample data -->\n\
    <Loop nrepeat=\"100\">\n\
      <Loop nrepeat=\"100\">\n\
        <ChooseSection nlevel=\"7\">\n\
          <Sample npart=\"1\" mover=\"Free\" species=\"H\"/>\n\
          <Loop nrepeat=\"5\">\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
          </Loop>\n\
        </ChooseSection>\n\
        <Loop nrepeat=\"5\">\n\
          <ChooseSection nlevel=\"5\">\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Loop nrepeat=\"17\">\n\
          <ChooseSection nlevel=\"3\">\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Measure estimator=\"all\"/>\n\
      </Loop>\n\
      <Collect estimator=\"all\"/>\n\
    </Loop>\n\
    <WritePaths file=\"paths.out\"/>\n\
  </PIMC>\n\
</Simulation>";

  std::ofstream file("hatom.xml");
  file << xmlstring << std::endl;
}
