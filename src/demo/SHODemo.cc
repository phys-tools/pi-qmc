#include "SHODemo.h"
#include <fstream>
#include <string>

void SHODemo::generate() const {
  std::cout << 
    "Writing simple harmonic oscillator input file: sho.xml\n" << std::endl;

  std::string xmlstring="\
<?xml version=\"1.0\"?>\n\
<Simulation>\n\
  <SuperCell a=\"25 A\" x=\"1\" y=\"1\" z=\"1\"/>\n\
  <Species name=\"e\" count=\"1\" mass=\"1 m_e\" charge=\"-1\"/>\n\
  <Temperature value=\"1 Ha\" tau=\"0.1 Ha-1\"/>\n\
  <Action>\n\
    <SpringAction/>\n\
    <SHOAction/>\n\
  </Action>\n\
  <Estimators>\n\
    <ThermalEnergyEstimator/>\n\
    <VirialEnergyEstimator nwindow=\"500\"/>\n\
  </Estimators>\n\
  <PIMC>\n\
    <RandomGenerator/>\n\
    <!-- Thermalize -->\n\
    <Loop nrepeat=\"100\">\n\
      <ChooseSection nlevel=\"3\">\n\
         <Sample npart=\"1\" mover=\"Free\" species=\"e\"/>\n\
      </ChooseSection>\n\
    </Loop>\n\
    <!-- Sample data -->\n\
    <Loop nrepeat=\"100\">\n\
      <Loop nrepeat=\"100\">\n\
        <ChooseSection nlevel=\"3\">\n\
          <Sample npart=\"1\" mover=\"Free\" species=\"e\"/>\n\
        </ChooseSection>\n\
        <Measure estimator=\"all\"/>\n\
      </Loop>\n\
      <Collect estimator=\"all\"/>\n\
      <WritePaths file=\"paths.out\"/>\n\
    </Loop>\n\
  </PIMC>\n\
</Simulation>";

  std::ofstream file("sho.xml");
  file << xmlstring << std::endl;
}
