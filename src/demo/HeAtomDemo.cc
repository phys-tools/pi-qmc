#include "HeAtomDemo.h"
#include <fstream>
#include <string>

void HeAtomDemo::generate() const {
  std::cout << "Writing helium atom input file: heatom.xml\n" << std::endl;

  std::string xmlstring="\
<?xml version=\"1.0\"?>\n\
<Simulation>\n\
  <SuperCell a=\"25 A\" x=\"1\" y=\"1\" z=\"1\"/>\n\
  <Species name=\"eup\" count=\"1\" mass=\"1 m_e\" charge=\"-1\"/>\n\
  <Species name=\"edn\" count=\"1\" mass=\"1 m_e\" charge=\"-1\"/>\n\
  <Species name=\"He\" count=\"1\" mass=\"4.0004248 amu\" charge=\"+2\"\n\
           isStatic=\"true\"/>\n\
  <Temperature value=\"2400 K\" nslice=\"12500\"/>\n\
  <Action>\n\
    <SpringAction/>\n\
    <CoulombAction norder=\"2\" rmin=\"0.\" rmax=\"10.\" ngridPoints=\"2000\"\n\
                   dumpFiles=\"true\"/>\n\
    <SphereAction radius=\"5.\" species=\"He\"/>\n\
  </Action>\n\
  <Estimators>\n\
    <ThermalEnergyEstimator/>\n\
    <VirialEnergyEstimator nwindow=\"500\"/>\n\
    <CoulombEnergyEstimator/>\n\
    <PairCFEstimator name=\"rho_up(r)\" species1=\"eup\" species2=\"He\">\n\
      <Radial nbin=\"100\" max=\"5\"/>\n\
    </PairCFEstimator>\n\
    <PairCFEstimator name=\"rho_dn(r)\" species1=\"edn\" species2=\"He\">\n\
      <Radial nbin=\"100\" max=\"5\"/>\n\
    </PairCFEstimator>\n\
  </Estimators>\n\
  <PIMC>\n\
    <RandomGenerator/>\n\
    <SetCubicLattice nx=\"2\" ny=\"2\" nz=\"2\" a=\"1.\"/>\n\
    <!-- Thermalize -->\n\
    <Loop nrepeat=\"10000\">\n\
      <ChooseSection nlevel=\"7\">\n\
         <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
         <Sample npart=\"1\" mover=\"Free\" species=\"edn\"/>\n\
      </ChooseSection>\n\
      <ChooseSection nlevel=\"5\">\n\
         <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
         <Sample npart=\"1\" mover=\"Free\" species=\"edn\"/>\n\
      </ChooseSection>\n\
    </Loop>\n\
    <!-- Sample data -->\n\
    <Loop nrepeat=\"100\">\n\
      <Loop nrepeat=\"50\">\n\
        <Loop nrepeat=\"3\">\n\
          <ChooseSection nlevel=\"9\">\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"eup\" nrepeat=\"5\"/>\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"edn\" nrepeat=\"5\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Loop nrepeat=\"25\">\n\
          <ChooseSection nlevel=\"7\">\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"eup\" nrepeat=\"2\"/>\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"edn\" nrepeat=\"2\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Loop nrepeat=\"71\">\n\
          <ChooseSection nlevel=\"5\">\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"edn\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Loop nrepeat=\"257\">\n\
          <ChooseSection nlevel=\"3\">\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"eup\"/>\n\
            <Sample npart=\"1\" mover=\"Free\" species=\"edn\"/>\n\
          </ChooseSection>\n\
        </Loop>\n\
        <Measure estimator=\"all\"/>\n\
      </Loop>\n\
      <Collect estimator=\"all\"/>\n\
      <WritePaths file=\"paths.out\"/>\n\
    </Loop>\n\
  </PIMC>\n\
</Simulation>";

  std::ofstream file("heatom.xml");
  file << xmlstring << std::endl;
}
