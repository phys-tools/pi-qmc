#include "parser/MainParser.h"
#include "util/startup/CommandLineParser.h"
#include "util/startup/MPILifecycle.h"
#include "util/startup/Greeting.h"
#include <string>

/** @mainpage pi: Path integral quantum Monte Carlo
 * @section usage Usage
 * To run the program, supply a pimc.xml file and any other input files,
 * then run the executable.
 * If you don't have an input file, try
 * @code
 * pi --demo FreeParticles
 * @endcode
 * This will generate a pimc.xml file for free particles.
 * Use "--demo" without any name to see a list of other demos.
 * @section overview Overview
 * <p>This program implements the path integral Monte Carlo technique.
 * The goal of this program is to develop PIMC methods for applied
 * physics.
 * This code is hosted at http://code.google.com/p/pi-qmc ,
 * where you can also find some documentation and user support. </p>
 * @section structure Code structure
 * <p>To get started understanding the code, take a look at the
 * <a href="hierarchy.html">Class Hierarchy</a>.
 * The objects needed for the simulation are constructed in the
 * XMLParser classes: SimInfoParser, ActionParser, EstimatorParser
 * and PIMCParser. The main Monte Carlo algorithm is implemented in
 * MultiLevelSampler, and path coordinates are stored in Paths and Beads.
 * @author John Shumway */
int main(int argc, char** argv) {
    MPILifecycle::initialize(argc, argv);
    Greeting::print();

    std::string xmlFileName = CommandLineParser::parse(argc, argv);

    MainParser parser(xmlFileName);
    parser.parse();

    MPILifecycle::finalize();
    return 0;
}

