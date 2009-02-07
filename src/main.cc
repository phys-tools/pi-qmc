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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MainParser.h"
#include "Help.h"
#include <demo/Demo.h>
#include <getopt.h>
void usage();
/** @mainpage PI: Path Integral Monte Carlo Program
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
  * physics.</p>
  * @section goals Immediate Goals
  * <p>There are three immediate goals to keep in mind:</p>
  * <ul>
  * <li>Simulate transport in quantum point contacts (John, Matthew).</li>
  * <li>Create practical all-electron or all-valence-electron calculations
        of atoms and molecules (John).</li>
  * <li>Try spin-state sampling techniques. (Daejin, Hosam) </li>
  * </ul>
  * @section spin Spin sampling and and semiconductor spintronics
  * <p>We should try spin sampling soon. Here is an outline of the approach.</p>
  * <ul>
  * <li>Add spin variables (done, in testing).</li>
  * <li>Add spin sampling (done, in testing).</li>
  * <li>Add magnetization estimator (done, in testing).</li>
  * <li>Add magnetization estimator (done, in testing).</li>
  * </ul>
  * @section plan Implementation Plan
  * <p> The code is now in its second phase of development, where we focus
  *     on production runs. We need the following improvements:</p>
  * <ol>
  * <li>Test and develop new restart features (perhaps some functionality 
        like qmcpack?).</li>
  * <li>Implement built-in coulomb action for 2 and 3 dimensions 
        (partially done).</li>
  * <li>Add long range coulomb action with Ewald sums. (partially done)</li>
  * <li>Improve estimators:
  * <ol>
  * <li>Design and implement general estimator manager.(done)</li>
  * <li>Implement pair corelation functions, structure factors, and 
  *     polarization estimators.(in progress)</li>
  * <li>Collect statistics on accept/reject rates. (done)</li>
  * <li>Develop new nodal models for atoms and molecules.</li>
  * <li>Develop a demo feature and help option (started).</li>
  * <li>Implement superfluid estimator.</li>
  * </ol>
  * </li>
  * <li>Implement position dependent mass.</li>
  * <li>Implement and test spin sampling. (in progress)</li>
  * <li>Implement and test pseudopotential or pseudohamiltonian.</li>
  * </ol>
  * @section outline Outline of code
  * <p>The classes can be viewed in the Class Hierarchy or the Compound
  * Members list.  The program runs by calling the MainParser
  * (a subclass of XMLParser), which
  * constructs and runs a SerialPIMC simulation object, using
  * input from an XML file "pimc.xml".</p>
  * <p>The PIMC algorithm is implemented in SerialPIMC, which also stores
  * the Beads and global permutation.
  * Sampling is implemented in MultiLevelSampler, which has a local copy
  * of moving Beads.  Trial paths are grown using FreeMover for free 
  * particle sampling.</p>
  * @todo Add built-in coulomb interaction.
  * @todo Add local pseudopotentials.
  * @todo Debug spinors.
  * @todo Add help and demo features.
  * @todo Implement atomic and molecular nodes.
  * @author John Shumway */
int main(int argc, char** argv) {
  int rank=0;
#ifdef ENABLE_MPI
  MPI::Init(argc,argv);
  rank=MPI::COMM_WORLD.Get_rank();
#endif
  if (rank==0) {
    std::cout << std::endl;
    std::cout << "PIMC program: "  << PACKAGE_STRING 
#ifdef ENABLE_MPI
              << " (MPI enabled, "
#else
              << " (MPI disabled, "
#endif
#ifdef ENABLE_SPRNG
              << " SPRNG enabled, "
#else
              << " SPRNG disabled, "
#endif
              <<  NDIM << "-d)" << std::endl;
    std::cout << std::endl;
  }
  // Parse command line.
  //bool verbose=false;
  static struct option longopts[] = {
    { "help",     no_argument,  NULL, 'h' },
    { "demo",     optional_argument,  NULL, 'd' },
    { "version",  no_argument,  NULL, 'V' },
    { NULL,       0,            NULL, 0 }
  };
  char ch;
  while (( ch = getopt_long(argc, argv, "hdV", longopts, NULL)) != -1) {
    std::string demoName;
    switch (ch) {
      case 'd': {
        if (optarg!=NULL) {
          demoName=std::string(optarg);
          std::cout << "Requested demo: " << demoName << "\n" << std::endl;
          Demo *demo=Demo::getDemo(demoName);
          if (demo) {
            demo->generate();
            exit(0);
          } else {
            Demo::listDemos(std::cout);
            exit(-1);
          }
        } else {
          Demo::listDemos(std::cout);
          exit(0);
        } }
        break;
      case 'V':
        exit(-1);
      case 'h':
        Help::printHelp();
      default:
        usage();
    }
  }
  argc -= optind;
  argv += optind;

  std::string xmlFileName="pimc.xml";
  if (argc==1) xmlFileName=argv[0];
  MainParser parser(xmlFileName);
  parser.parse();

#ifdef ENABLE_MPI
  MPI::Finalize();
#endif
  return 0;
};

void usage() {
  std::cout << "usage: pi [OPTIONS] [pimc.xml]" << std::endl;
  std::cout << "  OPTIONS" << std::endl;
  std::cout << "     -h, --help           Print a usage message and exit"
            << std::endl;
  std::cout << "     -V, --version        Print version number and exit"
            << std::endl;
  std::cout << "     -d, --demo=name      Write input file for a demo\n";
  std::cout << "                          (omit name to get a list of "
                                      << "available demos)"
            << std::endl;
  exit(0);
}
