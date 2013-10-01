#include "config.h"
#ifdef ENABLE_MPI 
#include <mpi.h>
#endif
#include "action/Action.h"
#include "algorithm/Algorithm.h"
#include "base/SimulationInfo.h"
#include "parser/MainParser.h"
#include "parser/ActionParser.h"
#include "parser/SimInfoParser.h"
#include "parser/PIMCParser.h"
#include "parser/EstimatorParser.h"
#include "spin/MainSpinParser.h"
#include "stats/EstimatorManager.h"
#include "stats/MPIManager.h"
#include "util/shiny/Shiny.h"
#include <iostream>
#include <ctime>

class ActionChoiceBase;

MainParser::MainParser(const std::string& filename) 
  : filename(filename), context (0) {
  doc = xmlParseFile((char*)filename.c_str());
  context = xmlXPathNewContext(doc);
}

MainParser::~MainParser() {
  xmlXPathFreeContext(context);
  xmlFreeDoc(doc);
}

void MainParser::parse() {
  parse(context);
}

void MainParser::parse(const xmlXPathContextPtr& ctxt) {
  PROFILE_BEGIN(Parser);

  MPIManager *mpi=0;
#ifdef ENABLE_MPI
  { xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//PIMC",ctxt);
    xmlNodePtr& pimcNode=obj->nodesetval->nodeTab[0]; ctxt->node=pimcNode;
    int nclone=getIntAttribute(pimcNode,"nclone");
    int nworker=getIntAttribute(pimcNode,"nworker");
    xmlXPathFreeObject(obj);
    if (nworker==0) nworker=1;
    if (nclone==0) nclone=MPI::COMM_WORLD.Get_size()/nworker;
    mpi=new MPIManager(nworker,nclone);
  }
  //print date
  std::time_t rawtime;
  std::time ( &rawtime );
  if (mpi){ 
    if ( mpi->isMain()) {
      std::cout << 
        "Start Simulation at current local time and date: " <<
        std::ctime (&rawtime) <<std ::endl ;
    } 
  }else {
    std::cout << 
      "Start Simulation at current local time and date: " <<
      std::ctime (&rawtime)<<std ::endl ;
  }
#else
  std::time_t rawtime;
  std::time(&rawtime);
  std :: cout << "Start Simulation at current local time and date: "<< std::ctime (&rawtime)<<std ::endl ;
#endif
  

  // Find the maximum level for any sampling.
  int maxlevel=1;
  { // Find all tags that start with Choose and end with Section. 
    xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//*[starts-with(name(),'Choose') and (substring(name(),string-length(name())-6) = 'Section')]",ctxt);
    //xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//*[starts-with(name(),'Choose') and ends-with(name(),'Section')]",ctxt);
    int nsample=obj->nodesetval->nodeNr;
    for (int i=0; i<nsample; ++i) {
      xmlNodePtr& node=obj->nodesetval->nodeTab[i];
      int nlevel=getIntAttribute(node,"nlevel");
std::cout << "Level = " << nlevel << std::endl;
      if (nlevel>maxlevel) maxlevel=nlevel;
    }
  }
  std::cout << "  maxlevel = " << maxlevel << std::endl;
  // Get the simulation info.
  SimInfoParser simInfoParser;
  simInfoParser.parse(ctxt);
  SimulationInfo &simInfo=simInfoParser.getSimInfo();
  if (!mpi || mpi->isMain()) std::cout << simInfo << std::endl;

  // Check if we need spin.
  MainSpinParser spinParser(simInfo);
  spinParser.parse(ctxt);

  // Set the action.
  const double tau=simInfo.getTau();
  ActionParser actionParser(simInfo, maxlevel,mpi);
  actionParser.parse(ctxt);
  Action* action=actionParser.getAction();
  DoubleAction* doubleAction=actionParser.getDoubleAction();
  ActionChoiceBase* actionChoice = actionParser.getActionChoice();
  // Set the estimators.
  EstimatorParser estimatorParser(simInfo,tau,action,doubleAction,
    actionChoice,mpi);
  estimatorParser.parse(ctxt);
  EstimatorManager* estimators=estimatorParser.getEstimatorManager();
  estimators->setInputFilename(filename);
//  estimators->recordInputDocument(filename);
  // Setup the PIMC method.
  PIMCParser pimcParser(simInfo, action, doubleAction, actionChoice,
      estimators, simInfo.getBeadFactory(),mpi);
  pimcParser.parse(ctxt);
  Algorithm* algorithm=pimcParser.getAlgorithm();
  PROFILE_END();


  // Run the simulation.
  PROFILE_BEGIN(Run);
  algorithm->run();
  PROFILE_END();

  //print date
#ifdef ENABLE_MPI
  std::time (&rawtime );
  if (mpi){ 
    if (mpi->isMain()) {
      std::cout << 
          "\n\n********** Simulation ended successfully :: " <<
          std::ctime (&rawtime)<<std ::endl ;
      std::cout << "\n\n** Timing Profile Data **\n" << std::endl;
      PROFILE_UPDATE(); // uptate all timing profiles
      PROFILE_OUTPUT("Shiny.out"); // print profile to file
      PROFILE_OUTPUT(NULL); // print profile to terminal
    } 
  } else {
    std::cout << 
        "\n\n********** Simulation ended successfully :: " <<
        std::ctime (&rawtime)<<std ::endl ;
    std::cout << "\n\n** Timing Profile Data **\n" << std::endl;
    PROFILE_UPDATE(); // uptate all timing profiles
    PROFILE_OUTPUT("Shiny.out"); // print profile to file
    PROFILE_OUTPUT(NULL); // print profile to terminal
  }
#else
  std::cout << 
      "\n\n********** Simulation ended successfully :: " <<
      std::ctime(&rawtime) <<std::endl;
  std::cout << "\n\n** Timing Profile Data **\n" << std::endl;
  PROFILE_UPDATE(); // uptate all timing profiles
  PROFILE_OUTPUT("Shiny.out"); // print profile to file
  PROFILE_OUTPUT(NULL); // print profile to terminal
#endif
}
