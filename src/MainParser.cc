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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI 
#include <mpi.h>
#endif

extern int irank;

#include "MainParser.h"
#include <iostream>
#include "Action.h"
#include "ActionParser.h"
#include "Algorithm.h"
#include "stats/MPIManager.h"
#include "SimInfoParser.h"
#include "SimulationInfo.h"
#include "PIMCParser.h"
#include "stats/EstimatorManager.h"
#include "EstimatorParser.h"
#include "spin/MainSpinParser.h"
#include <ctime>
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
    irank = MPI::COMM_WORLD.Get_rank();
  }
  //print date
  std::time_t rawtime;
  std::time ( &rawtime );
  if (mpi){ 
    if ( mpi->isMain()) {
      std :: cout << "Start Simulation at current local time and date: "<< std::ctime (&rawtime)<<std ::endl ;
    } 
  }else {
    std :: cout << "Start Simulation at current local time and date: "<< std::ctime (&rawtime)<<std ::endl ;
  }
#else
  std::time_t rawtime;
  std::time(&rawtime);
  std :: cout << "Start Simulation at current local time and date: "<< std::ctime (&rawtime)<<std ::endl ;
#endif
  

  // Find the maximum level for any sampling.
  int maxlevel=1;
  { xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//ChooseSection",ctxt);
    int nsample=obj->nodesetval->nodeNr;
    for (int i=0; i<nsample; ++i) {
      xmlNodePtr& node=obj->nodesetval->nodeTab[i];
      int nlevel=getIntAttribute(node,"nlevel");
      if (nlevel>maxlevel) maxlevel=nlevel;
    }
    // in case you have displace moves and fixed node action. Bug. need to fix for no double action case
        int nslices=1;
    obj = xmlXPathEval(BAD_CAST"//Temperature",ctxt);
    nsample=obj->nodesetval->nodeNr;
    for (int i=0; i<nsample; ++i) {
      xmlNodePtr& node=obj->nodesetval->nodeTab[i];
      nslices = getIntAttribute(node,"nslice");
    }
    obj = xmlXPathEval(BAD_CAST"//SampleDisplaceMove",ctxt);
    nsample=obj->nodesetval->nodeNr;
    for (int i=0; i<nsample; ++i) {
      xmlNodePtr& node=obj->nodesetval->nodeTab[i];
      int nrepeat =getIntAttribute(node,"nrepeat");
      if (mpi){
//	if (nrepeat>0) maxlevel= log(nslices/mpi->getNWorker())/log(2);
      }else{
//	  if (nrepeat>0) maxlevel=log(nslices)/log(2);
      }
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
  // Set the estimators.
  EstimatorParser estimatorParser(simInfo,tau,action,doubleAction,mpi);
  estimatorParser.parse(ctxt);
  EstimatorManager* estimators=estimatorParser.getEstimatorManager();
  estimators->recordInputDocument(filename);
  // Setup a serial PIMC method.
  PIMCParser pimcParser(simInfo,action,doubleAction,estimators,
                        simInfo.getBeadFactory(),mpi);
  pimcParser.parse(ctxt);
  Algorithm* algorithm=pimcParser.getAlgorithm();


  // Run the simulation.
  algorithm->run();

  //print date
#ifdef ENABLE_MPI
  std::time (&rawtime );
  if (mpi){ 
    if ( mpi->isMain()) {
      std :: cout << "\n\n********** Simulation ended successfully :: "<< std::ctime (&rawtime)<<std ::endl ;
    } 
  }else {
    std :: cout << "\n\n********** Simulation ended successfully :: "<< std::ctime (&rawtime)<<std ::endl ;
  }
#else
  std :: cout << "\n\n********** Simulation ended successfully :: "<< std::ctime (&rawtime)<<std ::endl ;
#endif
}
