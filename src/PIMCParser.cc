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
#include "PIMCParser.h"
#include "SerialPaths.h"
#include "ParallelPaths.h"
#include "DoubleParallelPaths.h"
#include "stats/AccRejEstimator.h"
#include "SimulationInfo.h"
#include "FreeMover.h"
#include "spin/SpinMover.h"
#include "spin/FreeSpinMover.h"
#include "DampedFreeTensorMover.h"
#include "FreeTensorMover.h"
#include "Algorithm.h"
#include "Loop.h"
#include "CompositeAlgorithm.h"
#include "SectionChooser.h"
#include "DoubleSectionChooser.h"
#include "HyperbolicMover.h"
#include "Collect.h"
#include "CubicLattice.h"
#include "spin/SpinSetter.h"
#include "SeedRandom.h"
#include "Measure.h"
#include "stats/MPIManager.h"
#include "MultiLevelSampler.h"
#include "DoubleMLSampler.h"
#include "Action.h"
#include "DoubleAction.h"
#include "stats/EstimatorManager.h"
#include "ProbDensityGrid.h"
#include "ConditionalDensityGrid.h"
#include "BinProbDensity.h"
#include "SimpleParticleChooser.h"
#include "SpeciesParticleChooser.h"
#include "PathReader.h"
#include "StructReader.h"
#include "WorkerShifter.h"
#include "WriteProbDensity.h"
#include "WritePaths.h"
#include <iostream>
#include "RandomPermutationChooser.h"
#include "WalkingChooser.h"
#include "NodeTester.h"
#include "FreeParticleNodes.h"

PIMCParser::PIMCParser(const SimulationInfo& simInfo, Action* action,
  DoubleAction* doubleAction, EstimatorManager* estimators,
  const BeadFactory &beadFactory, MPIManager *mpi)
  : XMLUnitParser(simInfo.getUnits()),
    paths(0), algorithm(0), simInfo(simInfo), action(action), 
    doubleAction(doubleAction), estimators(estimators), probDensityGrid(0),
    beadFactory(beadFactory), mpi(mpi) {
}

PIMCParser::~PIMCParser() {
  delete paths;
  delete action;
  delete doubleAction;
  delete probDensityGrid;
  delete algorithm;
}

void PIMCParser::parse(const xmlXPathContextPtr& ctxt) {
  // Read timestep and temperature.
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//PIMC",ctxt);
  xmlNodePtr& pimcNode=obj->nodesetval->nodeTab[0]; ctxt->node=pimcNode;
  bool useDoublePaths=getBoolAttribute(pimcNode,"useDoublePaths");
  xmlXPathFreeObject(obj);
  //bool useDoublePaths=false;
  double t=simInfo.getTemperature();
  double tau=simInfo.getTau();
  const int nslice=(int)(1.0/(tau*t)+0.5);
  // Construct the Paths object.
  if (mpi && (doubleAction || useDoublePaths) ) {
    if (mpi->isMain()) std::cout << "Need double sampling." << std::endl;
    //mpi->setDoubleSampling(true);
  }
  if (mpi && mpi->getNWorker()>1) {
    if (doubleAction) {
      paths=new DoubleParallelPaths(simInfo.getNPart(),nslice,tau,
                                   *simInfo.getSuperCell(),*mpi,beadFactory);
    } else {
      paths=new ParallelPaths(simInfo.getNPart(),nslice,tau,
                             *simInfo.getSuperCell(),*mpi,beadFactory);
    }
  } else {
    paths=new SerialPaths(simInfo.getNPart(),
                          nslice,tau,*simInfo.getSuperCell(),beadFactory);
  }
  // Parse the algorithm.
  algorithm=parseAlgorithm(ctxt);
}

Algorithm* PIMCParser::parseAlgorithm(const xmlXPathContextPtr& ctxt) {
  Algorithm* algorithm(0);
  std::string name((char*)ctxt->node->name);
  //std::cout << name << std::endl;
  if (name=="PIMC"||name=="Composite") {
    CompositeAlgorithm *composite=new CompositeAlgorithm(0);
    parseBody(ctxt,composite);
    algorithm=composite;
  } else if (name=="Loop") {
    int nrepeat=getIntAttribute(ctxt->node,"nrepeat");
    Loop *loop=new Loop(nrepeat);
    parseBody(ctxt,loop);
    algorithm=loop;
  } else if (name=="ChooseSection") {
    nlevel=getIntAttribute(ctxt->node,"nlevel");
    if (doubleAction==0) {
      sectionChooser=new SectionChooser(nlevel,*paths,*action,beadFactory);
      parseBody(ctxt,sectionChooser);
      algorithm=sectionChooser;
    } else {
      doubleSectionChooser=
        new DoubleSectionChooser(nlevel,*paths,*action,
                                 *doubleAction,beadFactory);
      parseBody(ctxt,doubleSectionChooser);
      algorithm=doubleSectionChooser;
    }
  } else if (name=="ShiftWorkers") {
    int maxShift=getIntAttribute(ctxt->node,"maxShift");
    WorkerShifter *shifter=new WorkerShifter(maxShift,*paths,mpi);
    parseBody(ctxt,shifter);
    algorithm=shifter;
  } else if (name=="Measure") {
    std::string estName=getStringAttribute(ctxt->node,"estimator");
    algorithm=new Measure(*paths,estimators->getEstimatorSet(estName));
  } else if (name=="Collect") {
    std::string estName=getStringAttribute(ctxt->node,"estimator");
    algorithm=new Collect(estName,*estimators,getLoopCount(ctxt));
  } else if (name=="RandomGenerator") {
    int iseed=getIntAttribute(ctxt->node,"iseed");
    algorithm=new SeedRandom(iseed);
  } else if (name=="Sample") {
    Mover* mover(0);
    std::string moverName=getStringAttribute(ctxt->node,"mover");
    if (moverName=="Free") mover = new FreeMover(simInfo,nlevel,10.0);
    else if (moverName=="Spin") mover = new SpinMover(simInfo,nlevel,10.0);
    else if (moverName=="FreeSpin") 
      mover = new FreeSpinMover(simInfo,nlevel,10.0);
    else if (moverName=="FreeTensor") mover = new FreeTensorMover(simInfo);
    else if (moverName=="DampedFreeTensor") {
      int saturationLevel=getIntAttribute(ctxt->node,"saturationLevel");
      mover = new DampedFreeTensorMover(simInfo,saturationLevel);
    } else if (moverName=="Hyperbolic") {
      mover = new HyperbolicMover(simInfo,nlevel);
    }   
    int nmoving=getIntAttribute(ctxt->node,"npart");
    std::string speciesName=getStringAttribute(ctxt->node,"species");
    ParticleChooser* particleChooser;
    WalkingChooser* walkingChooser=0;
    bool noPerm=getBoolAttribute(ctxt->node,"noPermutation");
    if (speciesName=="" || speciesName=="all") {
      particleChooser=new SimpleParticleChooser(simInfo.getNPart(),nmoving);
    } else if (noPerm) {
      particleChooser=new SpeciesParticleChooser(
                            simInfo.getSpecies(speciesName),nmoving);
    } else {
      walkingChooser=new WalkingChooser(nmoving,
                          simInfo.getSpecies(speciesName),nlevel,simInfo);
      particleChooser=walkingChooser;
    }
    if (doubleAction==0) {
      PermutationChooser *permutationChooser=0;
      if (walkingChooser) {
        permutationChooser = walkingChooser;
      } else if (noPerm) {
        permutationChooser = new PermutationChooser(nmoving);
      } else {
        permutationChooser = new RandomPermutationChooser(nmoving);
      } 
      int nrepeat=getIntAttribute(ctxt->node,"nrepeat");
      if (nrepeat==0) nrepeat=1;
      algorithm=new MultiLevelSampler(nmoving, *paths, *sectionChooser, 
                     *particleChooser, *permutationChooser,
                     *mover, action, nrepeat, beadFactory);
      if (walkingChooser)
        walkingChooser->setMLSampler((MultiLevelSampler*)algorithm);
    } else {
      bool both=getBoolAttribute(ctxt->node,"both");
      int nrepeat=getIntAttribute(ctxt->node,"nrepeat");
      if (nrepeat==0) nrepeat=1;
      if (walkingChooser) {
        WalkingChooser *walkingChooser2=new WalkingChooser(nmoving,
                          simInfo.getSpecies(speciesName),nlevel,simInfo);
        algorithm=new DoubleMLSampler(nmoving, *paths, *doubleSectionChooser, 
                        *walkingChooser, *walkingChooser,
                        *walkingChooser2, *walkingChooser2, *mover, action,
                        doubleAction, both, nrepeat, beadFactory);
        walkingChooser->setMLSampler((MultiLevelSampler*)algorithm);
        walkingChooser2->setMLSampler((MultiLevelSampler*)algorithm);
      } else {
        PermutationChooser *permutationChooser=0;
        if (noPerm) {
          permutationChooser = new PermutationChooser(nmoving);
        } else {
          permutationChooser = new RandomPermutationChooser(nmoving);
        }
        algorithm=new DoubleMLSampler(nmoving, *paths, *doubleSectionChooser, 
                     *particleChooser, *permutationChooser,
                     *particleChooser, *permutationChooser, *mover, action,
                     doubleAction, both, nrepeat, beadFactory);
      }
    }
    std::string accRejName="MLSampler";
    estimators->add(((MultiLevelSampler*)algorithm)->
                      getAccRejEstimator(accRejName));
  } else if (name=="ProbDensityGrid") {
    double a=getLengthAttribute(ctxt->node,"a");
    IVec n;
    for (int idim=0; idim<NDIM; ++idim) {
      n[idim]=getIntAttribute(ctxt->node,std::string("n")+dimName[idim]);
    }
    probDensityGrid=new ProbDensityGrid(n,a,simInfo,paths),
    algorithm=probDensityGrid;
  } else if (name=="BinProbDensity") {
    algorithm=new BinProbDensity(simInfo,probDensityGrid);
  } else if (name=="WriteProbDensity") {
    std::string filename=getStringAttribute(ctxt->node,"file");
    algorithm=new WriteProbDensity(simInfo,probDensityGrid,filename);
  } else if (name=="ConditionalDensityGrid") {
    double a=getLengthAttribute(ctxt->node,"a");
    std::string speciesName=getStringAttribute(ctxt->node,"species");
    double r=getLengthAttribute(ctxt->node,"radius");
    IVec n;
    Vec center;
    for (int idim=0; idim<NDIM; ++idim) {
      n=getIntAttribute(ctxt->node,std::string("n")+dimName[idim]);
      center[idim]=getLengthAttribute(ctxt->node,
                                      std::string(dimName).substr(idim,1));
    }
    condDensityGrid.push_back(new ConditionalDensityGrid(n,a,simInfo,paths,
      simInfo.getSpecies(speciesName),center,r));
    algorithm=*(--condDensityGrid.end());
  } else if (name=="BinConditionalDensity") {
    int i=getIntAttribute(ctxt->node,"index");
    algorithm=new BinProbDensity(simInfo,condDensityGrid[i]);
  } else if (name=="WriteConditionalDensity") {
    int i=getIntAttribute(ctxt->node,"index");
    std::string filename=getStringAttribute(ctxt->node,"file");
    algorithm=new WriteProbDensity(simInfo,condDensityGrid[i],filename);
  } else if (name=="WritePaths") {
    std::string filename=getStringAttribute(ctxt->node,"file");
    algorithm=new WritePaths(*paths,filename,mpi,beadFactory);
  } else if (name=="SetSpin") {
    algorithm=new SpinSetter(*paths,mpi);
  } else if (name=="SetCubicLattice") {
    double a=getLengthAttribute(ctxt->node,"a");
    Vec aa;
    IVec n;
    double scatter=getDoubleAttribute(ctxt->node,"scatter");
    for (int idim=0; idim<NDIM; ++idim) {
      n[idim]=getIntAttribute(ctxt->node,std::string("n")+dimName[idim]);
      aa[idim]=getLengthAttribute(ctxt->node,
                                  std::string(dimName).substr(idim,1));
    }
    std::string speciesName=getStringAttribute(ctxt->node,"species");
    if (speciesName=="") {
      algorithm=new CubicLattice(*paths,a,scatter,n,aa,mpi);
    } else {
      algorithm=new CubicLattice(*paths,a,scatter,n,aa,
                                 simInfo.getSpecies(speciesName),mpi);
    }
  } else if (name=="ReadPaths") {
    std::string filename=getStringAttribute(ctxt->node,"file");
    int bfactor=getIntAttribute(ctxt->node,"bfactor");
    if (bfactor==0) bfactor=1;
    algorithm=new PathReader(*paths,filename,beadFactory,bfactor,mpi);
  } else if (name=="ReadStruct") {
    std::string filename=getStringAttribute(ctxt->node,"file");
    algorithm=new StructReader(*paths,filename,mpi);
  //} else if (name=="TestNodes") {
  //  std::string filename=getStringAttribute(ctxt->node,"file");
  //  if (filename=="") filename="nodetest.out";
  //  std::string nodeName=getStringAttribute(ctxt->node,"nodeModel");
  //  std::string speciesName=getStringAttribute(ctxt->node,"species");
  //  const Species &species=simInfo.getSpecies(speciesName);
  //  NodeModel* nodes = new FreeParticleNodes(simInfo,
  //                           species,simInfo.getTemperature(),1);
  }
  return algorithm;
}

void PIMCParser::parseBody(const xmlXPathContextPtr& ctxt,
  CompositeAlgorithm* composite) {
  xmlNodePtr node = ctxt->node;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"*",ctxt);
  int nstep=obj->nodesetval->nodeNr;
  composite->resize(nstep);
  for (int i=0; i<nstep; ++i) {
    ctxt->node=obj->nodesetval->nodeTab[i];
    composite->set(i,parseAlgorithm(ctxt));
  }
  xmlXPathFreeObject(obj);
  ctxt->node=node;
}

int PIMCParser::getLoopCount(const xmlXPathContextPtr& ctxt) {
  int count=1;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"ancestor::Loop",ctxt);
  int nloop=obj->nodesetval->nodeNr;
  for (int i=0; i<nloop; ++i) {
    count *= getIntAttribute(obj->nodesetval->nodeTab[i],"nrepeat");
  }
  xmlXPathFreeObject(obj);
  return count;
}
