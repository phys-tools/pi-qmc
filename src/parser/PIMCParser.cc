// $Id$
/*  Copyright (C) 2004-2011 John B. Shumway, Jr.

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
#include "SuperCell.h"
#include "SerialPaths.h"
#include "ParallelPaths.h"
#include "DoubleParallelPaths.h"
#include "stats/AccRejEstimator.h"
#include "SimulationInfo.h"
#include "sampler/FreeMover.h"
#include "sampler/FreeMoverPBC.h"
#include "sampler/UniformMover.h"
#include "sampler/ExchangeMover.h"
#include "sampler/CollectiveMover.h"
#include "sampler/CollectiveSectionMover.h"
#include "sampler/DisplaceMoveSampler.h"
#include "sampler/DoubleDisplaceMoveSampler.h"
#include "sampler/ModelSampler.h"
#include "sampler/SpinModelSampler.h"
#include "emarate/EMARateMover.h"
#include "EnumeratedModelState.h"
#include "spin/SpinMover.h"
#include "spin/FreeSpinMover.h"
#include "sampler/DampedFreeTensorMover.h"
#include "sampler/FreeTensorMover.h"
#include "Algorithm.h"
#include "Loop.h"
#include "CompositeAlgorithm.h"
#include "sampler/SectionChooser.h"
#include "sampler/DoubleSectionChooser.h"
#include "sampler/MiddleSectionChooser.h"
#include "sampler/NonZeroSectionChooser.h"
#include "sampler/HyperbolicMover.h"
#include "Collect.h"
#include "CubicLattice.h"
#include "RingLattice.h"
#include "spin/SpinSetter.h"
#include "SeedRandom.h"
#include "Measure.h"
#include "stats/MPIManager.h"
#include "sampler/MultiLevelSampler.h"
#include "sampler/DoubleMLSampler.h"
#include "sampler/CollectiveSectionSampler.h"
#include "sampler/DoubleCollectiveSectionSampler.h"
#include "Action.h"
#include "ActionChoice.h"
#include "DoubleAction.h"
#include "stats/EstimatorManager.h"
#include "ProbDensityGrid.h"
#include "ConditionalDensityGrid.h"
#include "BinProbDensity.h"
#include "sampler/SimpleParticleChooser.h"
#include "sampler/SpeciesParticleChooser.h"
#include "sampler/MultiSpeciesParticleChooser.h"
#include "sampler/AlternatingParticleChooser.h"
#include "PathReader.h"
#include "StructReader.h"
#include "WorkerShifter.h"
#include "WriteProbDensity.h"
#include "WritePaths.h"
#include <iostream>
#include "sampler/RandomPermutationChooser.h"
#include "sampler/SpinStatePermutationChooser.h"
#include "sampler/WalkingChooser.h"
#include "sampler/PairChooser.h"
#include "sampler/TwoPairChooser.h"
#include "fixednode/NodeTester.h"
#include "fixednode/FreeParticleNodes.h"


PIMCParser::PIMCParser(const SimulationInfo &simInfo, Action *action,
  DoubleAction *doubleAction, ActionChoiceBase *actionChoice,
  EstimatorManager *estimators,
  const BeadFactory &beadFactory, MPIManager *mpi)
  : XMLUnitParser(simInfo.getUnits()),
    paths(0), algorithm(0), simInfo(simInfo), action(action), 
    doubleAction(doubleAction), actionChoice(actionChoice),
    estimators(estimators), probDensityGrid(0),
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
  double t=simInfo.getTemperature();
  double tau=simInfo.getTau();
  const int nslice=(int)(1.0/(tau*t)+0.5);
  // Construct the Paths object.
  if (mpi && (doubleAction || useDoublePaths) ) {
    if (mpi->isMain()) std::cout << "Need double sampling." << std::endl;
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
  if (actionChoice) {
    paths->setModelState(&actionChoice->getModelState());
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
    //The timer for the main loop is called Main. Sub loops can be named arbitrarily
    std::string timer = getStringAttribute(ctxt->node,"timer");
    int hour=getIntAttribute(ctxt->node,"hours");
    int min=getIntAttribute(ctxt->node,"minutes");
    int sec=getIntAttribute(ctxt->node,"seconds");
    int totalSimTime = hour*3600 + min*60 + sec; 
    if (timer == "Main" && totalSimTime ==0 ) totalSimTime = 12*3600;
    Loop *loop(0);      
    loop=new Loop(nrepeat,totalSimTime,timer,mpi);
    
    parseBody(ctxt,loop);
    algorithm=loop;
  } else if (name=="ChooseMiddleSection") {
    nlevel = getIntAttribute(ctxt->node,"nlevel");
    sectionChooser
      = new MiddleSectionChooser(nlevel,*paths,*action,beadFactory);
    parseBody(ctxt,sectionChooser);
    algorithm = sectionChooser;
  } else if (name=="ChooseNonZeroSection") {
    nlevel = getIntAttribute(ctxt->node,"nlevel");
    sectionChooser
      = new NonZeroSectionChooser(nlevel,*paths,*action,beadFactory);
    parseBody(ctxt,sectionChooser);
    algorithm = sectionChooser;
  } else if (name=="ChooseSection") {
    nlevel=getIntAttribute(ctxt->node,"nlevel");
    if (doubleAction==0) {
std::cout << "doubleAction==0" << std::endl;
      sectionChooser = new SectionChooser(nlevel,paths->getNPart(),*paths,
              *action,beadFactory);
      parseBody(ctxt,sectionChooser);
      algorithm = sectionChooser;
    } else {
std::cout << "doubleAction!=0" << std::endl;
      doubleSectionChooser=
        new DoubleSectionChooser(nlevel, paths->getNPart(), *paths, *action,
                                 *doubleAction,beadFactory);
      parseBody(ctxt,doubleSectionChooser);
      algorithm=doubleSectionChooser;
    }
  } else if (name=="SampleDisplaceMove") {
    std::string moverName=getStringAttribute(ctxt->node, "mover");
    UniformMover* mover(0);
    ParticleChooser* particleChooser=0;
    Vec dist;
    for (int idim=0; idim<NDIM; ++idim) {
      dist[idim] =getLengthAttribute(ctxt->node,std::string("d")+dimName[idim]);
    }
    if (dot(dist,dist)==0)
      dist =getLengthAttribute(ctxt->node,"dist");
    int nmoving=getIntAttribute(ctxt->node,"npart");
    if (nmoving==0) nmoving=1;
    std::string speciesName;
    if (moverName=="Exchange") {
      nmoving = 2;
      int iFirstSlice=0;
      if (doubleAction)
        iFirstSlice=paths->getLowestOwnedSlice(true);
      else
        iFirstSlice=paths->getLowestOwnedSlice(false);
      mover = new ExchangeMover(*paths, iFirstSlice, dist, mpi);
      std::string species1Name, species2Name;
      species1Name = getStringAttribute(ctxt->node,"species1");
      species2Name = getStringAttribute(ctxt->node,"species2");
      particleChooser = new AlternatingParticleChooser(
                             simInfo.getSpecies(species1Name),
                             simInfo.getSpecies(species2Name), nmoving);
      std::cout<<"Picked species1 "<<species1Name<<" and species2 "
               <<species2Name<<" for ExchangeMover."<<std::endl;
    } else if (moverName=="Collective") {
      int iFirstSlice=0;
      if (doubleAction)
        iFirstSlice=paths->getLowestOwnedSlice(true);
      else
        iFirstSlice=paths->getLowestOwnedSlice(false);
      Vec k;
      SuperCell* cell = simInfo.getSuperCell();
      for (int idim=0; idim<NDIM; ++idim) {
        k[idim] = getDoubleAttribute(ctxt->node,std::string("k")+dimName[idim]);
	k[idim] *= (2*3.1415926535897/(cell->a[idim]));
      }
      Vec center;
      for (int idim=0; idim<NDIM; ++idim) {
        center[idim]=getLengthAttribute(ctxt->node,dimName.substr(idim,1));
      }
      IVec random;
      for (int idim=0; idim<NDIM; ++idim) {
        random[idim]=getBoolAttribute(ctxt->node,
                                      std::string("rand")+dimName[idim])?1:0;
      }
      mover = new CollectiveMover(*paths, iFirstSlice, k, dist, center, 
                                  random, mpi);
      speciesName = getStringAttribute(ctxt->node,"species");
   // particleChooser = new SimpleParticleChooser(simInfo.getNPart(),nmoving);
      if (speciesName=="" || speciesName=="all") {
          nmoving = simInfo.getNPart();
          particleChooser = new SimpleParticleChooser(simInfo.getNPart(),nmoving);    
      }else{
      nmoving = simInfo.getSpecies(speciesName).count; 
          particleChooser
            = new SpeciesParticleChooser(simInfo.getSpecies(speciesName),nmoving);
      }

      std::cout<<"Picked "<<nmoving<<" species "<<speciesName<<" particles for CollectiveMover."<<std::endl;
    } else {
      mover = new UniformMover(dist,mpi);
      int nspecies = getIntAttribute(ctxt->node,"nspecies");
      if (nspecies<=1) {                                               
        speciesName=getStringAttribute(ctxt->node,"species");
        std::cout << "Picked species "<< speciesName 
                  << " for displacement." << std::endl;
        if (speciesName=="" || speciesName=="all") {
 	  particleChooser = new SimpleParticleChooser(simInfo.getNPart(),nmoving);
        }else{ 
	  particleChooser 
            = new SpeciesParticleChooser(simInfo.getSpecies(speciesName),nmoving);
        }
      } /*else {
      Species *speciesList = new Species [nspecies];
      for (int ispec=0; ispec<nspecies; ispec++){
	std::stringstream sispec;
	sispec << "species"<<(ispec+1);
	std::string speciesName=getStringAttribute(ctxt->node,sispec.str());
	std::cout << "Picked species " << speciesName << " for displacement."
                  << std::endl;
	speciesList[ispec]=simInfo.getSpecies(speciesName);
      }
      particleChooser = new MultiSpeciesParticleChooser(
                              speciesList, nspecies, nmoving);
      delete [] speciesList;
      }*/
    }
    bool delayedRejection=getBoolAttribute(ctxt->node,"delayedRejection"); 
    double newDistFactor=getDoubleAttribute(ctxt->node,"distFactor");
    if (newDistFactor==0 && delayedRejection) newDistFactor=0.5;
    
    int nrepeat=getIntAttribute(ctxt->node,"nrepeat");
    if (doubleAction) {
      algorithm = new DoubleDisplaceMoveSampler(nmoving, nrepeat, 
						*paths, *particleChooser, 
						*mover, action, doubleAction, mpi);
      std :: cout <<"Using DoubleDisplaceMoveSampler for species: "<<speciesName<<std :: endl;
    } else {    
      algorithm = new DisplaceMoveSampler(
                        nmoving, nrepeat, *paths, *particleChooser, 
                        *mover, action, mpi);
      std :: cout <<"Using DisplaceMoveSampler."<<std :: endl;
    } 
    std::string accRejName="DisplaceMoveSampler";
    estimators->add(((DisplaceMoveSampler*)algorithm)->
		    getAccRejEstimator(accRejName));

  } else if (name=="SampleModel") {
    int target = getIntAttribute(ctxt->node,"target");
    if (target==0) target=-1;
    algorithm = new ModelSampler(*paths, action, actionChoice, target, mpi);
    estimators->add(((ModelSampler*)algorithm)->
		    getAccRejEstimator("ModelSampler"));
  } else if (name=="SampleSpinModel") {
    algorithm = new SpinModelSampler(*paths, action, actionChoice, mpi);
    estimators->add(((SpinModelSampler*)algorithm)->
		    getAccRejEstimator("ModelSampler"));
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
    int iseed = getIntAttribute(ctxt->node,"iseed");
    algorithm=new SeedRandom(iseed);
  } else if (name=="Sample") {
    // delayed rejection stuff
    bool delayedRejection=getBoolAttribute(ctxt->node,"delayedRejection");  
    double defaultFactor=getDoubleAttribute(ctxt->node,"defaultSigmaFactor");
    if (defaultFactor==0) defaultFactor=1;
    double newFactor=getDoubleAttribute(ctxt->node,"newSigmaFactor");
    if (delayedRejection!=0 && newFactor==0) newFactor=0.75;
    Mover* mover(0);
    ParticleChooser* particleChooser=0;
    ParticleChooser* particleChooser2=0;
    PermutationChooser* permutationChooser=0;
    PermutationChooser* permutationChooser2=0;
    int nmoving=0;
    bool noPerm=true;
    std::string moverName=getStringAttribute(ctxt->node,"mover");
    if (moverName=="Free" || moverName=="")
      mover = new FreeMover(simInfo,nlevel,10.0);
    if (moverName=="FreePBC" || moverName=="")
      mover = new FreeMoverPBC(simInfo,nlevel,10.0);
    else if (moverName=="Spin") mover = new SpinMover(simInfo,nlevel,10.0);
    else if (moverName=="FreeSpin") 
      mover = new FreeSpinMover(simInfo,nlevel,10.0);
    else if (moverName=="FreeTensor") mover = new FreeTensorMover(simInfo);
    else if (moverName=="DampedFreeTensor") {
      int saturationLevel=getIntAttribute(ctxt->node,"saturationLevel");
      mover = new DampedFreeTensorMover(simInfo,saturationLevel);
    } else if (moverName=="Hyperbolic") {
      mover = new HyperbolicMover(simInfo,nlevel);
    } else if (moverName=="EMARate") {
      double c=getDoubleAttribute(ctxt->node,"c");
      std::string speciesName1=getStringAttribute(ctxt->node,"species1");
      std::string speciesName2=getStringAttribute(ctxt->node,"species2");
      EMARateMover* emaRateMover = new EMARateMover(simInfo,nlevel,c);
      mover = emaRateMover;
      particleChooser = emaRateMover;
      permutationChooser = emaRateMover;
      nmoving = 2;
    }   
    if (!particleChooser) {
    std::string chooserName=getStringAttribute(ctxt->node,"chooser");
    nmoving=getIntAttribute(ctxt->node,"npart");
    std::string speciesName=getStringAttribute(ctxt->node,"species");
    bool noPerm=getBoolAttribute(ctxt->node,"noPermutation");
    if (speciesName=="" || speciesName=="all") {
      particleChooser=new SimpleParticleChooser(simInfo.getNPart(),nmoving);
      if (doubleAction) {
        particleChooser2=new SimpleParticleChooser(simInfo.getNPart(),nmoving);
      }
    } else if (noPerm) {
      particleChooser=new SpeciesParticleChooser(
                            simInfo.getSpecies(speciesName),nmoving);
      if (doubleAction) {
        particleChooser2=new SpeciesParticleChooser(
                               simInfo.getSpecies(speciesName),nmoving);
      }
    } else if (chooserName=="twoPair") {
      std::string species2Name=getStringAttribute(ctxt->node,"species2");
      if (species2Name=="") species2Name=speciesName;
      TwoPairChooser* chooser
        = new TwoPairChooser(simInfo.getSpecies(speciesName),
                             simInfo.getSpecies(species2Name));
      particleChooser = chooser;
      permutationChooser = chooser;
      if (doubleAction) {
        TwoPairChooser* chooser
          = new TwoPairChooser(simInfo.getSpecies(speciesName),
                               simInfo.getSpecies(species2Name));
        particleChooser2 = chooser;
        permutationChooser2 = chooser;
      }
    } else if (chooserName=="pair") {
      std::string species2Name=getStringAttribute(ctxt->node,"species2");
      PairChooser* chooser
        = new PairChooser(nmoving,simInfo.getSpecies(speciesName),
                          simInfo.getSpecies(species2Name),nlevel,simInfo);
      particleChooser = chooser;
      permutationChooser = chooser;
      if (doubleAction) {
        PairChooser* chooser
          = new PairChooser(nmoving,simInfo.getSpecies(speciesName),
                            simInfo.getSpecies(species2Name),nlevel,simInfo);
        particleChooser2 = chooser;
        permutationChooser2 = chooser;
      }
      nmoving*=2;
    } else {
      if (actionChoice && (paths->getModelState()->isSpinModelState())) {
	  SpinStatePermutationChooser* chooser 
	    = new SpinStatePermutationChooser(nmoving,
		simInfo.getSpecies(speciesName),nlevel,simInfo,
		actionChoice->getModelState());
          particleChooser = chooser;
          permutationChooser = chooser;
      } else {
        WalkingChooser* chooser
          = new WalkingChooser(nmoving,simInfo.getSpecies(speciesName),
                               nlevel,simInfo);
        particleChooser = chooser;
        permutationChooser = chooser;
      }
      if (doubleAction) {
        if (actionChoice && (paths->getModelState()->isSpinModelState())) {
	  SpinStatePermutationChooser* chooser 
	    = new SpinStatePermutationChooser(nmoving,
		simInfo.getSpecies(speciesName),nlevel,simInfo,
		actionChoice->getModelState());
          particleChooser2 = chooser;
          permutationChooser2 = chooser;
	} else {
          WalkingChooser* chooser
            = new WalkingChooser(nmoving,simInfo.getSpecies(speciesName),
                                 nlevel,simInfo);
          particleChooser2 = chooser;
          permutationChooser2 = chooser;
	}
      }
    }
    }
    if (doubleAction==0) {
      if (!permutationChooser) {
        if (noPerm) {
          permutationChooser = new PermutationChooser(nmoving);
        } else {
          permutationChooser = new RandomPermutationChooser(nmoving);
        }
      } 
      int nrepeat=getIntAttribute(ctxt->node,"nrepeat");
      if (nrepeat==0) nrepeat=1;
      algorithm=new MultiLevelSampler(nmoving, *paths, *sectionChooser, 
                     *particleChooser, *permutationChooser,
                     *mover, action, nrepeat, beadFactory, delayedRejection,
				      defaultFactor,newFactor);
      permutationChooser->setMLSampler((MultiLevelSampler*)algorithm);
    } else {
      bool both=getBoolAttribute(ctxt->node,"both");
      int nrepeat=getIntAttribute(ctxt->node,"nrepeat");
      if (nrepeat==0) nrepeat=1;
      if (!permutationChooser) {
        if (noPerm) {
          permutationChooser = new PermutationChooser(nmoving);
          permutationChooser2 = permutationChooser;
        } else {
          permutationChooser = new RandomPermutationChooser(nmoving);
          permutationChooser2 = permutationChooser;
        }
      }
      algorithm = new DoubleMLSampler(nmoving, *paths, *doubleSectionChooser, 
				      *particleChooser, *permutationChooser,
				      *particleChooser2, *permutationChooser2, *mover, action,
				      doubleAction, both, nrepeat, beadFactory, delayedRejection,
				      defaultFactor,newFactor);
      permutationChooser->setMLSampler((MultiLevelSampler*)algorithm);
      permutationChooser2->setMLSampler((MultiLevelSampler*)algorithm);
    }
    std::string accRejName="MLSampler";
    estimators->add(((MultiLevelSampler*)algorithm)->
                      getAccRejEstimator(accRejName));
  } else if (name=="SampleCollective") {
    int nrepeat=getIntAttribute(ctxt->node,"nrepeat");
    if (nrepeat==0) nrepeat=1;
    int npart=simInfo.getNPart();
    double radius=getLengthAttribute(ctxt->node,"radius");
    SuperCell *cell=simInfo.getSuperCell();
    Vec amplitude, min, max;
    for (int idim=0; idim<NDIM; ++idim) {
      amplitude[idim]=getLengthAttribute(ctxt->node,
	                                 std::string("d")+dimName[idim]);
      min[idim]=-cell->a[idim]/2.;
      max[idim]=cell->a[idim]/2.;
    }
    parseBoundary(ctxt,min,max);
    CollectiveSectionMover *mover=new CollectiveSectionMover(radius,amplitude,
	                                                   npart,min,max,cell);
    bool both=getBoolAttribute(ctxt->node,"both");
    if (doubleAction==0) {
      algorithm = new CollectiveSectionSampler(npart, *sectionChooser, action,
	                                          nrepeat, beadFactory, mover);
    } else {
      algorithm = new DoubleCollectiveSectionSampler(npart, 
	                     *doubleSectionChooser, action, doubleAction,
			     nrepeat, beadFactory, mover, both);
    }
    std::string accRejName="CollectiveSampler";
    estimators->add(((CollectiveSectionSampler*)algorithm)->
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
      n[idim]=getIntAttribute(ctxt->node,std::string("n")+dimName[idim]);
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
    int dumpFreq=getIntAttribute(ctxt->node,"freq");
    int maxConfigs=getIntAttribute(ctxt->node,"configs");
    maxConfigs = (maxConfigs==0)?500:maxConfigs;
    bool writeMovie=0;
    if (mpi){
      if (mpi->isMain()) writeMovie=getBoolAttribute(ctxt->node,"movie");     
    } else {
      writeMovie=getBoolAttribute(ctxt->node,"movie");
    }
    algorithm=new WritePaths(*paths,filename,dumpFreq,maxConfigs,writeMovie,simInfo,mpi,beadFactory);
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
      if (n[idim]==0) n[idim]=1;
    }
    std::string speciesName=getStringAttribute(ctxt->node,"species");
    if (speciesName=="") {
      algorithm=new CubicLattice(*paths,a,scatter,n,aa,mpi);
    } else {
      algorithm=new CubicLattice(*paths,a,scatter,n,aa,
                                 simInfo.getSpecies(speciesName),mpi);
    }
  } else if (name=="SetRingLattice") {
    double radius = getLengthAttribute(ctxt->node,"radius");
    double angle0 = getDoubleAttribute(ctxt->node,"angle0");
    double anglef = getDoubleAttribute(ctxt->node,"anglef");
    double anglex = getDoubleAttribute(ctxt->node,"anglex");
    if (anglef == angle0 && angle0 != 0) {
      std::cout<<"The ending angle cannot be the same as the beginning angle."<<std::endl;
      exit(-1);
    }
    std::string speciesName=getStringAttribute(ctxt->node,"species");
    int numpart = paths->getNPart();
    if (anglef == 0 && angle0 == 0) {
      if (speciesName=="") {
      } else {
        Species temp = simInfo.getSpecies(speciesName);
        numpart = temp.count;
      }
      anglef = 2 * 3.141592653589793 * (numpart - 1) / numpart;
    }
    std::cout<<"angle0 = "<<angle0<<"; anglef = "<<anglef<<std::endl;
    if (speciesName=="") {
      algorithm = new RingLattice(*paths,radius,angle0,anglef,anglex,mpi);
    } else {
      algorithm = new RingLattice(*paths,radius,angle0,anglef,anglex,simInfo.getSpecies(speciesName),mpi);
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

void PIMCParser::parseBoundary(const xmlXPathContextPtr& ctxt, 
			       Vec& min, Vec& max) {
//  ctxt->node = estNode;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"*",ctxt);
  int N=obj->nodesetval->nodeNr;
  int idir=0;
  for (int idist=0; idist<N; ++idist) {
    xmlNodePtr distNode=obj->nodesetval->nodeTab[idist];
    std::string name=getName(distNode);
    if (name=="Cartesian") {
      std::string dirName = getStringAttribute(distNode,"dir");
      for (int i=0; i<NDIM; ++i) {
	if (dirName==dimName.substr(i,1)) {
	  min[i] = getLengthAttribute(distNode,"min");
	  max[i] = getLengthAttribute(distNode,"max");
	}
      }
    }
  }
}
