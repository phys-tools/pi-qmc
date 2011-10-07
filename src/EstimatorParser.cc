// $Id$
/*  Copyright (C) 2004-2009 John B. Shumway, Jr.

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
#include "EstimatorParser.h"
#include "stats/EstimatorManager.h"
#include "SimulationInfo.h"
#include "Action.h"
#include "ActionChoice.h"
#include "AngularMomentumEstimator.h"
#include "CoulombAction.h"
#include "ConductivityEstimator.h"
#include "ConductivityEstimator2D.h"
#include "ConductanceEstimator.h"
#include "DensDensEstimator.h"
#include "DensityEstimator.h"
#include "DensCountEstimator.h"
#include "DiamagneticEstimator.h"
#include "DynamicPCFEstimator.h"
#include "CountCountEstimator.h"
#include "Distance.h"
#include "EMARateEstimator.h"
#include "FreeEnergyEstimator.h"
#include "PairDistance.h"
#include "CoulombEnergyEstimator.h"
#include "EwaldCoulombEstimator.h"
#include "DoubleAction.h"
#include "ThermoEnergyEstimator.h"
#include "VirialEnergyEstimator.h"
#include "DipoleMomentEstimator.h"
#include "BondLengthEstimator.h"
#include "FrequencyEstimator.h"
#include "spin/SpinEstimator.h"
#include "PositionEstimator.h"
#include "BoxEstimator.h"
#include "PairCFEstimator.h"
#include "ZeroVarDensityEstimator.h"
#include "PermutationEstimator.h"
#include "JEstimator.h"
#include "stats/MPIManager.h"
#include "SuperCell.h"
#include "SpinChargeEstimator.h"
#include "SHOPhase.h"
#include "VIndEstimator.h"
#include "EIndEstimator.h"
#include "SimInfoWriter.h"
#include "SKOmegaEstimator.h"
#include "WindingEstimator.h"
#include "stats/Units.h"
#include "ModelState.h"

EstimatorParser::EstimatorParser(const SimulationInfo& simInfo,
    const double tau, const Action* action, const DoubleAction* doubleAction,
    const ActionChoiceBase *actionChoice, MPIManager *mpi)
  : XMLUnitParser(simInfo.getUnits()), manager(0),
    simInfo(simInfo), tau(tau), action(action), doubleAction(doubleAction),
    actionChoice(actionChoice), mpi(mpi) {
  manager=new EstimatorManager("pimc.h5", mpi, new SimInfoWriter(simInfo));
}

EstimatorParser::~EstimatorParser() {delete manager;}

void EstimatorParser::parse(const xmlXPathContextPtr& ctxt) {
  // First see if we need a FreeEnergyEstimator for ActionChoice.
  if (actionChoice) {
    manager->add(new FreeEnergyEstimator(simInfo,
       actionChoice->getModelState().getModelCount(), mpi));
  }
  // Then parse the xml estimator list.
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Estimators/*",ctxt);
  int nest=obj->nodesetval->nodeNr;
  for (int iest=0; iest<nest; ++iest) {
    xmlNodePtr estNode=obj->nodesetval->nodeTab[iest];
    std::string name=getName(estNode);
    if (name=="ThermalEnergyEstimator") {
      std::string unitName=getStringAttribute(estNode,"unit");
      double shift=getEnergyAttribute(estNode,"shift");
      double perN=getDoubleAttribute(estNode,"perN");
      double scale = 1.;
      if (perN>0) scale/=perN;
      if (unitName!="") scale/=simInfo.getUnits()->getEnergyScaleIn(unitName);
      manager->add(new ThermoEnergyEstimator(simInfo,action,doubleAction,mpi,
                                             unitName,scale,shift));
    }
    if (name=="SpinEstimator") {
      double gc=getDoubleAttribute(estNode,"gc");
      manager->add(new SpinEstimator(simInfo,0,gc));
      manager->add(new SpinEstimator(simInfo,1,gc));
      manager->add(new SpinEstimator(simInfo,2,gc));
    }
    if (name=="VirialEnergyEstimator") {
      int nwindow=getIntAttribute(estNode,"nwindow");
      if (nwindow==0) nwindow=1;
      std::string unitName=getStringAttribute(estNode,"unit");
      double shift=getEnergyAttribute(estNode,"shift");
      int perN=getIntAttribute(estNode,"perN");
      double scale = 1.;
      if (perN>0) scale/=perN;
      if (unitName!="") scale/=simInfo.getUnits()->getEnergyScaleIn(unitName);
      manager->add(
        new VirialEnergyEstimator(simInfo,action,doubleAction,nwindow,mpi,
                                  unitName,scale,shift));
    }
    if (name=="CoulombEnergyEstimator") {
      double epsilon=getDoubleAttribute(estNode,"epsilon");
      if (epsilon==0) epsilon=1;
      std::string unitName=getStringAttribute(estNode,"unit");
      double shift=getEnergyAttribute(estNode,"shift");
      int perN=getIntAttribute(estNode,"perN");
      double scale = 1.;
      if (perN>0) scale/=perN;
      if (unitName!="") scale/=simInfo.getUnits()->getEnergyScaleIn(unitName);
      bool useEwald=getBoolAttribute(estNode,"useEwald");
      if (useEwald) {
        double rcut=getLengthAttribute(estNode,"rcut");
        double kcut=getDoubleAttribute(estNode,"kcut");
	std::string ewaldType=getStringAttribute(estNode,"ewaldType");
	if (ewaldType=="") ewaldType="optEwald";
	if (ewaldType=="opt") ewaldType="optEwald";
	if (ewaldType=="trad") ewaldType="tradEwald";
	
	if (ewaldType=="" || ewaldType=="optEwald"){
	  manager->add(new EwaldCoulombEstimator(simInfo,action,
						 epsilon,rcut,kcut,mpi,unitName,scale,shift));
	} else 
	  if ( ewaldType=="tradEwald"){
	    int nimages=getIntAttribute(estNode, "ewaldImages");
	    if (nimages ==0) nimages=1;
	    double kappa=getDoubleAttribute(estNode,"kappa");
	    if (kappa==0) {
	      double Lside=1000000000000000.0 ;
	      for (int i=0; i<NDIM; i++) {
		if (Lside > (*simInfo.getSuperCell()).a[i] )
		  Lside = (*simInfo.getSuperCell()).a[i];
	      }
	      kappa = sqrt( pow(3.6*simInfo.getNPart(),(1.0/6.0))*sqrt(3.1415926535897931)/(Lside) );
	    }
	    std :: cout <<"Estimator Parser :: CoulombEnergyEstimator :: kappa :: "<< kappa<<std :: endl;
	    bool testEwald=getBoolAttribute(estNode,"testEwald");

	    manager->add(new EwaldCoulombEstimator(simInfo,action,
						   epsilon,rcut,kcut,mpi,unitName,scale,shift, kappa, nimages,testEwald));
	  }
      } else {
        manager->add(new CoulombEnergyEstimator(simInfo,action,epsilon,mpi,
                                                unitName,scale,shift));
      }
    }
    if (name=="AngularMomentumEstimator") {
      double omega=getEnergyAttribute(estNode,"omega");
      double b=getEnergyAttribute(estNode,"b");
      PhaseModel* phaseModel=new SHOPhase(simInfo,simInfo.getSpecies(0),omega,
                                          simInfo.getTemperature(),b,0);
      manager->add(new AngularMomentumEstimator(simInfo,phaseModel));
    }
    if (name=="DipoleMomentEstimator") {
      std::string component=getStringAttribute(estNode,"component");
      int index=2; //default use z component
      if (component=="x") index=0;
      else if (component=="y") index=1;
      manager->add(new DipoleMomentEstimator(simInfo, index));
    }
    if (name=="BondLengthEstimator") {
      std::string species1=getStringAttribute(estNode,"species1");
      std::string species2=getStringAttribute(estNode,"species2");
      std::string unitName=getStringAttribute(estNode,"unit");
      manager->add(new BondLengthEstimator(simInfo,simInfo.getSpecies(species1)
			      ,simInfo.getSpecies(species2),unitName));
    }
    if (name=="FrequencyEstimator") {
      std::string species1=getStringAttribute(estNode,"species1");
      std::string species2=getStringAttribute(estNode,"species2");
      int nfreq=getIntAttribute(estNode,"nfreq");
      if (nfreq==0) nfreq=simInfo.getNSlice();
      int nstride=getIntAttribute(estNode,"nstride");
      if (nstride==0) nstride=1;
      manager->add(new FrequencyEstimator(simInfo,simInfo.getSpecies(species1),
	simInfo.getSpecies(species2), nfreq, nstride, mpi));
    }
    if (name=="BoxEstimator") {
      std::string component=getStringAttribute(estNode,"component");
      int index=2; //default use z component
      if (component=="x") index=0;
      else if (component=="y") index=1;
      double lower=getDoubleAttribute(estNode,"lower");
      double upper=getDoubleAttribute(estNode,"upper");
      std::string species=getStringAttribute(estNode,"species");
      manager->add(new BoxEstimator(simInfo, simInfo.getSpecies(species), 
                                    lower, upper, index));
    }
    if (name=="PositionEstimator") {
      std::string component=getStringAttribute(estNode,"component");
      int index=2; //default use z component
      if (component=="x") index=0;
      else if (component=="y") index=1;
      std::string species=getStringAttribute(estNode,"species");
      manager->add(new PositionEstimator(simInfo, 
                                         simInfo.getSpecies(species), index));
    }
    if (name=="ConductivityEstimator") {
      int nbin=getIntAttribute(estNode,"nbin");
      if (nbin==0) nbin=1;
      int ndbin=getIntAttribute(estNode,"ndbin");
      if (ndbin==0) ndbin=1;
      int nfreq=getIntAttribute(estNode,"nfreq");
      if (nfreq==0) nfreq=simInfo.getNSlice();
      int nstride=getIntAttribute(estNode,"nstride");
      if (nstride==0) nstride=1;
      manager->add(
        new ConductivityEstimator(simInfo,nfreq,nbin,ndbin,nstride,mpi));
      if (getBoolAttribute(estNode,"calcInduced")) {
        std::cout << "Also calculating induced voltage" << std::endl;
        const Action *coulAction=action->getAction(typeid(CoulombAction));
        if (!coulAction) {
          std::cout << "ERROR, can't find CoulombAction" << std::endl;
        } else {
          //manager->add(new VIndEstimator(simInfo,
          //    dynamic_cast<const CoulombAction*>(coulAction),
          //    nfreq,nbin,ndbin,nstride,mpi));
          manager->add(new VIndEstimator(simInfo,
              dynamic_cast<const CoulombAction*>(coulAction),
              nfreq,nbin,nstride,mpi));
        }
      }
      if (getBoolAttribute(estNode,"calcEInduced")) {
        std::cout << "Also calculating induced efield" << std::endl;
        const Action *coulAction=action->getAction(typeid(CoulombAction));
        if (!coulAction) {
          std::cout << "ERROR, can't find CoulombAction" << std::endl;
        } else {
          manager->add(new EIndEstimator(simInfo,
              dynamic_cast<const CoulombAction*>(coulAction),
              nfreq,nbin,ndbin,nstride,mpi));
        }
      }
    }
    if (name=="ConductivityEstimator2D") {
      std::string dim1 = getStringAttribute(estNode, "dim1");
      if(dim1.length()==0) dim1="xy";
      std::string dim2 = getStringAttribute(estNode, "dim2");
      if(dim2.length()==0) dim2="xy";
      int nxbin=getIntAttribute(estNode,"nxbin");
      if (nxbin==0) nxbin=1;
      int nybin=getIntAttribute(estNode,"nybin");
      if (nybin==0) nybin=1;
      int nxdbin=getIntAttribute(estNode,"nxdbin");
      if (nxdbin==0) nxdbin=1;
      int nydbin=getIntAttribute(estNode,"nydbin");
      if (nydbin==0) nydbin=1;
      double xmin=getLengthAttribute(estNode,"xmin");
      double xmax=getLengthAttribute(estNode,"xmax");
      double ymin=getLengthAttribute(estNode,"ymin");
      double ymax=getLengthAttribute(estNode,"ymax");
      int nfreq=getIntAttribute(estNode,"nfreq");
      if (nfreq==0) nfreq=simInfo.getNSlice();
      int nstride=getIntAttribute(estNode,"nstride");
      if (nstride==0) nstride=1;
      manager->add(new ConductivityEstimator2D(simInfo, xmin, xmax, ymin, ymax, nfreq,
                       dim1, dim2, nxbin, nybin, nxdbin, nydbin, nstride ,mpi));
    }
    if (name=="SpinChargeEstimator") {
      int nbin=getIntAttribute(estNode,"nbin");
      if (nbin==0) nbin=1;
      int ndbin=getIntAttribute(estNode,"ndbin");
      if (ndbin==0) ndbin=1;
      int nfreq=getIntAttribute(estNode,"nfreq");
      if (nfreq==0) nfreq=simInfo.getNSlice();
      int nstride=getIntAttribute(estNode,"nstride");
      if (nstride==0) nstride=1;
      std::string speciesUp=getStringAttribute(estNode,"speciesUp");
      std::string speciesDown=getStringAttribute(estNode,"speciesDown");
      const Species &sup(simInfo.getSpecies(speciesUp));
      const Species &sdn(simInfo.getSpecies(speciesDown));
      manager->add(new SpinChargeEstimator(
                         simInfo,sup,sdn,nfreq,nbin,ndbin,nstride,mpi));
    }
    if (name=="ConductanceEstimator") {
      int nfreq=getIntAttribute(estNode,"nfreq");
      if (nfreq==0) nfreq=simInfo.getNSlice();
      std::string component=getStringAttribute(estNode,"component");
      int idim=-1; // Default use all components in tensor.
      if (component=="all") idim=-1;
      else if (component=="x") idim=0;
      else if (component=="y") idim=1;
      else if (component=="z") idim=2;
      bool useCharge=getBoolAttribute(estNode,"useCharge");
      bool useSpeciesTensor=getBoolAttribute(estNode,"useSpeciesTensor");
      int norder=getIntAttribute(estNode,"norder");
      if (norder==0) norder=1;
      manager->add(new ConductanceEstimator(simInfo,nfreq,0,
                           useSpeciesTensor,idim,useCharge,mpi,norder));
    }
    if (name=="DensityEstimator" || name=="DensCountEstimator"||
        name=="DensDensEstimator" || name=="CountCountEstimator") {
      //bool useCharge=getBoolAttribute(estNode,"useCharge");
      std::string species=getStringAttribute(estNode,"species");
      const Species *spec = 0;
      if (species!="" && species!="all") spec=&simInfo.getSpecies(species);
      std::string species1=getStringAttribute(estNode,"species1");
      const Species *spec1 = 0;
      if (species1!="" && species1!="all") spec1=&simInfo.getSpecies(species1);
      std::string species2=getStringAttribute(estNode,"species2");
      const Species *spec2 = 0;
      if (species2!="" && species2!="all") spec2=&simInfo.getSpecies(species2);
      if (spec1==0) spec1=spec;
      if (spec2==0) spec2=spec;
      std::string estName=getStringAttribute(estNode,"name");
      if (estName=="") {
        if (name=="DensityEstimator") {
          estName = "rho";
        } else if (name=="DensDensEstimator") {
          estName = "nn";
        } else {
          estName = "count";
        }
        if (spec) estName += species;
      }
      DensityEstimator::DistArray dist;
      std::vector<double> min_,max_;
      std::vector<int> nbin_;
      parseDistance(estNode,ctxt,dist,min_,max_,nbin_);
      IVec nbin;
      Vec min,max;
      if (dist.size()==0) {
        nbin = getIVecAttribute(estNode,"n");
        double a = getLengthAttribute(estNode,"a");
        min=-0.5*a*nbin;
        max=0.5*a*nbin;
        for (int i=0; i<NDIM; ++i) dist.push_back(new Cart(i));
      } else {
        for (unsigned int i=0; i<NDIM; ++i) {
          if (i<dist.size()) {
            min[i]=min_[i]; max[i]=max_[i]; nbin[i]=nbin_[i];
          } else {
            min[i]=0; max[i]=1; nbin[i]=1; dist.push_back(new Distance());
          }
        }
      }
      if (name=="DensityEstimator") {
        manager->add(new DensityEstimator(simInfo,estName,spec,
                                          min,max,nbin,dist, mpi));
      } else if (name=="DensDensEstimator") {
        DensDensEstimator::IVecN nbinN;
        int nfreq=getIntAttribute(estNode,"nfreq");
        if (nfreq==0) nfreq=simInfo.getNSlice();
        int nstride=getIntAttribute(estNode,"nstride");
        if (nstride==0) nstride=1;
        for (int i=0; i<NDIM; ++i) {
          nbinN[i]=nbin[i];
          nbinN[i+NDIM]=nbin[i];
        }
        nbinN[2*NDIM]=nfreq;
        manager->add(new DensDensEstimator(simInfo,estName,spec1,spec2,
                                  min,max,nbin, nbinN,dist,nstride,mpi));
      } else {
        int maxCount=getIntAttribute(estNode,"maxCount");
        if (maxCount==0) maxCount=1;
        if (name=="DensCountEstimator") {
          DensCountEstimator::IVecN nbinN;
          for (int i=0; i<NDIM; ++i) nbinN[i]=nbin[i];
          nbinN[NDIM]=maxCount+1;
          manager->add(new DensCountEstimator(simInfo,estName,spec,
                                              min,max,nbin,nbinN,dist,mpi));
        } else {
          CountCountEstimator::IVecN nbinN;
          for (int i=0; i<NDIM; ++i) {
            nbinN[i]=nbin[i];
            nbinN[i+NDIM]=nbin[i];
          }
          nbinN[2*NDIM]=nbinN[2*NDIM+1]=maxCount+1;
          int nfreq=getIntAttribute(estNode,"nfreq");
          if (nfreq==0) nfreq=simInfo.getNSlice();
          int nstride=getIntAttribute(estNode,"nstride");
          if (nstride==0) nstride=1;
          nbinN[2*NDIM+2]=nfreq;
          manager->add(new CountCountEstimator(simInfo,estName,spec,
                             min,max,nbin,nbinN,dist,nstride,mpi));
        }
      }
    }
    if (name=="PairCFEstimator") {
      ctxt->node = estNode;
      xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"*",ctxt);
      int N=obj->nodesetval->nodeNr;
      switch (N) {
        case 1: manager->add(parsePairCF<1>(estNode,obj)); break;
        case 2: manager->add(parsePairCF<2>(estNode,obj)); break;
        case 3: manager->add(parsePairCF<3>(estNode,obj)); break;
        case 4: manager->add(parsePairCF<4>(estNode,obj)); break;
        case 5: manager->add(parsePairCF<5>(estNode,obj)); break;
      }
    }
    if (name=="ZeroVarDensityEstimator") {
      
      std::string name=getStringAttribute(estNode,"name");
      int nspecies = getIntAttribute(estNode,"nspecies");  

      if (nspecies == 0) nspecies=2;
      Species *speciesList = new Species [nspecies];
      for (int ispec=0; ispec<nspecies; ispec++){
	std::stringstream sispec;
	sispec << "species"<<(ispec+1);
	std::string speciesName=getStringAttribute(estNode,sispec.str());
	std :: cout<<"Picked species "<< speciesName <<" for Zero Variance pairCF."<<std :: endl;
	speciesList[ispec]=simInfo.getSpecies(speciesName);
      }
      
      int nbin=1; double max=1.; double min=0.;
      ctxt->node = estNode;
      xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"*",ctxt);
      int N=obj->nodesetval->nodeNr; 
      for (int idist=0; idist<N; ++idist) {
	xmlNodePtr distNode=obj->nodesetval->nodeTab[idist];
	std::string name=getName(distNode);
	if (name=="Radial") {
	  min=getDoubleAttribute(distNode,"min");
	  max=getDoubleAttribute(distNode,"max");
	  nbin=getIntAttribute(distNode,"nbin");
	}
      }
   
      // delete [] speciesList;
      manager->add(new ZeroVarDensityEstimator(simInfo,name,speciesList,nspecies,min,max,nbin,action,doubleAction,mpi));
    }

    if (name=="DynamicPCFEstimator") {
      std::string name=getStringAttribute(estNode,"name");
      std::string species1=getStringAttribute(estNode,"species1");
      std::string species2=getStringAttribute(estNode,"species2");
      const Species &s1(simInfo.getSpecies(species1));
      const Species &s2(simInfo.getSpecies(species2));
      int nfreq=getIntAttribute(estNode,"nfreq");
      if (nfreq==0) nfreq=simInfo.getNSlice();
      int nstride=getIntAttribute(estNode,"nstride");
      if (nstride==0) nstride=1;
      std::vector<PairDistance*> dist;
      std::vector<double> min,max;
      std::vector<int> nbin;
      parsePairDistance(estNode,ctxt,dist,min,max,nbin);
      if (dist.size()==1) {
        manager->add(new DynamicPCFEstimator(simInfo,name,&s1,&s2,
             min[0], max[0], nbin[0], nfreq, nstride, dist[0], mpi));
      }
    }
    if (name=="PermutationEstimator") {
      std::string name=getStringAttribute(estNode,"name");
      std::string species1=getStringAttribute(estNode,"species1");
      if (species1=="") 
        species1=getStringAttribute(estNode,"species");
      const Species &s1(simInfo.getSpecies(species1));
      if (species1=="") 
        species1=s1.name;
      if (name=="") name = "perm_" + species1;
      manager->add(new PermutationEstimator(simInfo, name, s1, mpi));
    }
    if (name=="JEstimator") {
      int nBField=getIntAttribute(estNode,"nBField");
      double bmax=getDoubleAttribute(estNode,"bmax");
      manager->add(new JEstimator(simInfo,nBField,bmax,mpi));
    }
    if (name=="DiamagneticEstimator") {
      std::string unitName=getStringAttribute(estNode,"unit");
      double perN=getDoubleAttribute(estNode,"perN");
      double scale = getDoubleAttribute(estNode,"scale");
      if (perN > 0.) scale/=perN;
      if (unitName!="") scale/=simInfo.getUnits()->getLengthScaleIn(unitName,3);
      manager->add(new DiamagneticEstimator(simInfo,simInfo.getTemperature(),
                                            unitName,scale));
    }
    if (name=="WindingEstimator") {
      int nmax=getIntAttribute(estNode,"nmax");
      bool isChargeCoupled = getBoolAttribute(estNode,"isChargeCoupled");
      std::string name = getStringAttribute(estNode,"name");
      if (name=="") name = isChargeCoupled ? "charge_winding" : "winding";
      manager->add(new WindingEstimator(simInfo,nmax,name,isChargeCoupled,mpi));
    }
    if (name=="SKOmegaEstimator") {
      IVec nbin = getIVecAttribute(estNode,"n");
      blitz::TinyVector<int,NDIM+3> nbinN;
      nbinN[0] = nbinN[1] = simInfo.getNSpecies();
      for (int i=0; i<NDIM; ++i) nbinN[i+2] = (nbin[i]==0)?1:nbin[i];
      nbinN[NDIM+2] = getIntAttribute(estNode,"nfreq"); 
      int nstride = getIntAttribute(estNode,"nstride"); 
      if (nstride==0) nstride=1;
      std::string name=getStringAttribute(estNode,"name");
      if (name=="") name="skomega";
      manager->add(new SKOmegaEstimator(simInfo,name,nbin,nbinN,nstride,mpi));
    }
    if (name=="EMARateEstimator") {
      double C = getDoubleAttribute(estNode,"c");
      manager->add(new EMARateEstimator(simInfo,C));
    }
  }
  xmlXPathFreeObject(obj);
}

template<int N>
PairCFEstimator<N>* EstimatorParser::parsePairCF(xmlNodePtr estNode,
    xmlXPathObjectPtr obj) {
  std::string name=getStringAttribute(estNode,"name");
  // std::string species1=getStringAttribute(estNode,"species1");
  //const Species &s1(simInfo.getSpecies(species1));

  int nspecies = getIntAttribute(estNode,"nspecies");
  if (nspecies == 0) nspecies=2;
  Species *speciesList = new Species [nspecies];
  for (int ispec=0; ispec<nspecies; ispec++){
    std::stringstream sispec;
    sispec << "species"<<(ispec+1);
    std::string speciesName=getStringAttribute(estNode,sispec.str());
    std :: cout<<"Picked species "<< speciesName <<" for pairCF."<<std :: endl;
    speciesList[ispec]=simInfo.getSpecies(speciesName);
    }
  
  typename PairCFEstimator<N>::VecN min(0.), max(1.);
  typename PairCFEstimator<N>::IVecN nbin(1);
  typename PairCFEstimator<N>::DistN dist(N,(PairDistance*)0);
  for (int idist=0; idist<N; ++idist) {
    xmlNodePtr distNode=obj->nodesetval->nodeTab[idist];
    std::string name=getName(distNode);
    if (name=="Cartesian" || name=="Cartesian1" || name=="Cartesian2") {
      std::string dirName = getStringAttribute(distNode,"dir");
      int idir=0; if (dirName=="y") idir=1; else if (dirName=="z") idir=2;
      if (name=="Cartesian") {
        dist[idist]=new PairCart(idir);
      } else if (name=="Cartesian1") {
        dist[idist]=new PairCart1(idir);
      } else if (name=="Cartesian2") {
        dist[idist]=new PairCart2(idir);
      }
      min[idist] = getLengthAttribute(distNode,"min");
      max[idist] = getLengthAttribute(distNode,"max");
      nbin[idist] = getIntAttribute(distNode,"nbin");
    } else if (name=="Radial"||name=="Radial1"||name=="Radial2") {
      std::string dirName = getStringAttribute(distNode,"dir");
      int idir=-1;
      if (dirName=="x") idir=0;
      else if (dirName=="y") idir=1;
      else if (dirName=="z") idir=2;
      if (name=="Radial") {
        dist[idist] = new PairRadial(idir);
      } else if (name=="Radial1") {
        dist[idist] = new PairRadial1(idir);
      } else if (name=="Radial2") {
        dist[idist] = new PairRadial2(idir);
      }
      min[idist] = getLengthAttribute(distNode,"min");
      max[idist] = getLengthAttribute(distNode,"max");
      nbin[idist] = getIntAttribute(distNode,"nbin");
    } else if (name=="Angle" || name=="Angle1" || name=="Angle2") {
      int idim=0, jdim=1;
      if (NDIM==1) {
        jdim=0;
      } else if (NDIM>2) {
        std::string dirName = getStringAttribute(distNode,"dir");
        if (dirName=="x") {idim=1; jdim=2;}
        else if (dirName=="y") {idim=0; jdim=2;}
      }
      if (name=="Angle") {
        dist[idist]=new PairAngle(idim,jdim);
      } else if (name=="Angle1") {
        dist[idist]=new PairAngle1(idim,jdim);
      } else if (name=="Angle2") {
        dist[idist]=new PairAngle2(idim,jdim);
      }
      double minv = getDoubleAttribute(distNode,"min");
      double maxv = getDoubleAttribute(distNode,"max");
      const double PI=3.14159265358793;
      if (fabs(minv-maxv)<1e-9) {minv=-PI; maxv=+PI;}
      min[idist]=minv;
      max[idist]=maxv;
      nbin[idist]=getIntAttribute(distNode,"nbin");
std::cout << name << min[idist] << " - " << max[idist] << "  " << nbin << std::endl;
    }
  }
  return new PairCFEstimator<N>(simInfo,name,speciesList,nspecies,min,max,nbin,dist,mpi);
  delete [] speciesList;
}

void EstimatorParser::parseDistance(xmlNodePtr estNode, 
    const xmlXPathContextPtr& ctxt,
    std::vector<Distance*> &darray, std::vector<double> &min,
    std::vector<double> &max, std::vector<int>& nbin) {
  ctxt->node = estNode;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"*",ctxt);
  int N=obj->nodesetval->nodeNr;
  int idir=0;
  for (int idist=0; idist<N; ++idist) {
    xmlNodePtr distNode=obj->nodesetval->nodeTab[idist];
    std::string name=getName(distNode);
    if (name=="Cartesian") {
      std::string dirName = getStringAttribute(distNode,"dir");
      for (int i=0; i<NDIM; ++i) if (dirName==dimName.substr(i,1)) idir=i;
      darray.push_back(new Cart(idir));
      min.push_back(getLengthAttribute(distNode,"min"));
      max.push_back(getLengthAttribute(distNode,"max"));
      nbin.push_back(getIntAttribute(distNode,"nbin"));
      idir++;
    } else if (name=="Radial") {
      int idir=-1;
      std::string dirName = getStringAttribute(distNode,"dir");
      for (int i=0; i<NDIM; ++i) if (dirName==dimName.substr(i,1)) idir=i;
      darray.push_back(new Radial(idir));
      min.push_back(getLengthAttribute(distNode,"min"));
      max.push_back(getLengthAttribute(distNode,"max"));
      nbin.push_back(getIntAttribute(distNode,"nbin"));
    }
  }
}

void EstimatorParser::parsePairDistance(xmlNodePtr estNode, 
    const xmlXPathContextPtr& ctxt,
    std::vector<PairDistance*> &darray, std::vector<double> &min,
    std::vector<double> &max, std::vector<int>& nbin) {
  ctxt->node = estNode;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"*",ctxt);
  int N=obj->nodesetval->nodeNr;
  int idir=0;
  for (int idist=0; idist<N; ++idist) {
    xmlNodePtr distNode=obj->nodesetval->nodeTab[idist];
    std::string name=getName(distNode);
    if (name=="Cartesian" || name=="Cartesian1" || name=="Cartesian2") {
      std::string dirName = getStringAttribute(distNode,"dir");
      for (int i=0; i<NDIM; ++i) if (dirName==dimName.substr(i,1));
      if (name=="Cartesian") {
        darray.push_back(new PairCart(idir));
      } else if (name=="Cartesian1") {
        darray.push_back(new PairCart1(idir));
      } else {
        darray.push_back(new PairCart2(idir));
      }
      min.push_back(getLengthAttribute(distNode,"min"));
      max.push_back(getLengthAttribute(distNode,"max"));
      nbin.push_back(getIntAttribute(distNode,"nbin"));
      idir++;
    } else if (name=="Radial"||name=="Radial1"||name=="Radial2") {
      int idir=-1;
      std::string dirName = getStringAttribute(distNode,"dir");
      for (int i=0; i<NDIM; ++i) if (dirName==dimName.substr(i,1)) idir=i;
      if (name=="Radial") {
        darray.push_back(new PairRadial(idir));
      } else if (name=="Radial1") {
        darray.push_back(new PairRadial1(idir));
      } else if (name=="Radial2") {
        darray.push_back(new PairRadial2(idir));
      }
      min.push_back(getLengthAttribute(distNode,"min"));
      max.push_back(getLengthAttribute(distNode,"max"));
      nbin.push_back(getIntAttribute(distNode,"nbin"));
    } else if (name=="Angle" || name=="Angle1" || name=="Angle2") {
      int idim=0, jdim=1;
      if (NDIM==1) {
        jdim=0;
      } else if (NDIM>2) {
        std::string dirName = getStringAttribute(distNode,"dir");
        if (dirName=="x") {idim=1; jdim=2;}
        else if (dirName=="y") {idim=0; jdim=2;}
      }
      if (name=="Angle") {
        darray.push_back(new PairAngle(idim,jdim));
      } else if (name=="Angle1") {
        darray.push_back(new PairAngle1(idim,jdim));
      } else if (name=="Angle2") {
        darray.push_back(new PairAngle2(idim,jdim));
      }
      double minv = getDoubleAttribute(distNode,"min");
      double maxv = getDoubleAttribute(distNode,"max");
      const double PI=3.14159265358793;
      if (fabs(minv-maxv)<1e-9) {minv=-PI; maxv=+PI;}
      min.push_back(minv);
      max.push_back(maxv);
      nbin.push_back(getIntAttribute(distNode,"nbin"));
    }
  }
}
