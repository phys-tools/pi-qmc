// $Id$
/*  Copyright (C) 2004-2008 John B. Shumway, Jr.

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
#include "AngularMomentumEstimator.h"
#include "CoulombAction.h"
#include "ConductivityEstimator.h"
#include "ConductanceEstimator.h"
#include "DensDensEstimator.h"
#include "DensityEstimator.h"
#include "DensCountEstimator.h"
#include "Distance.h"
#include "CoulombEnergyEstimator.h"
#include "EwaldCoulombEstimator.h"
#include "DoubleAction.h"
#include "ThermoEnergyEstimator.h"
#include "VirialEnergyEstimator.h"
#include "DipoleMomentEstimator.h"
#include "GreenDipoleMEstimator.h"
#include "QuadrupoleMomentEstimator.h"
#include "VVCorrelationEstimator.h"
#include "BondLengthEstimator.h"
#include "BondLengthTwoEstimator.h"
#include "FrequencyEstimator.h"
#include "spin/SpinEstimator.h"
#include "PositionEstimator.h"
#include "BoxEstimator.h"
#include "PairCFEstimator.h"
#include "PermutationEstimator.h"
#include "JEstimator.h"
#include "stats/MPIManager.h"
#include "SuperCell.h"
#include "SpinChargeEstimator.h"
#include "SHOPhase.h"
#include "VIndEstimator.h"
#include "EIndEstimator.h"
#include "SimInfoWriter.h"
#include "stats/Units.h"

EstimatorParser::EstimatorParser(const SimulationInfo& simInfo,
    const double tau, const Action* action, const DoubleAction* doubleAction,
    MPIManager *mpi)
  : XMLUnitParser(simInfo.getUnits()), manager(0),
    simInfo(simInfo), tau(tau), action(action), doubleAction(doubleAction),
    mpi(mpi) {
  manager=new EstimatorManager("pimc.h5", mpi, new SimInfoWriter(simInfo));
}

EstimatorParser::~EstimatorParser() {delete manager;}

void EstimatorParser::parse(const xmlXPathContextPtr& ctxt) {
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
        manager->add(new EwaldCoulombEstimator(simInfo,action,
                       epsilon,rcut,kcut,mpi,unitName,scale,shift));
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
    if (name=="GreenDipoleMEstimator") {
      std::string component=getStringAttribute(estNode,"component");
      int index=2; //default use z component
      if (component=="x") index=0;
      else if (component=="y") index=1;
      manager->add(new GreenDipoleMEstimator(simInfo, index, mpi));
    }	
    if (name=="QuadrupoleMomentEstimator") {
      std::string component=getStringAttribute(estNode,"component");
      int index=5; //default use zz component (Q33)
      if (component=="xx") index=0;
      else if (component=="xy") index=1;
      else if (component=="xz") index=2;
      else if (component=="yy") index=3;
      else if (component=="yz") index=4;
      else if (component=="zz") index=5;
      manager->add(new QuadrupoleMomentEstimator(simInfo, index));
    }
    if (name=="VVCorrelationEstimator") {
      manager->add(new VVCorrelationEstimator(simInfo));
    }
    if (name=="BondLengthEstimator") {
      std::string species1=getStringAttribute(estNode,"species1");
      std::string species2=getStringAttribute(estNode,"species2");
      manager->add(new BondLengthEstimator(simInfo,simInfo.getSpecies(species1)
			      ,simInfo.getSpecies(species2)));
    }
    if (name=="BondLengthTwoEstimator") {
      std::string species1=getStringAttribute(estNode,"species1");
      std::string species2=getStringAttribute(estNode,"species2");
      manager->add(new BondLengthTwoEstimator(simInfo,
                               simInfo.getSpecies(species1),
			       simInfo.getSpecies(species2)));
    }
    if (name=="FrequencyEstimator") {
      std::string species1=getStringAttribute(estNode,"species1");
      std::string species2=getStringAttribute(estNode,"species2");
      manager->add(new FrequencyEstimator(simInfo,simInfo.getSpecies(species1)
			      ,simInfo.getSpecies(species2), mpi));
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
    if (name=="SpinChargeEstimator") {
      int nbin=getIntAttribute(estNode,"nbin");
      if (nbin==0) nbin=1;
      int ndbin=getIntAttribute(estNode,"ndbin");
      if (ndbin==0) ndbin=1;
      int nfreq=getIntAttribute(estNode,"nfreq");
      if (nfreq==0) nfreq=simInfo.getNSlice();
      std::string speciesUp=getStringAttribute(estNode,"speciesUp");
      std::string speciesDown=getStringAttribute(estNode,"speciesDown");
      const Species &sup(simInfo.getSpecies(speciesUp));
      const Species &sdn(simInfo.getSpecies(speciesDown));
      manager->add(new SpinChargeEstimator(
                         simInfo,sup,sdn,nfreq,nbin,ndbin,mpi));
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
    if (name=="DensityEstimator" || name=="DensCountEstimator") {
      //bool useCharge=getBoolAttribute(estNode,"useCharge");
      std::string species=getStringAttribute(estNode,"species");
      const Species *spec = 0;
      if (species!="" && species!="all") spec=&simInfo.getSpecies(species);
      std::string estName=getStringAttribute(estNode,"name");
      if (estName=="") {
        if (name=="DensityEstimator") {
          estName = "rho";
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
      } else {
        int maxCount=getIntAttribute(estNode,"maxCount");
        if (maxCount==0) maxCount=1;
        DensCountEstimator::IVecN nbinN;
        for (int i=0; i<NDIM; ++i) nbinN[i]=nbin[i];
        nbinN[NDIM]=maxCount+1;
        manager->add(new DensCountEstimator(simInfo,estName,spec,
                                            min,max,nbin,nbinN,dist,mpi));
      }
    }
    if (name=="DensDensEstimator") {
      int nbin=getIntAttribute(estNode,"nbin");
      if (nbin==0) nbin=1;
      int ndbin=getIntAttribute(estNode,"ndbin");
      if (ndbin==0) ndbin=1;
      manager->add(new DensDensEstimator(simInfo,action,doubleAction,
                                         nbin,ndbin,mpi));
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
    if (name=="PermutationEstimator") {
      manager->add(new PermutationEstimator(simInfo,mpi));
    }
    if (name=="JEstimator") {
      int nBField=getIntAttribute(estNode,"nBField");
      double bmax=getDoubleAttribute(estNode,"bmax");
      manager->add(new JEstimator(simInfo,nBField,bmax,mpi));
    }
  }
  xmlXPathFreeObject(obj);
}

template<int N>
PairCFEstimator<N>* EstimatorParser::parsePairCF(xmlNodePtr estNode,
    xmlXPathObjectPtr obj) {
  std::string name=getStringAttribute(estNode,"name");
  std::string species1=getStringAttribute(estNode,"species1");
  std::string species2=getStringAttribute(estNode,"species2");
  const Species &s1(simInfo.getSpecies(species1));
  const Species &s2(simInfo.getSpecies(species2));
  typename PairCFEstimator<N>::VecN min(0.), max(1.);
  typename PairCFEstimator<N>::IVecN nbin(1);
  typename PairCFEstimator<N>::DistN 
    dist(N,(typename PairCFEstimator<N>::Dist*)0);
  for (int idist=0; idist<N; ++idist) {
    xmlNodePtr distNode=obj->nodesetval->nodeTab[idist];
    std::string name=getName(distNode);
    if (name=="Cartesian" || name=="Cartesian1" || name=="Cartesian2") {
      std::string dirName = getStringAttribute(distNode,"dir");
      int idir=0; if (dirName=="y") idir=1; else if (dirName=="z") idir=2;
      if (name=="Cartesian") {
        dist[idist]=new typename PairCFEstimator<N>::Cart(idir);
      } else if (name=="Cartesian1") {
        dist[idist]=new typename PairCFEstimator<N>::Cart1(idir);
      } else if (name=="Cartesian2") {
        dist[idist]=new typename PairCFEstimator<N>::Cart2(idir);
      }
      min[idist] = getLengthAttribute(distNode,"min");
      max[idist] = getLengthAttribute(distNode,"max");
      nbin[idist] = getIntAttribute(distNode,"nbin");
    } else if (name=="Radial") {
      std::string dirName = getStringAttribute(distNode,"dir");
      int idir=-1;
      if (dirName=="x") idir=0;
      else if (dirName=="y") idir=1;
      else if (dirName=="z") idir=2;
      dist[idist] = new typename PairCFEstimator<N>::Radial(idir);
      min[idist] = getLengthAttribute(distNode,"min");
      max[idist] = getLengthAttribute(distNode,"max");
      nbin[idist] = getIntAttribute(distNode,"nbin");
    }
  }
  return new PairCFEstimator<N>(simInfo,name,s1,s2,min,max,nbin,dist,mpi);
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
      for (int i=0; i<NDIM; ++i) if (dirName==dimName.substr(i,1));
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
