// $Id: ActionParser.cc,v 1.86 2008/11/23 15:55:21 jshumwa Exp $
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

#include<Beads.h>


#include "ActionParser.h"
#include "Action.h"
#include "SimulationInfo.h"
#include "SphereAction.h"
#include "CompositeAction.h"
#include "CompositeDoubleAction.h"
#include "spin/SpinAction.h"
#include "SpringAction.h"
#include "SpringTensorAction.h"
#include "CoulombAction.h"
#include "GaussianAction.h"
#include "GaussianDotAction.h"
#include "OpticalLatticeAction.h"
#include "FixedNodeAction.h"
#include "FixedPhaseAction.h"
#include "spin/SpinFixedPhaseAction.h"
#include "spin/LatticeSpinPhase.h"
#include "spin/SpinPhase.h"
#include "FreeParticleNodes.h"
#include "ExcitonNodes.h"
#include "FreePartNodesNoUpdate.h"
#include "AnisotropicNodes.h"
#include "PairAction.h"
#include "EmpiricalInteraction.h"
#include "GridPotential.h"
#include "GroundStateSNode.h"
#include "GroundStateWFNodes.h"
#include "HyperbolicAction.h"
#include "QPCAction.h"
#include "TimpQPC.h"
#include "TwoQDAction.h"
#include "JelliumSlab.h"
#include "SmoothedGridPotential.h"
#include "SHOAction.h"
#include "SHODotAction.h"
#include "SHONodes.h"
#include "WireNodes.h"
#include "SHOPhase.h"
#include "Spin4DPhase.h"
#include "PrimSHOAction.h"
#include "EFieldAction.h"
#include "EwaldAction.h"
#include "StillWebAction.h"
#include "stats/MPIManager.h"
//#include "GrapheneAction.h"
#include "WellImageAction.h"
#include <iostream>

ActionParser::ActionParser(const SimulationInfo& simInfo, const int maxlevel,
  MPIManager *mpi)
  : XMLUnitParser(simInfo.getUnits()), 
    action(0), doubleAction(0), simInfo(simInfo), tau(simInfo.getTau()),
    maxlevel(maxlevel), mpi(mpi) {
}

void ActionParser::parse(const xmlXPathContextPtr& ctxt) {
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Action/*",ctxt);
  int naction=obj->nodesetval->nodeNr;
  CompositeAction* composite=new CompositeAction(naction); action=composite;
  CompositeDoubleAction* doubleComposite=new CompositeDoubleAction(naction);
  doubleAction=doubleComposite;
  for (int iaction=0; iaction<naction; ++iaction) {
    xmlNodePtr actNode=obj->nodesetval->nodeTab[iaction];
    std::string name=getName(actNode);
    //std::cout << "Action node name: " << name << std::endl;
    if (name=="SpringAction") {
      double pgDelta=getDoubleAttribute(actNode,"pgDelta");
      if (pgDelta==0) pgDelta=0.02;
      composite->addAction(new SpringAction(simInfo,maxlevel,pgDelta));
      continue;
    } else if (name=="HyperbolicAction") {
      composite->addAction(new HyperbolicAction(simInfo,maxlevel));
      continue;
    } else if (name=="CoulombAction") {
      double epsilon=getDoubleAttribute(actNode,"epsilon");
      if (epsilon==0) epsilon=1;
      int norder=getIntAttribute(actNode,"norder");
      double rmin=getLengthAttribute(actNode,"rmin");
      double rmax=getLengthAttribute(actNode,"rmax");
      int ngpts=getIntAttribute(actNode,"ngridPoints");
      bool dumpFiles=getBoolAttribute(actNode,"dumpFiles");
      bool useEwald=getBoolAttribute(actNode,"useEwald");
      int ndim=getIntAttribute(actNode,"ewaldNDim");
      if (ndim==0) ndim=NDIM;
      double kcut=getDoubleAttribute(actNode,"ewaldKcut");
      double screenDist=getLengthAttribute(actNode,"screenDist");
      composite->addAction(
        new CoulombAction(epsilon,simInfo,norder,rmin,rmax,ngpts,dumpFiles,
                          useEwald,ndim,kcut,screenDist));
      continue;
    } else if (name=="GaussianAction") {
      double v0=getEnergyAttribute(actNode,"v0");
      double alpha=getDoubleAttribute(actNode,"alpha");
      composite->addAction(new GaussianAction(v0,alpha,simInfo));
      continue;
    } else if (name=="GridPotential") {
      std::string fileName=getStringAttribute(actNode,"file");
      if (fileName=="") fileName="emagrids.h5";
      composite->addAction(new GridPotential(simInfo,fileName));
      continue;
    } else if (name=="SmoothedGridPotential") {
      std::string fileName=getStringAttribute(actNode,"file");
      if (fileName=="") fileName="emagrids.h5";
      int maxLevel=getIntAttribute(actNode,"level");
      if (maxLevel==0) maxLevel=12; //maybe this should be 8
      composite->addAction(new SmoothedGridPotential(simInfo,maxLevel,fileName));
      continue;
    } else if (name=="OpticalLatticeAction") {
      Vec v0;
      for (int idim=0; idim<NDIM; ++idim) {
        v0[idim]=getEnergyAttribute(actNode,std::string("v0")+
                                    std::string(dimName).substr(idim,1));
      }
      Vec length;
      for (int idim=0; idim<NDIM; ++idim) {
        length[idim]=getLengthAttribute(actNode,std::string("l")+
                                        std::string(dimName).substr(idim,1));
      }
      Vec max;
      for (int idim=0; idim<NDIM; ++idim) {
        max[idim]=getLengthAttribute(actNode,std::string("max")+
                                     std::string(dimName).substr(idim,1));
      }
      composite->addAction(new OpticalLatticeAction(v0,length,max,simInfo));
      continue;
    } else if (name=="PrimSHOAction") {
      double a=getDoubleAttribute(actNode,"a");
      double b=getDoubleAttribute(actNode,"b");
      double omega=getEnergyAttribute(actNode,"omega");
      if (a==0) a=0.5*omega*omega;
      int ndim=getIntAttribute(actNode,"ndim");
      if (ndim==0) ndim=NDIM;
      composite->addAction(new PrimSHOAction(a,b,simInfo,ndim));
      continue;
    } else if (name=="StillingerWeberAction") {
      composite->addAction(new StillWebAction(simInfo,"struct.h5"));
//    } else if (name=="GrapheneActionAction") {
//      composite->addAction(new GrapheneAction(simInfo));
    } else if (name=="SpinAction") {
      const double bx=getEnergyAttribute(actNode,"bx");
      const double by=getEnergyAttribute(actNode,"by");
      const double bz=getEnergyAttribute(actNode,"bz");
      const double gc=getEnergyAttribute(actNode,"gc");
      const double tau=simInfo.getTau();
      const double mass=simInfo.getSpecies(0).mass;
      composite->addAction(
        new SpinAction(tau,mass,bx,by,bz,simInfo.spinOmega,gc));
      continue;
    } else if (name=="SHOAction") {
      double omega=getEnergyAttribute(actNode,"omega");
      const double mass=simInfo.getSpecies(0).mass;
      if (omega==0) omega=1;
      int ndim=getIntAttribute(actNode,"ndim");
      if (ndim==0) ndim=NDIM;
      composite->addAction(new SHOAction(simInfo.getTau(),omega,mass,ndim));
      continue;
    } else if (name=="TwoQDAction") {
      const double omega=getEnergyAttribute(actNode,"omega");
      const double mass=simInfo.getSpecies(0).mass;
      const double d=getLengthAttribute(actNode,"d");
      double alpha=getDoubleAttribute(actNode,"alpha");
      if (alpha==0.) alpha=1.0;
      composite->addAction(
        new TwoQDAction(simInfo.getTau(),omega,mass,d,alpha));
      continue;
    } else if (name=="QPCAction") {
      double d=getLengthAttribute(actNode,"d");
      if (d==0) d=780.46;
      double v0=getEnergyAttribute(actNode,"v0");
      if (v0==0) v0=0.00011025;
      bool isWireOnly=getBoolAttribute(actNode,"isWireOnly");
      if (isWireOnly) {d=v0=0;}
      double omega=getEnergyAttribute(actNode,"omega");
      if (omega==0) omega=0.0000735;
      const double mass=simInfo.getSpecies(0).mass;
      composite->addAction(new QPCAction(simInfo.getTau(),d,v0,omega,mass));
      continue;
    } else if (name=="TimpQPC") {
      double w=getLengthAttribute(actNode,"w"); if (w==0) w=189.;
      double l=getLengthAttribute(actNode,"l"); if (l==0) l=587.;
      double vG=getEnergyAttribute(actNode,"vG"); if (vG==0) vG=-0.03675;
      double z=getLengthAttribute(actNode,"z"); if (z==0) z=189.;
      composite->addAction(new TimpQPC(*simInfo.getSuperCell(),
                                        simInfo.getTau(),w,l,vG,z,mpi));
    } else if (name=="JelliumSlab") {
      double qtot=getDoubleAttribute(actNode,"qtot");
      double thick=getLengthAttribute(actNode,"thickness");
      bool isSlabInMiddle=getBoolAttribute(actNode,"isSlabInMiddle");
      composite->addAction(new JelliumSlab(simInfo,qtot,thick,isSlabInMiddle));
    } else if (name=="SphereAction") {
      double radius=getLengthAttribute(actNode,"radius");
      std::string specName=getStringAttribute(actNode,"species");
      const Species& species(simInfo.getSpecies(specName));
      composite->addAction(new SphereAction(simInfo.getTau(),radius,species));
      continue;
    } else if (name=="GaussianDotAction") {
      double R=getLengthAttribute(actNode,"radius");
      double alpha=1./(R*R);
      Vec center;
      for (int idim=0; idim<NDIM; ++idim) {
        center[idim]=getLengthAttribute(actNode,
                                  std::string(dimName).substr(idim,1));
      }
      double v0=getEnergyAttribute(actNode,"v0");
      composite->addAction(new GaussianDotAction(v0,alpha,center,simInfo));
      continue;
    } else if (name=="SHODotAction") {
      double t=getLengthAttribute(actNode,"thickness");
      double v0=getEnergyAttribute(actNode,"v0");
      double k=getDoubleAttribute(actNode,"k");
      composite->addAction(new SHODotAction(simInfo.getTau(),t,v0,k));
      continue;
    } else if (name=="SpringTensorAction") {
      composite->addAction(new SpringTensorAction(simInfo));
      continue;
    }  else if (name=="FixedNodeAction") {
      ctxt->node=actNode;
      std::string specName=getStringAttribute(ctxt->node,"species");
      std::string modelName=getStringAttribute(ctxt->node,"model");
      const Species& species(simInfo.getSpecies(specName));
      double t=getEnergyAttribute(actNode,"temperature");
      bool noNodalAction=getBoolAttribute(actNode,"noNodalAction");
      bool useDistDerivative=getBoolAttribute(actNode,"useDistDerivative");
      if (t==0) t=simInfo.getTemperature();
      NodeModel *nodeModel=0;
      if (modelName=="SHONodes") {
        const double omega=getEnergyAttribute(ctxt->node,"omega");
        nodeModel=new SHONodes(simInfo,species,omega,t,maxlevel);
      } else if (modelName=="GSSNode" || modelName=="GroundStateSNode") {
        std::string centerName=getStringAttribute(ctxt->node,"center");
        const int icenter=simInfo.getSpecies(centerName).ifirst;
        nodeModel=new GroundStateSNode(species,icenter,simInfo.getTau());
      } else if (modelName=="WireNodes") {
        const double omega=getEnergyAttribute(ctxt->node,"omega");
        const bool updates=getBoolAttribute(ctxt->node,"useUpdates");
        int maxMovers=0;
        if (updates) maxMovers=3;
        nodeModel=new WireNodes(simInfo,species,omega,t,maxlevel,updates,
                                maxMovers);
      } else if (modelName=="ExcitonNodes") {
        const double radius=getLengthAttribute(ctxt->node,"radius");
        const bool updates=getBoolAttribute(ctxt->node,"useUpdates");
        int maxMovers=3;
        nodeModel=new ExcitonNodes(simInfo,species,t,maxlevel,radius,updates,
                                        maxMovers);
      } else {
        const bool updates=getBoolAttribute(ctxt->node,"useUpdates");
        int maxMovers=0;
        if (updates) maxMovers=3;
        nodeModel=new FreeParticleNodes(simInfo,species,t,maxlevel,updates,
                                        maxMovers);
      }
      doubleComposite->addAction(
               new FixedNodeAction(simInfo,species,nodeModel,!noNodalAction,
                                   useDistDerivative));
      continue;
    }  else if (name=="FixedPhaseAction") {
      ctxt->node=actNode;
      std::string specName=getStringAttribute(ctxt->node,"species");
      std::string modelName=getStringAttribute(ctxt->node,"model");
      const Species& species(simInfo.getSpecies(specName));
      double t=getEnergyAttribute(actNode,"temperature");
      if (t==0) t=simInfo.getTemperature();
      PhaseModel *phaseModel=0;
      if (modelName=="SHOPhase") {
        double omega=getEnergyAttribute(actNode,"omega");
	double b=getEnergyAttribute(actNode,"b");
        phaseModel=new SHOPhase(simInfo,species,omega,t,b,maxlevel);
      }else if(modelName=="SpinPhase"){
        double t=getEnergyAttribute(actNode,"t");
        if (t==0) t=simInfo.getTemperature();
        double bx=getEnergyAttribute(actNode,"bx");
        double by=getEnergyAttribute(actNode,"by");
        double bz=getEnergyAttribute(actNode,"bz");
        double gmubs=getDoubleAttribute(actNode,"gmubs");
        if (gmubs==0) gmubs=1.0;
        phaseModel=new Spin4DPhase(simInfo,species,t,bx,by,bz,gmubs,maxlevel);
      }  	      
        doubleComposite->addAction(
                              new FixedPhaseAction(simInfo,species,phaseModel));
    }  else if (name=="SpinFixedPhaseAction") {
      ctxt->node=actNode;
      std::string specName=getStringAttribute(ctxt->node,"species");
      std::string modelName=getStringAttribute(ctxt->node,"model");
      const Species& species(simInfo.getSpecies(specName));
      double t=getEnergyAttribute(actNode,"temperature");
      if (t==0) t=simInfo.getTemperature();
      double bx=getEnergyAttribute(actNode,"bx");
      double by=getEnergyAttribute(actNode,"by");
      double bz=getEnergyAttribute(actNode,"bz");
      double gmubs=getDoubleAttribute(actNode,"gmubs");
      if (gmubs==0) gmubs=1.0;
      SpinPhaseModel *phaseModel=0;
      if (modelName=="LatticeSpinPhase") {
        phaseModel=
          new LatticeSpinPhase(simInfo,species,t,bx,by,bz,gmubs,maxlevel);
      }else if(modelName=="SpinPhase") {
        phaseModel=new SpinPhase(simInfo,species,t,bx,by,bz,gmubs,maxlevel);
      }  	      
        doubleComposite->addAction(
                         new SpinFixedPhaseAction(simInfo,species,phaseModel));
    } else if (name=="FreePartNodesNoUpdate") {
      ctxt->node=actNode;
      std::string specName=getStringAttribute(ctxt->node,"species");
      const Species& species(simInfo.getSpecies(specName));
      double t=getEnergyAttribute(actNode,"temperature");
      if (t==0) t=simInfo.getTemperature();
      doubleComposite->addAction(new FreePartNodesNoUpdate(simInfo,species,t));
      continue;
    } else if (name=="PairAction") {
      ctxt->node=actNode;
      std::string specName=getStringAttribute(ctxt->node,"species1");
      const Species& species1(simInfo.getSpecies(specName));
      specName=getStringAttribute(ctxt->node,"species2");
      const Species& species2(simInfo.getSpecies(specName));
      std::string filename=getStringAttribute(ctxt->node,"file");
      int norder=getIntAttribute(actNode,"norder");
      bool isDMD=getBoolAttribute(actNode,"isDMD");
      composite->addAction(new PairAction(species1,species2,filename,
                                          simInfo,norder,isDMD));
      continue;
    } else if (name=="EmpiricalInteraction") {
      ctxt->node=actNode;
      std::string specName=getStringAttribute(ctxt->node,"species1");
      const Species& species1(simInfo.getSpecies(specName));
      specName=getStringAttribute(ctxt->node,"species2");
      const Species& species2(simInfo.getSpecies(specName));
      std::string modelName=getStringAttribute(ctxt->node,"model");
      EmpiricalInteraction::Potential* pot=0;
      int norder=getIntAttribute(actNode,"norder");
      double rmin=getLengthAttribute(actNode,"rmin");
      double rmax=getLengthAttribute(actNode,"rmax");
      int ngpts=getIntAttribute(actNode,"ngridPoints");
      if (modelName=="cosh2") {
        double v0=getEnergyAttribute(actNode,"v0");
        double kappa=getInvLengthAttribute(actNode,"kappa");
        pot = new EmpiricalInteraction::Cosh2Potential(v0,kappa); 
      }
      EmpiricalInteraction empAction(*pot,simInfo.getTau());
      composite->addAction(new PairAction(species1,species2,empAction,
                                          simInfo,norder,rmin,rmax,ngpts));
      delete pot;
      continue;
    } else if (name=="GroundStateWFNodes") {
      ctxt->node=actNode;
      std::string specName=getStringAttribute(ctxt->node,"species");
      const Species& species(simInfo.getSpecies(specName));
      specName=getStringAttribute(ctxt->node,"refSpecies");
      const Species& refSpecies(simInfo.getSpecies(specName));
      composite->addAction(new GroundStateWFNodes(species,refSpecies));
    } else if (name=="EFieldAction") {
      double scale=getDoubleAttribute(actNode,"scale");
      std::string component=getStringAttribute(actNode,"component");
      int index=2; //default use z component
      if (component=="x") index=0;
      else if (component=="y") index=1;
      composite->addAction(new EFieldAction(simInfo, scale, index));
      continue;
    } else if (name=="EwaldAction" && NDIM==3) {
      ctxt->node=actNode;
      composite->addAction(parseEwaldActions(ctxt));
      continue;
    } else if (name=="WellImageAction") {
      double epsIn=getDoubleAttribute(actNode,"epsIn");
      double epsOut=getDoubleAttribute(actNode,"epsOut");
      double width=getDoubleAttribute(actNode,"width");
      double z0=getDoubleAttribute(actNode,"z0");
      double delta=getDoubleAttribute(actNode,"delta");
      composite->addAction(new WellImageAction(simInfo,epsIn,epsOut,width,
                                 z0,delta)); 
      continue;
    }
  }
  xmlXPathFreeObject(obj);
  if (doubleComposite->getCount()==0) {
     doubleAction=0; delete doubleComposite;
  }
}

Action* ActionParser::parseEwaldActions(const xmlXPathContextPtr& ctxt) {
  EwaldAction* ewald=0;
#if NDIM==3
  double rcut=getLengthAttribute(ctxt->node,"rcut");
  double kcut=getDoubleAttribute(ctxt->node,"kcut");
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"PairAction",ctxt);
  int naction=obj->nodesetval->nodeNr;
  ewald=new EwaldAction(simInfo,rcut,kcut,naction);
  for (int iaction=0; iaction<naction; ++iaction) {
    xmlNodePtr actNode=obj->nodesetval->nodeTab[iaction];
    std::string specName=getStringAttribute(actNode,"species1");
    const Species& species1(simInfo.getSpecies(specName));
    specName=getStringAttribute(actNode,"species2");
    const Species& species2(simInfo.getSpecies(specName));
    std::string filename=getStringAttribute(actNode,"file");
    int norder=getIntAttribute(actNode,"norder");
    ewald->addAction(new PairAction(species1,species2,filename,simInfo,norder));
  }
  xmlXPathFreeObject(obj);
  ewald->setup();
#endif
  return ewald;
}

const std::string ActionParser::dimName="xyzklmnopqrstuv";
