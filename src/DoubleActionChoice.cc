// $Id: DoubleActionChoice.cc 185 2009-10-13 06:00:07Z john.shumwayjr $
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
#include "DoubleActionChoice.h"
#include "DoubleMLSampler.h"
#include "DoubleSectionChooser.h"
#include "Paths.h"

DoubleActionChoice::DoubleActionChoice(const int n)
  : CompositeDoubleAction(n) {
}

DoubleActionChoice::~DoubleActionChoice() {
}

double DoubleActionChoice::getActionDifference(
    const DoubleMLSampler& sampler, const int level) {
  imodel = sampler.getPaths().getModelState();
  double diff = actions[imodel]->getActionDifference(sampler,level);
  return diff;
}

double DoubleActionChoice::getActionDifference(const Paths &paths, 
    const VArray &displacement, int nmoving, const IArray &movingIndex, 
    int iFirstSlice, int nslice) {
  imodel = paths.getModelState();
  double diff = actions[imodel]->getActionDifference(paths,displacement,nmoving,
                                           movingIndex,iFirstSlice,nslice);
  return diff;
}

double DoubleActionChoice::getTotalAction(const Paths& paths, int level) const {
  double total = actions[imodel]->getTotalAction(paths,level);
  return total;
}

void DoubleActionChoice::getBeadAction(const Paths& paths, 
       int ipart, int islice, double& u, double& utau, double& ulambda,
       Action::Vec& fm, Action::Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  actions[imodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
}

void DoubleActionChoice::initialize(const DoubleSectionChooser& 
                                             sectionChooser) {
  imodel = sectionChooser.getPaths().getModelState();
  actions[imodel]->initialize(sectionChooser);
}

void DoubleActionChoice::acceptLastMove() {
  actions[imodel]->acceptLastMove();
}

double DoubleActionChoice::getActionDifference(const Paths &paths, int j) {
 jmodel = j;
 paths.sumOverLinks(*this);
 return actionDifference;
}

void DoubleActionChoice::initCalc(const int nslice, const int firstSlice) {
  actionDifference = 0.;
}

void DoubleActionChoice::handleLink(const LinkSummable::Vec &start,
    const LinkSummable::Vec &end, const int ipart, const int islice,
    const Paths &paths) {
  double u=0., utau=0, ulambda=0;
  LinkSummable::Vec fm=0.; LinkSummable::Vec fp=0.;
  actions[imodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
  actionDifference -= u;
  actions[jmodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
  actionDifference += u;
}

void DoubleActionChoice::endCalc(const int nslice) {
}

