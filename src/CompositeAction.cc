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
#include "CompositeAction.h"

CompositeAction::CompositeAction(const int n) : actions(0) {
  if (n>0) actions.reserve(n);
}

CompositeAction::~CompositeAction() {
  for (ActionIter a=actions.begin(); a<actions.end(); ++a) delete *a;
}

double CompositeAction::getActionDifference(
    const MultiLevelSampler& sampler, const int level) {
  double diff=0;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action) diff+=(*action)->getActionDifference(sampler,level);
  }
  return diff;
}

double CompositeAction::getActionDifference(const Paths &paths, 
    const VArray &displacement, int nmoving, const IArray &movingIndex, 
    int iFirstSlice, int nslice) {
  double diff=0;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action)
      diff+=(*action)->getActionDifference(paths,displacement,nmoving,
                                           movingIndex,iFirstSlice,nslice);
  }
  return diff;
}

double CompositeAction::getTotalAction(const Paths& paths, int level) const {
  double total=0;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action) total+=(*action)->getTotalAction(paths,level);
  }
  return total;
}

void CompositeAction::getBeadAction(const Paths& paths, int ipart, int islice,
       double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    double ui=0, utaui=0, ulambdai=0; 
    Vec fmi=0.0, fpi=0.0; 
    if (*action) (*action)->getBeadAction(paths,ipart,islice,
                                          ui,utaui,ulambdai,fmi,fpi);
    u+=ui; utau+=utaui; ulambda+=ulambdai; fm+=fmi; fp+=fpi;
  }
}

void CompositeAction::initialize(const SectionChooser& sectionChooser) {
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    (*action)->initialize(sectionChooser);
  }
}

void CompositeAction::acceptLastMove() {
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    (*action)->acceptLastMove();
  }
}

const Action* CompositeAction::getAction(const std::type_info &type) const {
  const Action* actionPtr=0;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    actionPtr=(*action)->getAction(type);
    if (actionPtr) return actionPtr;
  }
  return 0;
}
