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
#include "CompositeDoubleAction.h"

CompositeDoubleAction::CompositeDoubleAction(const int n) : actions(0) {
  if (n>0) actions.reserve(n);
}

CompositeDoubleAction::~CompositeDoubleAction() {
  for (ActionIter a=actions.begin(); a<actions.end(); ++a) delete *a;
}

double CompositeDoubleAction::getActionDifference(
    const DoubleMLSampler& sampler, const int level) {
  double diff=0;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action) diff+=(*action)->getActionDifference(sampler,level);
  }
  return diff;
}

double CompositeDoubleAction::getTotalAction(const Paths& paths, int level) const {
  double total=0;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    if (*action) total+=(*action)->getTotalAction(paths,level);
  }
  return total;
}

void CompositeDoubleAction::getBeadAction(const Paths& paths, 
       int ipart, int islice,
       double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    double ui=0, utaui=0, ulambdai=0; 
    Vec fmi=0., fpi=0.;
    if (*action) (*action)->getBeadAction(paths,ipart,islice,
                                          ui,utaui,ulambdai,fmi,fpi);
    u+=ui; utau+=utaui; ulambda+=ulambdai; fm+=fmi; fp+=fpi;
  }
}

void CompositeDoubleAction::initialize(const DoubleSectionChooser& 
                                             sectionChooser) {
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    (*action)->initialize(sectionChooser);
  }
}

void CompositeDoubleAction::acceptLastMove() {
  for (ConstActionIter action=actions.begin(); action<actions.end(); ++action) {
    (*action)->acceptLastMove();
  }
}
