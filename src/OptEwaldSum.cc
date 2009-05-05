// $Id: OptEwaldSum.cc 38 2009-04-09 20:01:17Z john.shumwayjr $
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#include "OptEwaldSum.h"
#include <cmath>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Paths.h"
#include "Beads.h"
#include "MultiLevelSampler.h"

OptEwaldSum::OptEwaldSum(const SuperCell& cell, int npart,
  double rcut, double kcut, double khalo, int npoly, int ncts)
  : EwaldSum(cell, npart, rcut, kcut),
    npoly(npoly), coef(npoly) {
}

OptEwaldSum::~OptEwaldSum() {
}

double OptEwaldSum::evalVShort(double r) const {
  double v=0.;
  double r2 = r*r;
  for (int i=npoly-1; i>=0; --i) {
    v *= r2;
    v += coef(i);
  }
  return v;
}

double OptEwaldSum::evalVLong(double k2) const {
  return 0;
}

void OptEwaldSum::evalSelfEnergy() const {
}
