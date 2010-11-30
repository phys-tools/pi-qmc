// $Id$
/*  Copyright (C) 2008-2009 John B. Shumway, Jr.

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
#include "TradEwaldSum.h"
#include <cmath>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Paths.h"
#include "Beads.h"
#include "MultiLevelSampler.h"

TradEwaldSum::TradEwaldSum(const SuperCell& cell, const int npart,
			   const double rcut, const double kcut, const double kappa)
  : EwaldSum(cell,npart,rcut,kcut), kappa(kappa) {
  setLongRangeArray();
  evalSelfEnergy();
}

TradEwaldSum::~TradEwaldSum() {
}

double TradEwaldSum::evalFR(const double r) const {
 return erf(kappa*r)/r;
}

double TradEwaldSum::evalFR0() const {
 return 2*kappa/sqrt(PI);
}

double TradEwaldSum::evalFK(const double k) const {
 return 4*PI*exp(-(k*k)/(4*kappa*kappa))/(k*k);
}

double TradEwaldSum::evalFK0() const {
 return -PI/(kappa*kappa);
}

