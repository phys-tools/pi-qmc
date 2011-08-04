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
#include "Species.h"
#include <cstdlib>
#include <blitz/tinyvec.h>

Species::Species() 
  : name("?"), count(0), ifirst(0), mass(1.0), charge(0.0), 
    twos(0), isFermion(false), anMass(0) {
}

Species::Species(const Species& s)  
  : name(s.name), count(s.count), ifirst(s.ifirst), mass(s.mass), 
    charge(s.charge), twos(s.twos),
    isFermion(s.isFermion), 
    anMass((s.anMass)?new Vec(*s.anMass):0) {
}

Species& Species::operator=(const Species& s) {
  name=s.name; count=s.count; ifirst=s.ifirst; mass=s.mass; charge=s.charge;
  twos=s.twos;
  isFermion=s.isFermion; 
  anMass=(s.anMass)?new Vec(*s.anMass):0;
  return *this;
}

Species::~Species() {delete anMass;}

Species::Species(const std::string& name, const int count, const double mass,
  const double charge, const int twos, const bool isFermion) 
  : name(name), count(count), ifirst(ifirst), mass(mass), charge(charge),
    twos(twos), isFermion(isFermion) {
}

std::ostream& operator<<(std::ostream& os, const Species& s) {
  os << s.count << " " << s.name << ", s=" 
     << (s.twos%2==0?s.twos/2:s.twos) 
     << (s.twos%2==0?"":"/2") << ", m=";
  if (s.anMass) os << *s.anMass; else os << s.mass;
  os << ", q=" << s.charge;
  if (fabs(s.displace)>1e-10) os << ", displace=" << s.displace;
  os << (s.isFermion?" (fermionic)":"") 
     << std::endl;
  return os;
}
