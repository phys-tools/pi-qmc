// $Id: Permutation.cc,v 1.4 2006/10/18 17:08:19 jshumwa Exp $
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
#include "Permutation.h"

Permutation::Permutation(const int n) : permutation(n) {
  for (int i=0; i<n; ++i) permutation(i)=i;
}

Permutation::Permutation(const Permutation& p)
  : permutation(p.permutation.copy()) {
}

Permutation& Permutation::operator=(const Permutation& p) {
  permutation=p.permutation;
  return *this;
}

void Permutation::reset() {
  for (int i=0; i<permutation.size(); ++i) permutation(i)=i;
}

Permutation& Permutation::prepend(const Permutation& p) {
  int temp[permutation.size()];
  for (int i=0; i<permutation.size(); ++i) {
    temp[i]=permutation(p.permutation(i));
  }
  for (int i=0; i<permutation.size(); ++i) permutation(i)=temp[i];
  return *this;
}

Permutation& Permutation::append(const Permutation& p) {
  int temp[permutation.size()];
  for (int i=0; i<permutation.size(); ++i) {
    temp[i]=p.permutation(permutation(i));
  }
  for (int i=0; i<permutation.size(); ++i) permutation(i)=temp[i];
  return *this;
}

void Permutation::setToInverse(const Permutation& p) {
  for (int i=0; i<permutation.size(); ++i) {
    permutation(p.permutation(i)) = i;
  }
}

bool Permutation::isIdentity() const {
  for (int i=0; i<permutation.size(); ++i) {
    if (permutation(i)!=i) return false;
  }
  return true;
}


std::ostream& operator<<(std::ostream &os, const Permutation &p) {
  return os << p.permutation;
}
