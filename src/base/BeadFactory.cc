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
#include "BeadFactory.h"

BeadFactory::BeadFactory() : protoBeads(*new Beads<NDIM>(1,1)) {
  beadNames.push_back("main");
}

BeadFactory::~BeadFactory() {
  delete &protoBeads;
}

Beads<NDIM>* BeadFactory::getNewBeads(const int npart, const int nslice) const {
  return dynamic_cast<Beads<NDIM>*>(protoBeads.cloneBase(npart,nslice));
}

int BeadFactory::addAuxBeads(BeadsBase* auxBeads,
                             const std::string& name) {
  auxBeads->auxBeads.pop_back();
  beadNames.push_back(name);
  int count = beadNames.size();
  // Push pointers to previous beads onto new auxBeads list.
  for (int i=0; i<count-1; ++i) {
    auxBeads->auxBeads.push_back(protoBeads.auxBeads[i]);
  }
  // Push pointer to new beads onto all beads' auxbeads list.
  for (int i=0; i<count; ++i) {
    protoBeads.auxBeads[i]->auxBeads.push_back(auxBeads);
    protoBeads.auxBeads[i]->nAuxBeads=count;
  }
  return count;
};
