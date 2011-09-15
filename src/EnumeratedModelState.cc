// $Id: EnumeratedModelState.cc 338 2010-11-30 18:56:16Z john.shumwayjr $
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
#include "EnumeratedModelState.h"
#include <cstdlib>

EnumeratedModelState::EnumeratedModelState(int modelCount)
  : modelCount(modelCount), modelState(0) {
}

void EnumeratedModelState::write(std::ostream &os) const {
  os << ", state " << modelState+1 << " of " << modelCount << ".";
}

bool EnumeratedModelState::read(const std::string &line) {
  std::cout << "Checking for model state..." << std::endl;
  int i = line.find("state");
  if (i != -1) {
    i += 6;
    int j = line.find("of");
    modelState = atoi(line.substr(i,j-i-1).c_str());
    //modelCount = atoi(line.substr(j+3).c_str());
    // Hack to handle old and new indexing of model states.
    // New convention is to start from 1 and end line with "."
    if (line[line.size()-1]=='.') --modelState;
    std::cout << "Model state " << modelState << "." << std::endl;
  }
}
