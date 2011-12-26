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
#include "NodeModel.h"

NodeModel::NodeModel(const std::string& name) :
#ifdef NODE_DIST_DEBUG
   logFile(std::string("logGradDist"+name+".test").c_str())
#endif
   updateObj(0), spinModelState(0) {
}

void NodeModel::testGradLogDist( const VArray &r1, const VArray &r2,
  const int islice, VMatrix &gradd1, VMatrix &gradd2,
  const Array &d1, const Array &d2, const int npart, const int ifirst) {
#ifdef NODE_DIST_DEBUG
  logFile << d1 << d2 << std::endl;
#endif
}
