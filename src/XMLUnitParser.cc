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
#include "XMLUnitParser.h"
#include "stats/Units.h"
#include <blitz/tinyvec-et.h>

XMLUnitParser::XMLUnitParser(const Units* units) : XMLParser(), units(units) {
};

double XMLUnitParser::getLengthAttribute(const xmlNodePtr &node,
                                         const std::string &attName) {
  parseUnitAndValue(node,attName);
  double scale=1;
  if (unit!="") scale=units->getLengthScaleIn(unit);
  return scale*value;
}

double XMLUnitParser::getEnergyAttribute(const xmlNodePtr &node,
                                         const std::string &attName) {
  parseUnitAndValue(node,attName);
  double scale=1;
  if (unit!="") scale=units->getEnergyScaleIn(unit);
  return scale*value;
}

double XMLUnitParser::getMassAttribute(const xmlNodePtr &node,
                                         const std::string &attName) {
  parseUnitAndValue(node,attName);
  double scale=1;
  if (unit!="") scale=units->getMassScaleIn(unit);
  return scale*value;
}

double XMLUnitParser::getTimeAttribute(const xmlNodePtr &node,
                                       const std::string &attName) {
  parseUnitAndValue(node,attName);
  double scale=1;
  if (unit!="") {
    // Assume time units end in "-1".
    scale=1.0/units->getEnergyScaleIn(unit.substr(0,unit.length()-2));
  }
  return scale*value;
}

double XMLUnitParser::getInvLengthAttribute(const xmlNodePtr &node,
                                       const std::string &attName) {
  parseUnitAndValue(node,attName);
  double scale=1;
  if (unit!="") {
    // Assume units end in "-1".
    scale=1.0/units->getLengthScaleIn(unit.substr(0,unit.length()-2));
  }
  return scale*value;
}

void XMLUnitParser::parseUnitAndValue(const xmlNodePtr &node,
                                        const std::string &attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  if (temp==0) {value=0; unit=""; return;}
  std::string attr(temp);
  std::string::size_type loc = attr.find( " ", 0 );
  if( loc != std::string::npos ) {
    value = atof(attr.substr(0,loc).c_str());
    unit = attr.substr(loc+1);
  } else {
    value = atof(attr.c_str());
    unit = "";
  }
}
