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
#include "XMLParser.h"
#include <blitz/tinyvec-et.h>

double XMLParser::getDoubleAttribute(const xmlNodePtr& node,
                                     const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  double value = temp?atof((char*)temp):0;
  xmlFree(temp);
  return value;
}

int XMLParser::getIntAttribute(const xmlNodePtr& node,
                               const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  int value = temp?atoi((char*)temp):0;
  xmlFree(temp);
  return value;
}

std::string XMLParser::getName(const xmlNodePtr& node) {
  std::string value((char*)node->name);
  return value;
}

bool XMLParser::getBoolAttribute(const xmlNodePtr& node,
                                          const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string value(temp?(char*)temp:"");
  xmlFree(temp);
  return value!="" && (value[0]=='t' || value[0]=='T');
}

std::string XMLParser::getStringAttribute(const xmlNodePtr& node,
                                          const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string value(temp?(char*)temp:"");
  xmlFree(temp);
  return value;
}

const blitz::TinyVector<double,NDIM> 
XMLParser::getVecAttribute(const xmlNodePtr& node) {
  blitz::TinyVector<double,NDIM> v;
  v[0]=getDoubleAttribute(node,"x");
  v=v[0];
#if NDIM>1
  v[1]=getDoubleAttribute(node,"y");
#endif
#if NDIM>2
  v[2]=getDoubleAttribute(node,"z");
#endif
#if NDIM>3
  v[3]=getDoubleAttribute(node,"k");
#endif
  return v;
}

std::string XMLParser::getLinkBase(const xmlNodePtr& node,
                                   const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string href(temp?(char*)temp:"");
  xmlFree(temp);
  return href.substr(0,href.find('|'));
}

std::string XMLParser::getLinkPath(const xmlNodePtr& node,
                                   const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string href(temp?(char*)temp:"");
  xmlFree(temp);
  return href.substr(href.find('|')+1);
}
