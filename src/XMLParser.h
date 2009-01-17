// $Id: XMLParser.h,v 1.5 2006/10/18 17:08:19 jshumwa Exp $
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
#ifndef __XMLParser_h_
#define __XMLParser_h_
#include <string>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <blitz/tinyvec.h>

/** Base class for classes that parse XML data.  
  * Includes helper functions.
  * @version $Revision: 1.5 $
  * @author John Shumway */
class XMLParser {
public:
/** Destructor. */
virtual ~XMLParser() {};
/** Parse some xml. */
virtual void parse(const xmlXPathContextPtr& ctxt)=0;
protected:
/** Get the name of the node.*/
static std::string getName(const xmlNodePtr& node);
/** Get a double valued attribute.*/
static double getDoubleAttribute(const xmlNodePtr& node,
                          const std::string& attName);
/** Get an integer valued attribute.*/
static int getIntAttribute(const xmlNodePtr& node,
                       const std::string& attName);
/** Get a string valued attribute.*/
static std::string getStringAttribute(const xmlNodePtr& node,
                               const std::string& attName);
/** Get a bool valued attribute.*/
static bool getBoolAttribute(const xmlNodePtr& node,
                               const std::string& attName);
/** Get a TinyVector valued attribute.*/
static const blitz::TinyVector<double,NDIM> 
  getVecAttribute(const xmlNodePtr& node);

/** Get the base of a link.*/
static std::string getLinkBase(const xmlNodePtr& node,
                               const std::string& attName="href");

/** Get the path of a link.*/
static std::string getLinkPath(const xmlNodePtr& node,
                               const std::string& attName="href");
};
#endif
