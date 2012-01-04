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
#ifndef __XMLUnitParser_h_
#define __XMLUnitParser_h_
#include "parser/XMLParser.h"
#include <string>
class Units;

/** Base class for classes that parse XML data with units.  
  * Includes helper functions.
  * @version $Revision$
  * @author John Shumway */
class XMLUnitParser : public XMLParser{
public:
  /// Constructor.
  XMLUnitParser(const Units*);
  /// Destructor. 
  virtual ~XMLUnitParser() {};
protected:
  /// Information about the units and conversions.
  const Units* units;
  /// Parse a length with unit conversion.
  double getEnergyAttribute(const xmlNodePtr &node,
                            const std::string &attName);
  /// Parse a energy with unit conversion.
  double getLengthAttribute(const xmlNodePtr &node,
                            const std::string &attName);
  /// Parse an inverse energy with unit conversion.
  double getTimeAttribute(const xmlNodePtr &node,
                          const std::string &attName);
  /// Parse an inverse length with unit conversion.
  double getInvLengthAttribute(const xmlNodePtr &node,
                            const std::string &attName);
  /// Parse a density with unit conversion.
  double getDensityAttribute(const xmlNodePtr &node,
                             const std::string &attName);
  /// Parse a mass with unit conversion.
  /// Parse a mass with unit conversion.
  double getMassAttribute(const xmlNodePtr &node,
                          const std::string &attName);
  /// Parse a field strength with unit conversion.
  double getFieldStrengthAttribute(const xmlNodePtr &node,
                                   const std::string &attName);
  /// Internal method to parse a unit and value.
  void parseUnitAndValue(const xmlNodePtr &node, const std::string &attName);
  /// The most recently parsed value.
  double value;
  /// The most recently parsed unit name.
  std::string unit;
};
#endif
