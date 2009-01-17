// $Id: MainParser.h,v 1.2 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __MainParser_h_
#define __MainParser_h_
#include "XMLParser.h"
/** Main parser for PIMC code. 
  * @version $Revision: 1.2 $
  * @author John Shumway */
class MainParser : public XMLParser {
public:
  /// Construct by providing an input xml file.
  MainParser(const std::string& filename);
  /// Virtual destructor.
  ~MainParser();
  /// Parse from the context.
  void parse();
  /// Parse from a context.
  virtual void parse(const xmlXPathContextPtr& ctxt);
private:
  /// The xml document to parse from.
  xmlDocPtr doc;
  /// The context to parse from.
  xmlXPathContextPtr context;
};
#endif
