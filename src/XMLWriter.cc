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
#include "XMLWriter.h"
#include <libxml/parser.h>
#include <sstream>

void XMLWriter::dumpXMLToFile(const std::string& filename) const {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
  doc->children = xmlNewDocNode(doc,0,BAD_CAST filename.c_str(),0);
  addToXMLTree(doc->children);
  xmlSaveFormatFile((filename+".xml").c_str(),doc,1);
}
