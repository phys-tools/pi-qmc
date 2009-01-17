// $Id: Demo.h,v 1.2 2006/12/29 21:55:49 jshumwa Exp $
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
#ifndef __Demo_h_
#define __Demo_h_

#include <iostream>

/// Base class for demonstration input-file generators.
/// @version $Revision: 1.2 $
/// @author John Shumway.
class Demo {
public:
  /// Construct a demo by name.
  static Demo* getDemo(const std::string& name);
  /// List demos.
  static void listDemos(std::ostream& out);
  /// Method to generate demo.
  virtual void generate() const=0;
};
#endif
