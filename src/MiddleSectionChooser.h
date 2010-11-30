// $Id$
/*  Copyright (C) 2010 John B. Shumway, Jr.

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
#ifndef __MiddleSectionChooser_h_
#define __MiddleSectionChooser_h_

#include "SectionChooser.h"
#include <gsl/gsl_qrng.h>

/// Algorithm class for choosing a section.
/// @version $Revision$
/// @author John Shumway
class MiddleSectionChooser : public SectionChooser {
public:
  /// Constructor.
  MiddleSectionChooser(const int nlevel, Paths &paths, Action& action,
                 const BeadFactory&);
  /// Virtual destructor.
  virtual ~MiddleSectionChooser();
  /// Run the algorithm.
  virtual void run();
};
#endif
