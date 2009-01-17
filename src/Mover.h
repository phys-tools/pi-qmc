// $Id: Mover.h,v 1.3 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __Mover_h_
#define __Mover_h_
class MultiLevelSampler;
/** Virtual base class for routines to select trial moves for beads.
  * @version $Revision: 1.3 $
  * @author John Shumway. */
class Mover{
public:
  /// Virtual destructor.
  virtual ~Mover() {}
  /// Move the samplers moving beads for a given level, returning
  /// the probability for the old move divided by the probability for the 
  /// new move.
  virtual double makeMove(MultiLevelSampler&, const int level)=0;
};
#endif
