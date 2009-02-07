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
#ifndef __CompositeAlgorithm_h_
#define __CompositeAlgorithm_h_

#include "Algorithm.h"
#include <vector>

/** Algorithm class combining several steps.
 * @version $Revision$
 * @author John Shumway */
class CompositeAlgorithm : public Algorithm {
public:
  /// Constructor.
  CompositeAlgorithm(const int nstep) : step(nstep,(Algorithm*)0){}
  /// Virtual destructor.
  virtual ~CompositeAlgorithm() {
    for(unsigned int i=0;i<step.size();++i) delete step[i];
  }
  /// Run the algorithm.
  virtual void run() {
    for(unsigned int i=0;i<step.size();++i) if (step[i]) step[i]->run();
  }
  /// Resize the number of steps in the algorithm.
  void resize(const int n) {
    for(unsigned int i=0;i<step.size();++i) delete step[i];
    step.resize(n);
  }
  /// Set algorithm for step i.
  void set(const int i, Algorithm* a) {step[i]=a;}
private:
  /// Steps in the algorithm.
  std::vector<Algorithm*> step;
};
#endif
