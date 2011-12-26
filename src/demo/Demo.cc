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
#include "Demo.h"
#include "SHODemo.h"
#include "HAtomDemo.h"
#include "HeAtomDemo.h"
#include "SpinDemo.h"
#include "WireDemo.h"

Demo* Demo::getDemo(const std::string& name) {
  Demo *demo=0;
  if (name=="SHO") {
    demo = new SHODemo();
  } else if (name=="H_atom") {
    demo = new HAtomDemo();
  } else if (name=="He_atom") {
    demo = new HeAtomDemo();
  } else if (name=="spin") {
    demo = new SpinDemo();
  } else if (name=="wire") {
    demo = new WireDemo();
  } else {
    std::cout << "Could not find demo: " << name << "\n" << std::endl;
  }
  return demo;
}

void Demo::listDemos(std::ostream& out) {
  out << "Available demos:" << std::endl;
  out << "     SHO                  A simple harmonic oscillator" << std::endl;
  out << "     H_atom               A hydrogen atom" << std::endl;
  out << "     He_atom              A helium atom" << std::endl;
  out << "     spin                 A single spin in a magnetic field" 
      << std::endl;
  out << "     wire                 10 electrons in a quasi-1D GaAs (2D) wire" 
      << std::endl;
  out << std::endl;
}
