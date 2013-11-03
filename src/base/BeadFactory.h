#ifndef __BeadFactory_h_
#define __BeadFactory_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <iostream>
#include <vector>
#include <string>
#include "Beads.h"

/// Factory for Beads.
/// Uses cloneBase method of Beads to produce specialized Beads for the 
/// simulation.
/// @version $Revision$
/// @author John Shumway 
class BeadFactory {
public:
  /// Constructor.
  BeadFactory();
  /// Copy constructor.
  BeadFactory(const BeadFactory&);
  /// Destructor
  ~BeadFactory();
  /// Get new beads.
  Beads<NDIM>* getNewBeads(const int npart, const int nslice) const;
  /// Add prototype beads to factory.
  int addAuxBeads(BeadsBase*, const std::string& name);
private:
  /// Array of bead names.
  std::vector<std::string> beadNames;
  /// Prototype beads.
  Beads<NDIM>& protoBeads;
};
#endif
