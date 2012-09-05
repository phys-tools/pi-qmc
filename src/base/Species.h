#ifndef __Species_h_
#define __Species_h_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <iostream>
#include <cstdlib>
#include <blitz/tinyvec.h>

/** Class for storing species information.
 * @author John Shumway */
class Species {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  Species();
  Species(const std::string& name, const int count, const double mass,
          const double charge, const int twos, const bool isFermion);
  Species(const Species& s);
  Species& operator=(const Species& s);
  ~Species();

  /// The name of the species.
  std::string name;
  /// The number of particles.
  int count;
  /// The index of the first particle of this species.
  int ifirst;
  /// The mass of the species.
  double mass;
  /// The charge of the species.
  double charge;
  /// The spin of the species (in units of 1/2).
  int twos;
  /// Flag for whether the particles have fermi statistics.
  bool isFermion;
  /// Flag to indicate whether particle is static or dynamic.
  bool isStatic;
  /// Pointer to opitional anisotropic mass.
  Vec* anMass;
  /// An optional displacement (useful for parallel quantum wells in 2D).
  double displace;
};

std::ostream& operator<<(std::ostream& os, const Species& s);
#endif
