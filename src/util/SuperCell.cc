#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "SuperCell.h"
#include <iostream>
#include <blitz/tinyvec-et.h>

SuperCell::SuperCell(const Vec a) : a(a) {
    computeRecipricalVectors();
}

SuperCell::~SuperCell() {}

void SuperCell::computeRecipricalVectors() {
  b=1./a;
  /// Set the rcut2 to be the square of half the shortest side.
  double l2=a[0]*a[0]; rcut2=l2;
  for (int i=1; i<NDIM; ++i)  {
    l2=a[i]*a[i]; 
    rcut2=(rcut2>l2)?l2:rcut2;
  }
  rcut2 *= 0.25;
}

SuperCell::Vec& SuperCell::pbc(Vec& v) const {
  if (dot(v,v)>rcut2) {
    for (int idim=0; idim<NDIM; ++idim) {
      // int i = (int)(b[idim]*v[idim]+1000.5)-1000;
      double i = trunc (b[idim]*v[idim]+1000.5)-1000; 
      if (i!=0) v[idim]-=a[idim]*i;
    }
  }
  return v;
}

double SuperCell::pbc(double dist, int idim) const {
  if (dist>(0.5*a[idim])) {
    double i = trunc (b[idim]*dist+1000.5)-1000; 
    if (i!=0) dist-=a[idim]*i;
  }
  return dist;
}



std::ostream& SuperCell::write(std::ostream& os) const {
  os << "SuperCell:" << std::endl;
  for (int i=0; i<NDIM; ++i) {
    os << "    (" ;
    for (int j=0; j<i; ++j) os << "0,";
    os << a[i]; 
    for (int j=i+1; j<NDIM; ++j) os << ",0";
    os << ")" << std::endl;
  }
  return os;
}

std::ostream& operator<<(std::ostream &os, const SuperCell &cell) {
    cell.write(os);
    return os;
}

