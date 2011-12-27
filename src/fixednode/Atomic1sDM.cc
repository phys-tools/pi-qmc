#include "Atomic1sDM.h"

#include <blitz/tinyvec-et.h>

Atomic1sDM::Atomic1sDM(
  double Z, int ifirst, int npart, int nfermion, double weight)
  : AtomicOrbitalDM(ifirst, npart, nfermion, weight),
    nfermion(nfermion), coef(sqrt(weight*Z*Z*Z/PI)),
    work1(nfermion,npart), work2(nfermion,npart), Z(Z) {
}

Atomic1sDM::~Atomic1sDM() {
}

void
Atomic1sDM::evaluateValue(Matrix &mat, double scale) const {
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      work1(i,k) = sqrt(dot(dist1(i,k),dist1(i,k)));
      work2(i,k) = sqrt(dot(dist2(i,k),dist2(i,k)));
    }
  }
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      work1(i,k) = coef*exp(-Z*work1(i,k));
      work2(i,k) = coef*exp(-Z*work2(i,k));
    }
  }
  for (int i=0; i<nfermion; ++i) {
    for (int j=0; j<nfermion; ++j) {
      for (int k=0; k<npart; ++k) {
        mat(i,j) += scale*work1(j,k)*work2(i,k);
      }
    }
  }
}


void
Atomic1sDM::evaluateValueAndGrad() const {
  // First set the work arrays to the scalar distances.
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      work1(i,k) = sqrt(dot(dist1(i,k),dist1(i,k)));
      work2(i,k) = sqrt(dot(dist2(i,k),dist2(i,k)));
    }
  }
  // Then rescale the distance arrays to their unit vectors.
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      dist1(i,k) /= (work1(i,k)+1e-300);
      dist2(i,k) /= (work2(i,k)+1e-300);
    }
  }
  // Then evaluate the orbitals.
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      work1(i,k) = coef*exp(-Z*work1(i,k));
      work2(i,k) = coef*exp(-Z*work2(i,k));
    }
  }
  // Finally, set the distance arrays to the gradients.
  for (int i=0; i<nfermion; ++i) {
    for (int k=0; k<npart; ++k) {
      dist1(i,k) *= -Z*work1(i,k);
      dist2(i,k) *= -Z*work2(i,k);
    }
  }
}

AtomicOrbitalDM::ValAndGrad
Atomic1sDM::operator()(int i, int j) const {
  ValAndGrad result;
  result.val = 0.;
  result.grad1 = 0.;
  result.grad2 = 0.;
  for (int k=0; k<npart; ++k) {
    result.val += work1(j,k)*work2(i,k);
    result.grad1 += dist1(j,k)*work2(i,k);
    result.grad2 += work1(j,k)*dist2(i,k);
  }
  return result;
}


const double Atomic1sDM::PI = acos(-1.0);
