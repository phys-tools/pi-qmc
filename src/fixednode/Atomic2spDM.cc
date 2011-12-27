#include "Atomic2spDM.h"

#include <blitz/tinyvec-et.h>

const double Atomic2spDM::PI = acos(-1.0);

Atomic2spDM::Atomic2spDM(double Z, int ifirst,
        int npart, int nfermion, double pweight, double weight)
:   AtomicOrbitalDM(ifirst, npart, nfermion, weight),
    nfermion(nfermion), coef(sqrt(weight*Z*Z*Z*Z*Z/(32*PI))),
    work1(nfermion,npart), work2(nfermion,npart),
    work3(nfermion,npart), work4(nfermion,npart),
    work5(nfermion,npart), work6(nfermion,npart),
    pweight(pweight), Z(Z) {
}


void
Atomic2spDM::evaluateValue(Matrix &mat, double scale) const {
    // First set work1 and work2 to r and r'.
    for (int i=0; i<nfermion; ++i) {
        for (int k=0; k<npart; ++k) {
            work1(i,k) = sqrt(dot(dist1(i,k),dist1(i,k)));
            work2(i,k) = sqrt(dot(dist2(i,k),dist2(i,k)));
        }
    }
    // Next set unit vectors arrays work3 and work4.
    for (int i=0; i<nfermion; ++i) {
        for (int k=0; k<npart; ++k) {
            work3(i,k) = dist1(i,k)/(work1(i,k)+1e-300);
            work4(i,k) = dist2(i,k)/(work2(i,k)+1e-300);
        }
    }
    // Next set work1 and work2 to psi and psi*.
    for (int i=0; i<nfermion; ++i) {
        for (int k=0; k<npart; ++k) {
            work1(i,k) *= coef*exp(-0.5*Z*work1(i,k));
            work2(i,k) *= coef*exp(-0.5*Z*work2(i,k));
        }
    }
    for (int i=0; i<nfermion; ++i) {
        for (int j=0; j<nfermion; ++j) {
            for (int k=0; k<npart; ++k) {
                mat(i,j) += scale*work1(j,k)*work2(i,k)
                            *(1.+pweight*dot(work3(j,k),work4(i,k)));
            }
        }
    }
}

void
Atomic2spDM::evaluateValueAndGrad() const {
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
    // Next set work3 and work4 to unit vectors terms.
    for (int i=0; i<nfermion; ++i) {
        for (int k=0; k<npart; ++k) {
            work3(i,k) = dist1(i,k)*(1./(work1(i,k)+1e-300) - 0.5*Z);
            work4(i,k) = dist2(i,k)*(1./(work2(i,k)+1e-300) - 0.5*Z);
        }
    }
    // Then evaluate the orbitals.
    for (int i=0; i<nfermion; ++i) {
        for (int k=0; k<npart; ++k) {
            work5(i,k) = work1(i,k)*coef*exp(-Z*0.5*work1(i,k));
            work6(i,k) = work2(i,k)*coef*exp(-Z*0.5*work2(i,k));
        }
    }
}

AtomicOrbitalDM::ValAndGrad
Atomic2spDM::operator()(int i, int j) const {
    ValAndGrad result;
    result.val = 0.;
    result.grad1 = 0.;
    result.grad2 = 0.;
    for (int k=0; k<npart; ++k) {
        double dotprod=dot(dist1(j,k),dist2(i,k));
        double angle=(1.+pweight*dotprod);
        result.val += work5(j,k)*work6(i,k)*angle;
        result.grad1 += (pweight*(dist2(i,k) -
                dotprod*dist1(j,k))/work1(j,k) +
                work3(j,k)*angle)*work5(j,k)*work6(i,k);
        result.grad2 += (pweight*(dist1(j,k) -
                dotprod*dist2(i,k))/work2(i,k) +
                work4(i,k)*angle)*work5(j,k)*work6(j,k);
    }
    return result;
}
