#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PeriodicGaussian.h"
#include <complex>

PeriodicGaussian::PeriodicGaussian(double alpha, double length)
:   alpha(alpha),
    length(length),
    nome(exp(-PI*PI/(alpha*length*length))),
    nmax(numberOfTerms(alpha, length)),
    prefactor(sqrt(PI/(alpha*length*length))),
    k(2 * PI / length) {
}

int PeriodicGaussian::numberOfTerms(double alpha, double length) {
    return 1 + int(2.0 * sqrt(alpha*length*length));
}

double PeriodicGaussian::evaluate(double x) const {
    std::complex<double> w = exp(std::complex<double>(0, k * x));
    std::complex<double> sum0 = 0;
    std::complex<double> sum1 = 0;
    std::complex<double> sum2 = 0;
    std::complex<double> ik = std::complex<double>(0, k);
    double k2 = -k * k;
    std::complex<double> term = 1.0;
    double nomeN = 1.0;
    for (int n = 1; n < nmax; ++n) {
        term *= nomeN * nomeN * nome * w;
        sum0 += term;
        sum1 += term * ik * double(n);
        sum2 += term * k2 * double(n * n);
        nomeN *= nome;
    }
    value = prefactor * (1 + 2 * real(sum0));
    gradient = prefactor * 2 * real(sum1);
    secondDerivative = prefactor * 2 * real(sum2);
    return value;
}

const double PeriodicGaussian::PI = 3.141592653589793;

