#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PeriodicGaussian.h"
#include <complex>

PeriodicGaussian::PeriodicGaussian(double alpha, double length,
        int gridCount)
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

double PeriodicGaussian::operator()(double x) const {
    evaluate(x);
    return value;
}

double PeriodicGaussian::grad(double x) const {
    evaluate(x);
    return gradValue;
}


double PeriodicGaussian::d2(double x) const {
    evaluate(x);
    return d2Value;
}

void PeriodicGaussian::evaluate(double x) const {
    std::complex<double> w = exp(std::complex<double>(0, k * x));
    std::complex<double> sum0 = 0;
    std::complex<double> sum1 = 0;
    std::complex<double> sum2 = 0;
    std::complex<double> ik = std::complex<double>(0, k);
    double k2 = -k * k;
    for(int n = 1; n < nmax; ++n){
        std::complex<double> term = pow(nome, n*n) * pow(w, n);
        sum0 += term;
        sum1 += term * ik * double(n);
        sum2 += term * k2 * double(n * n) ;
    }
    value = prefactor * (1 + 2 * real(sum0));
    gradValue = prefactor * 2 * real(sum1);
    d2Value = prefactor * 2 * real(sum2);
}

const double PeriodicGaussian::PI = 3.141592653589793;

