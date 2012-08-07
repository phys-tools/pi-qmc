#ifndef ERF_H_
#define ERF_H_

/* Incomplete gamma function
   1 / Gamma(a) * Int_0^x exp(-t) t^(a-1) dt  */
double p_gamma(double a, double x, double loggamma_a);

/* Incomplete gamma function
   1 / Gamma(a) * Int_x^inf exp(-t) t^(a-1) dt  */
double q_gamma(double a, double x, double loggamma_a);

double erf(double x);

double erfc(double x);

#endif
