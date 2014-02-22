#ifndef PROPAGATORTESTUTIL_H_
#define PROPAGATORTESTUTIL_H_

class PropagatorTestUtil
{
public:
  /// Approximate grid propagator with continuum propagator.
  double approximateK0(double x1, double x2, double tau, double gridDeltaX)
  {
    const double PI = 3.141592653589793;
    double delta2 = (x1 - x2) * (x1 - x2);
    double prefactor = gridDeltaX / sqrt(2 * PI * tau / mass);
    return prefactor * exp(-0.5 * mass * delta2 / tau);
  }

  double K(double x1, double x2, double tau, double gridDeltaX)
  {
    const double PI = 3.141592653589793;
    double sinhwt = sinh(omega * tau);
    double coshwt = cosh(omega * tau);
    return gridDeltaX * sqrt(mass * omega / (2.0 * PI * sinhwt))
        * exp(
            -(mass * omega * (x1 * x1 + x2 * x2) * coshwt - 2 * x1 * x2)
                / (2.0 * sinhwt));
  }

  double mass;
  double omega;
};

#endif
