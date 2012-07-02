#include "Coulomb1DLinkAction.h"
#include <cmath>

Coulomb1DLinkAction::Coulomb1DLinkAction(double stau) :
		stau(stau),
		stauToMinus1(1.0 / stau),
		stauToMinus2(stauToMinus1 * stauToMinus1),
		stauToMinus3(stauToMinus2 * stauToMinus1),
		stauToMinus4(stauToMinus3 * stauToMinus1),
		stauToMinus5(stauToMinus4 * stauToMinus1) {
}

Coulomb1DLinkAction::~Coulomb1DLinkAction() {
}

double Coulomb1DLinkAction::calculateValueAtOrigin() const {
	double u0 =
			stau * (1.772453851
					+ stau * (-0.074137740
							+ stau * (0.005834805
									+ stau * (-0.000382686
											+ stau * (0.000008738
													+ stau * 0.000002138)))));
	return u0;
}

double Coulomb1DLinkAction::calculateU0(double uOrigin, double reff,
		double taueff) const {

	double a = 0.25300593 * stauToMinus1 + 0.01432126;
	double b = 0.07936898 * stauToMinus2 - 0.01634421 / stau;
	double c = 0.07263383 * stauToMinus3;
	double d = 0.00940013 * stauToMinus4;
	double e = 0.03160181 * stauToMinus5;
	double f = 0.18976335 * stauToMinus1;
	double g = 0.00343053 * stauToMinus2;

	double u0 = (uOrigin + reff * ((u0 * a - 1.)
							+ reff * (f + reff * (g + reff * e * taueff))));
	u0 /= (1.
			+ reff * (a + reff * (b + reff * (c + reff * (d + reff * e)))));
	return u0;
}

double Coulomb1DLinkAction::calculateU1(double stau, double reff) const {
	double a = 0.07776906 * stauToMinus1;
	double b = -0.04695226 * stauToMinus2;
	double c = 0.01291994 * stauToMinus3;
	double d = 0.60942687 * stauToMinus1;
	double e = -1.09588064 * stauToMinus2;
	double f = 1.15591565 * stauToMinus3;
	double g = -0.57485278 * stauToMinus4;
	double h = 0.15212765 * stauToMinus5;
	double u1 = reff * reff * (a + reff * (b + reff * c));
	u1 /= (1 + reff * (d + reff * (e + reff * (f + reff * (g + reff * h)))));
	return u1;
}

double Coulomb1DLinkAction::calculateU2(double stau, double reff) const {
	double a = 9.52130652e-04 * pow(stau, -3);
	double b = -6.39946071e-04 * pow(stau, -4);
	double c = 1.14359558e-04 * pow(stau, -5);
	double d = -7.59575215e-01 * pow(stau, -3);
	double e = 5.53798980e-01 * pow(stau, -4);
	double f = -7.15626226e-02 * pow(stau, -5);
	double g = -3.49081618e-02 * pow(stau, -6);
	double h = 8.40084881e-03 * pow(stau, -7);
	double u1_2 =
			pow(reff, 4) * (a + reff * (b + reff * c))
			/ (1
					+ pow(reff, 3)
					* (d
							+ reff
							* (e
									+ reff
									* (f
											+ reff
											* (g
													+ reff
													* h)))));
	return u1_2;
}

double Coulomb1DLinkAction::calculateU3(double stau, double reff) const {
	double a = -2.06984256e-1 * pow(stau, -5);
	double b = -9.51423947e-3 * pow(stau, -6);
	double c = 8.97532561e-2 * pow(stau, -7);
	double d = 6.10876797e+4 * pow(stau, -3);
	double e = -1.19016292e+5 * pow(stau, -4);
	double f = 9.29836599e+4 * pow(stau, -5);
	double g = -3.49919130e+4 * pow(stau, -6);
	double h = 6.42490539e+3 * pow(stau, -7);
	double i = -6.02729064e+2 * pow(stau, -8);
	double j = 6.01285550e+1 * pow(stau, -9);
	double u1_3 =
			pow(reff, 6) * (a + reff * (b + reff * c))
			/ (1
					+ pow(reff, 3)
					* (d
							+ reff
							* (e
									+ reff
									* (f
											+ reff
											* (g
													+ reff
													* (h
															+ reff
															* (i
																	+ reff
																	* j)))))));
	return u1_3;
}

double Coulomb1DLinkAction::calculateU4(double stau, double reff) const {
	double a = -2.84102296e-6 * pow(stau, -5);
	double b = 6.26961672e-7 * pow(stau, -7);
	double c = -2.20166473e-1 * pow(stau, -3);
	double d = 1.07981903e-1 * pow(stau, -5);
	double e = -1.94671111e-2 * pow(stau, -7);
	double f = 1.56217930e-3 * pow(stau, -9);
	double u1_4 =
			pow(reff, 6) * (a + b * reff * reff)
			/ (1
					+ reff * reff * reff
					* (c
							+ reff * reff
							* (d
									+ reff * reff
									* (e
											+ reff
											* reff
											* f))));
	return u1_4;
}


