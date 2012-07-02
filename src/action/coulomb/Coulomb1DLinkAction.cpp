#include "Coulomb1DLinkAction.h"

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

double Coulomb1DLinkAction::calculateU0(double u0, double reff,
		double taueff) const {

	double a = 0.25300593 * stauToMinus1 + 0.01432126;
	double b = 0.07936898 * stauToMinus2 - 0.01634421 / stau;
	double c = 0.07263383 * stauToMinus3;
	double d = 0.00940013 * stauToMinus4;
	double e = 0.03160181 * stauToMinus5;
	double f = 0.18976335 * stauToMinus1;
	double g = 0.00343053 * stauToMinus2;

	double u1_0 = (u0
			+ reff
					* ((u0 * a - 1.)
							+ reff * (f + reff * (g + reff * e * taueff))));
	u1_0 /= (1.
			+ reff * (a + reff * (b + reff * (c + reff * (d + reff * e)))));
	return u1_0;
}
