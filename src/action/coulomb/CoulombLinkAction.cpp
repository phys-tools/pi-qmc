#include "CoulombLinkAction.h"
#include "Coulomb1DLinkAction.h"
#include <cmath>
#include <blitz/tinyvec-et.h>

CoulombLinkAction::CoulombLinkAction(double q1q2, double epsilon, double mu,
		double deltaTau) :
		q1q2(q1q2),
		epsilon(epsilon),
		mu(mu),
		deltaTau(deltaTau),
		stau(q1q2 / epsilon * sqrt(2.0 * mu * deltaTau)) {
}

CoulombLinkAction::~CoulombLinkAction() {
}

double CoulombLinkAction::getValue(Vec delta1, Vec delta2) const {
	double r = calculateAverageSeparation(delta1, delta2);
	double s2 = calculateS2(delta1, delta2);
	double u0AtOrigin = Coulomb1DLinkAction::calculateValueAtOrigin(stau);
	double taueff = 2.0 * mu * q1q2 * q1q2 * deltaTau / (epsilon * epsilon);
	double reff = 2.0 * mu * q1q2 * r / epsilon;
	return Coulomb1DLinkAction::calculateU0(stau, u0AtOrigin, reff, taueff);
}

double CoulombLinkAction::calculateAverageSeparation(Vec delta1, Vec delta2) {
	return 0.5 * (sqrt(dot(delta1, delta1)) + sqrt(dot(delta2, delta2)));
}

double CoulombLinkAction::calculateS2(Vec delta1, Vec delta2) {
	Vec s = delta1 - delta2;
	return dot(s, s);
}

