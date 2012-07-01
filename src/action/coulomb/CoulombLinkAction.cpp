#include "CoulombLinkAction.h"
#include "Coulomb1DLinkAction.h"
#include <cmath>

CoulombLinkAction::CoulombLinkAction(double q1q2, double epsilon,
		double deltaTau) {
	double mu = 1.0;
	stau = q1q2 / epsilon * sqrt(2.0 * mu * deltaTau);

}

CoulombLinkAction::~CoulombLinkAction() {
}

double CoulombLinkAction::getValue(Vec delta1, Vec delta2) const {
	return Coulomb1DLinkAction::calculate1DValueAtOrigin(stau);
}
