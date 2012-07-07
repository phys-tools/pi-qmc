#include "CoulombLinkAction.h"
#include "Coulomb1DLinkAction.h"
#include "Coulomb3DLinkAction.h"
#include <cmath>
#include <blitz/tinyvec-et.h>

CoulombLinkAction::CoulombLinkAction(double q1q2, double epsilon, double mu,
		double deltaTau, int norder) :
		q1q2(q1q2),
		epsilon(epsilon),
		mu(mu),
		deltaTau(deltaTau),
		stau(q1q2 / epsilon * sqrt(2.0 * mu * deltaTau)),
		norder(norder) {
}

CoulombLinkAction::~CoulombLinkAction() {
}

double CoulombLinkAction::getValue(Vec delta1, Vec delta2) const {
	double r = calculateAverageSeparation(delta1, delta2);
	double s2 = calculateS2(delta1, delta2) / (r * r + 1e-200);
	Coulomb1DLinkAction coulomb1D(stau);
	Coulomb3DLinkAction coulomb3D(coulomb1D);
	double reff = 2.0 * mu * q1q2 * r / epsilon;
    double u = coulomb3D.calculateU0(reff);
    switch (norder)  {
    case 1:
        u += coulomb3D.calculateU1(reff) * s2;
        break;
    case 2:
        u += coulomb3D.calculateU1(reff) * s2;
        u += coulomb3D.calculateU2(reff) * s2 * s2;
        break;
    case 3:
        u += coulomb3D.calculateU1(reff) * s2;
        u += coulomb3D.calculateU2(reff) * s2 * s2;
        u += coulomb3D.calculateU3(reff) * s2 * s2 * s2;
        break;
    case 4:
        u += coulomb3D.calculateU1(reff) * s2;
        u += coulomb3D.calculateU2(reff) * s2 * s2;
        u += coulomb3D.calculateU3(reff) * s2 * s2 * s2;
        u += coulomb3D.calculateU4(reff) * s2 * s2 * s2 * s2;
        break;
    }
    return u;
}

double CoulombLinkAction::calculateAverageSeparation(Vec delta1, Vec delta2) {
	return 0.5 * (sqrt(dot(delta1, delta1)) + sqrt(dot(delta2, delta2)));
}

double CoulombLinkAction::calculateS2(Vec delta1, Vec delta2) {
	Vec s = delta1 - delta2;
	return dot(s, s);
}

