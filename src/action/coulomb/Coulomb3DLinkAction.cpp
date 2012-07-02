#include "Coulomb3DLinkAction.h"
#include "Coulomb1DLinkAction.h"
#include <cmath>
#include <iostream>

Coulomb3DLinkAction::Coulomb3DLinkAction(Coulomb1DLinkAction &coulomb1D) :
        coulomb1D(coulomb1D) {
}

Coulomb3DLinkAction::~Coulomb3DLinkAction() {
}

double Coulomb3DLinkAction::calculateU0(double taueff, double reff) const {
    double u0 = coulomb1D.calculateValueAtOrigin();
    double u1_0 = coulomb1D.calculateU0(u0, reff, taueff);
    double u1_1 = coulomb1D.calculateU1(reff);

    double a = 2. * taueff / (reff * reff + 1e-200);
    double b = exp(-2. / a);
    double c = 1. / (1 + 2 * a * (1 - b) * u1_1);
    double u = u1_0 + log(c);
    return u;
}

double Coulomb3DLinkAction::calculateU1(double taueff, double reff) const {
    double u1_1 = coulomb1D.calculateU1(reff);
    double u1_2 = coulomb1D.calculateU2(reff);

    double a = 2. * taueff / (reff * reff + 1e-200);
    double b = exp(-2. / a);
    double c = 1. / (1 + 2 * a * (1 - b) * u1_1);
    double u = u1_1 + c * (b * u1_1 - 4 * a * (1 - b) * u1_2);
    return u;
}

double Coulomb3DLinkAction::calculateU2(double taueff, double reff) const {
    double u1_1 = coulomb1D.calculateU1(reff);
    double u1_2 = coulomb1D.calculateU2(reff);
    double u1_3 = coulomb1D.calculateU3(reff);

    double a = 2. * taueff / (reff * reff + 1e-200);
    double b = exp(-2. / a);
    double c = 1. / (1 + 2 * a * (1 - b) * u1_1);
    double u = u1_2
            + 0.25 * c / a
            * (c* (b * u1_1 * (1 + 2 * a * u1_1) + 8 * a * b * u1_2
                    + 32 * a * a * a * (1 - b) * (1 - b) * u1_2
                    * u1_2)
                    - 24 * a * a * (1 - b) * u1_3);
    return u;
}

double Coulomb3DLinkAction::calculateU3(double taueff, double reff) const {
    double u1_1 = coulomb1D.calculateU1(reff);
    double u1_2 = coulomb1D.calculateU2(reff);
    double u1_3 = coulomb1D.calculateU3(reff);
    double u1_4 = coulomb1D.calculateU4(reff);

    double a = 2. * taueff / (reff * reff + 1e-200);
    double b = exp(-2. / a);
    double c = 1. / (1 + 2 * a * (1 - b) * u1_1);
    double u = u1_3 + c * c * c / (24. * a * a)
		* (128 * a * a * a * a * a * b * b * b * (4 * u1_2 * u1_2 * u1_2
		- 9 * u1_1 * u1_2 * u1_3+ 6 * u1_1 * u1_1 * u1_4)
		+ 64 * a * a * a * (a * u1_2 * (
		-8 * a * u1_2 * u1_2 + 9 * (1 + 2 * a * u1_1) * u1_3)
		- 3 * (1 + 2 * a * u1_1)* (1 + 2 * a * u1_1)
		* u1_4) + 2 * a * b * b
		* (2 * a * u1_1 * u1_1 * u1_1
		+ 96 * a * a * u1_2 * (u1_2 - 8 * a * a * u1_2 * u1_2
		+ 3 * a * u1_3) + 12 * a * u1_1 * (u1_2 - 6 * a * u1_3
		+ 144 * a * a * a * u1_2 * u1_3
		- 32 * a * a * u1_4) + u1_1 * u1_1
		* (1 - 1152 * a * a * a * a * u1_4))
		+ b * (4 * a * a * u1_1 * u1_1 * u1_1
		+ 12 * a * (u1_2 - 16 * a * a * u1_2 * u1_2
		+ 128 * a * a * a * a * u1_2 * u1_2 * u1_2+ 6 * a * u1_3- 96 * a * a* a
		* u1_2* u1_3+ 16 * a * a* u1_4)+ 4 * u1_1 * u1_1
		* (a+ 576 * a* a* a* a* a* u1_4)
		+ u1_1 * (1 + 24 * a * a * (u1_2 + 6 * a * u1_3
		        - 144 * a * a * a * u1_2 * u1_3 + 64 * a * a * u1_4))));
    return u;
}

double Coulomb3DLinkAction::calculateU4(double reff) const {
    double u1_4 = coulomb1D.calculateU4(reff);
    double u = u1_4; //Not exact
    return u;
}
