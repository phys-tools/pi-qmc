#include "Coulomb1DLinkAction.h"

Coulomb1DLinkAction::Coulomb1DLinkAction() {

}

Coulomb1DLinkAction::~Coulomb1DLinkAction() {
}

double Coulomb1DLinkAction::calculate1DValueAtOrigin(double stau) {
    double u0 =
            stau * (1.772453851
                    + stau * (-0.074137740
                            + stau * (0.005834805
                                    + stau * (-0.000382686
                                            + stau * (0.000008738
                                                    + stau * 0.000002138)))));
    return u0;
}

