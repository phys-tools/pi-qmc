#ifndef COULOMB1DLINKACTION_H_
#define COULOMB1DLINKACTION_H_

class Coulomb1DLinkAction {
public:
	Coulomb1DLinkAction();
	virtual ~Coulomb1DLinkAction();

    static double calculateValueAtOrigin(double stau);
    static double calculateU0(
    		double stau, double u0, double reff, double taueff);
};

#endif
