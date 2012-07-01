#ifndef COULOMB1DLINKACTION_H_
#define COULOMB1DLINKACTION_H_

class Coulomb1DLinkAction {
public:
	Coulomb1DLinkAction();
	virtual ~Coulomb1DLinkAction();

    static double calculate1DValueAtOrigin(double stau);
};

#endif
