#ifndef COULOMBLINKACTION_H_
#define COULOMBLINKACTION_H_

#include <config.h>
#include <blitz/array.h>

class CoulombLinkAction {
public:
    typedef blitz::TinyVector<double,NDIM> Vec;

	CoulombLinkAction(double q1q2, double epsilon, double deltaTau);
	virtual ~CoulombLinkAction();

	double getValue(Vec delta1, Vec delta2) const;

private:
	double stau;
};

#endif
