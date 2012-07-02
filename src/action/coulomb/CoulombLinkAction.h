#ifndef COULOMBLINKACTION_H_
#define COULOMBLINKACTION_H_

#include <config.h>
#include <blitz/array.h>

class CoulombLinkAction {
public:
    typedef blitz::TinyVector<double,NDIM> Vec;

	CoulombLinkAction(double q1q2, double epsilon, double mu,
			double deltaTau, int norder);
	virtual ~CoulombLinkAction();

	double getValue(Vec delta1, Vec delta2) const;

	static double calculateAverageSeparation(Vec delta1, Vec delta2);
	static double calculateS2(Vec delta1, Vec delta2);

private:
	const double q1q2;
	const double epsilon;
	const double mu;
	const double deltaTau;
	const double stau;
	const int norder;
};

#endif
