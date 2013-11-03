#ifndef MAGNETICFLUXCALCULATOR_H_
#define MAGNETICFLUXCALCULATOR_H_

#include "LinkSummable.h"
class Charges;
class Paths;

class MagneticFluxCalculator: public LinkSummable {
public:
    MagneticFluxCalculator(Charges*);
    virtual ~MagneticFluxCalculator();
    double calculate(Paths* paths);
    virtual void initCalc(const int nslice, const int firstSlice);
    virtual void handleLink(const Vec& start, const Vec& end, int ipart,
            int islice, const Paths&);
    virtual void endCalc(int nslice);
private:
    Charges* charges;
    double flux;
};

#endif
