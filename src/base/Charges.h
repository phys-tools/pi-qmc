#ifndef CHARGES_H_
#define CHARGES_H_

class SimulationInfo;

class Charges {
public:
    Charges(const SimulationInfo*);
    virtual ~Charges();
    double getValue(int index) const;
private:
    double *value;
};

#endif
