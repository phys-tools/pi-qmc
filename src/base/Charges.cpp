#include "Charges.h"

#include "SimulationInfo.h"

Charges::Charges(const SimulationInfo* simInfo) {
    int particleCount = simInfo->getNPart();
    value = new double[particleCount];
    for (int index = 0; index < particleCount; ++index) {
        value[index] = simInfo->getPartSpecies(index).charge;
    }
}

Charges::~Charges() {
    delete[] value;
}

double Charges::getValue(int index) const {
    return value[index];
}

