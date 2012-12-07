#include "FermionWeight.h"
#include "base/Paths.h"
#include "util/Permutation.h"

FermionWeight::FermionWeight()
    :   sign(1.0) {
}

FermionWeight::~FermionWeight() {
}

int FermionWeight::getPartitionCount() const {
    return 2;
}

void FermionWeight::evaluate(Paths* paths) {
    sign = paths->getGlobalPermutation().sign();
}

double FermionWeight::getValue(int i) const {
    return (i==0) ? 1.0 : sign;
}
