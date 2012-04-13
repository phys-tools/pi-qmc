#include "Hungarian.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define ASSNDX_F77 F77_FUNC(assndx,ASSNDX)
extern "C" void ASSNDX_F77(const int *mode, double *a, const int *n,
  const int *m, const int *ida, int *k, double *sum, int *iw, const int *idw);

Hungarian::Hungarian(int size)
    :   size(size), kindex(new int[size]) {
}

Hungarian::~Hungarian() {
    delete [] kindex;
}

int Hungarian::operator[](int index) const {
    return kindex[index];
}

double Hungarian::getSum() const {
    return sum;
}

void Hungarian::solve(double* matrix) {
    int work[6*size];
    ASSNDX_F77(&MODE,matrix,&size,&size,&size,kindex,&sum,work,&size);
    for (int index = 0; index < size; ++index) {
        kindex[index] -= 1;
    }
}

const int Hungarian::MODE = 1;
