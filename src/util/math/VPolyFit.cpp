#include "VPolyFit.h"
#include <config.h>

#define DCOPY_F77 F77_FUNC(dcopy,DCOPY)
extern "C" void DCOPY_F77(const int *n, const double *x, const int *incx,
                                        const double *y, const int *incy);

#define DSCAL_F77 F77_FUNC(dscal,DSCAL)
extern "C" void DSCAL_F77(const int *n, const double *da, const double *dx,
  const int *incx);

#define DAXPY_F77 F77_FUNC(daxpy,DAXPY)
extern "C" void DAXPY_F77(const int *n, const double *da, const double *dx,
  const int *incx, const double *dy, const int *incy);

VPolyFit::VPolyFit(int dataCount, int dimension, const double* xdata,
        const double* ydata)
    :   xdata(xdata),
        ydata(ydata),
        dataCount(dataCount),
        dimension(dimension) {
    solution = new double[dimension];
    worka = new double[dataCount * dimension];
    workc = new double[dataCount * dimension];
    workd = new double[dataCount * dimension];
    y0 = new double[dimension];
}

VPolyFit::~VPolyFit() {
    delete[] solution;
    delete[] worka;
    delete[] workc;
    delete[] workd;
    delete[] y0;
}

void VPolyFit::fit() {
}

const double* VPolyFit::getSolution() const {
    int ntot = dataCount * dimension;
    int one = 1;
    const double* lasty = ydata + (dataCount - 1) * dimension;

    DCOPY_F77(&dimension, lasty, &one, y0, &one);
    DCOPY_F77(&ntot, ydata, &one, workc, &one);
    DCOPY_F77(&ntot, ydata, &one, workd, &one);

    for (int j = 1; j < dataCount; ++j) {
        for (int i = 0; i < dataCount - j; ++i) {
            double denom = 1. / (xdata[i] - xdata[i + j]);
            double* ciplus1 = workc + (i + 1) * dimension;
            DCOPY_F77(&dimension, ciplus1, &one, worka, &one);

            DSCAL_F77(&dimension, &denom, worka, &one);
            denom *= -1;
            double* di = workd + i * dimension;
            double* ci = workc + i * dimension;
            DAXPY_F77(&dimension, &denom, di, &one, worka, &one);
            DCOPY_F77(&dimension, worka, &one, ci, &one);
            DSCAL_F77(&dimension, xdata + i, ci, &one);
            DCOPY_F77(&dimension, worka, &one, di, &one);
            DSCAL_F77(&dimension, xdata + i + j, di, &one);
        }
        double unity = 1.0;
        double* lastd = workd + (dataCount - j - 1) * dimension;
        DAXPY_F77(&dimension, &unity, lastd, &one, y0, &one);
    }
    //    diff = workd(0,all);
    return y0;
}

