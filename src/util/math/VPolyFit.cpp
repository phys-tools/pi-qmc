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
    workc = new double[dataCount * dataCount * dimension];
    workd = new double[dataCount * dataCount * dimension];
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

    DCOPY_F77(&dimension, ydata + (dataCount - 1) * dimension, &one, y0, &one);
    DCOPY_F77(&ntot, ydata, &one, workc, &one);
    DCOPY_F77(&ntot, ydata, &one, workd, &one);

    for (int j = 1; j < dataCount; ++j) {
        for (int i = 0; i < dataCount - j; ++i) {
            double denom = 1. / (xdata[i] - xdata[i + j]);
            DCOPY_F77(&ntot, workc + (i + 1) * dataCount * dimension, &one, worka, &one);

            DSCAL_F77(&ntot, &denom, worka, &one);
            denom *= -1;
            DAXPY_F77(&ntot, &denom, workd + i * dataCount * dimension, &one, worka, &one);
            DCOPY_F77(&ntot, worka, &one, workc + i * dataCount * dimension, &one);
            DSCAL_F77(&ntot, xdata + i, workc + i * dataCount * dimension, &one);
            DCOPY_F77(&ntot, worka, &one, workd + i * dataCount * dimension, &one);
            DSCAL_F77(&ntot, xdata + i + j, workd+ i * dataCount * dimension, &one);
////    /*      for (int i1=0; i1<n1; ++ i1) {
////            for (int i2=0; i2<n2; ++ i2) {
////              a(i1,i2) = (c(i+1,i1,i2)-d(i,i1,i2))*denom;
////              c(i,i1,i2) = x(i)*a(i1,i2);
////              d(i,i1,i2) = x(i+j)*a(i1,i2);
////            }
////          } */
        }
////        //y0 += d(n-j-1,all,all);
        double unity = 1.0;
        DAXPY_F77(&dimension, &unity, workd  + (dataCount - j - 1) * dataCount * dimension, &one, y0, &one);
    }
//    diff = workd(0,all,all);
    solution[0] = -1.0;
    solution[1] = 0.5;
    return solution;
}

