#ifndef BLAS_H_
#define BLAS_H_

#include <config.h>

#define DCOPY_F77 F77_FUNC(dcopy,DCOPY)
extern "C" void DCOPY_F77(const int *n, const double *x, const int *incx,
        const double *y, const int *incy);

#define DSCAL_F77 F77_FUNC(dscal,DSCAL)
extern "C" void DSCAL_F77(const int *n, const double *a, double *x,
        const int *incx);

#define DAXPY_F77 F77_FUNC(daxpy,DAXPY)
extern "C" void DAXPY_F77(const int *n, const double *a, const double *x,
        const int *incx, double *y, const int *incy);


class BLAS {
public:
    static void dcopy(int n, const double *fromX, int incx, double *toY, int incy) {
        DCOPY_F77(&n, fromX, &incx, toY, &incy);
    }

    static void dscal(int n, double a, double *x, int incx) {
        DSCAL_F77(&n, &a, x, &incx);
    }

    static void daxpy(int n, double a, const double *x, int incx, double *y,
            int incy) {
        DAXPY_F77(&n, &a, x, &incx, y, &incy);
    }
};
#endif
