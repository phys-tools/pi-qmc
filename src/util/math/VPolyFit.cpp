#include "VPolyFit.h"

VPolyFit::VPolyFit(int dataCount, int dimension, const double* xdata,
        const double* ydata)
    :   dataCount(dataCount),
        dimension(dimension) {
    solution = new double[dimension];
}

VPolyFit::~VPolyFit() {
    delete [] solution;
}

void VPolyFit::fit() {
}

const double* VPolyFit::getSolution() const {
    int ntot = dataCount * dimension;
//      int none = 2*shape[1]*shape[2];
//      int one = 1;
//      //y0 = y(n-1,all,all);
//      //c(blitz::Range(0,n-1),all,all) = y;
//      //d(blitz::Range(0,n-1),all,all) = y;
//      DCOPY_F77(&none,(double*)&y(n-1,0,0),&one,(double*)y0.data(),&one);
//      DCOPY_F77(&ntot,(double*)y.data(),&one,(double*)c.data(),&one);
//      DCOPY_F77(&ntot,(double*)y.data(),&one,(double*)d.data(),&one);
//      for (int j=1; j<n; ++j) {
//        for (int i=0; i<n-j; ++i) {
//          double denom=1./(x(i)-x(i+j));
//          DCOPY_F77(&none,(double*)&c(i+1,0,0),&one,(double*)a.data(),&one);
//          DSCAL_F77(&none,&denom,(double*)a.data(),&one);
//          denom *= -1;
//          DAXPY_F77(&none,&denom,(double*)&d(i,0,0),&one,(double*)a.data(),&one);
//          DCOPY_F77(&none,(double*)a.data(),&one,(double*)&c(i,0,0),&one);
//          DSCAL_F77(&none,&x(i),(double*)&c(i,0,0),&one);
//          DCOPY_F77(&none,(double*)a.data(),&one,(double*)&d(i,0,0),&one);
//          DSCAL_F77(&none,&x(i+j),(double*)&d(i,0,0),&one);
//    /*      for (int i1=0; i1<n1; ++ i1) {
//            for (int i2=0; i2<n2; ++ i2) {
//              a(i1,i2) = (c(i+1,i1,i2)-d(i,i1,i2))*denom;
//              c(i,i1,i2) = x(i)*a(i1,i2);
//              d(i,i1,i2) = x(i+j)*a(i1,i2);
//            }
//          } */
//        }
//        //y0 += d(n-j-1,all,all);
//        double unity=1.;
//        DAXPY_F77(&none,&unity,(double*)&d(n-j-1),&one,(double*)y0.data(),&one);
//      }
//      diff = d(0,all,all);
    solution[0] = -1.0;
    solution[1] = 0.5;
    return solution;
}


