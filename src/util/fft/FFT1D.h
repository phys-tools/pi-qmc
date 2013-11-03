#ifndef FFT1D_H_
#define FFT1D_H_

#include <complex>
#include <fftw3.h>

class FFT1D {
public:
    typedef std::complex<double> Complex;
    FFT1D(Complex *data, int dataCount);
    virtual ~FFT1D();

    void forward();
    void reverse();
private:
    Complex *data;
    const int dataCount;
    fftw_plan forwardPlan;
    fftw_plan reversePlan;
};

#endif
