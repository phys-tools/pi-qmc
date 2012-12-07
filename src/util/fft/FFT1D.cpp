#include "FFT1D.h"

FFT1D::FFT1D(Complex* data, int dataCount)
    :   data(data),
        dataCount(dataCount) {

    fftw_complex *pointer = (fftw_complex*)data;

    forwardPlan = fftw_plan_dft_1d(
            dataCount, pointer, pointer, FFTW_FORWARD, FFTW_ESTIMATE);

    reversePlan = fftw_plan_dft_1d(
            dataCount, pointer, pointer, FFTW_BACKWARD, FFTW_ESTIMATE);
}

FFT1D::~FFT1D() {
    fftw_destroy_plan(reversePlan);
    fftw_destroy_plan(forwardPlan);
}

void FFT1D::forward() {
    fftw_execute(forwardPlan);
}

void FFT1D::reverse() {
    fftw_execute(reversePlan);
}

