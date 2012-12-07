#include <gtest/gtest.h>
#include "util/fft/FFT1D.h"
#include <cmath>
#include <complex>

namespace {
typedef std::complex<double> Complex;

class FFT1DTest: public testing::Test {
protected:

    void SetUp() {
        dataCount = 16;
        twoPiOverDataCount = 2.0 * PI / dataCount;
        dataCountInverse = 1.0 / dataCount;
        data = new Complex[dataCount];
        for (int i = 0; i < dataCount; ++i) data[i] = 0.0;
        fft1D = new FFT1D(data, dataCount);
    }

    void TearDown() {
        delete fft1D;
        delete [] data;
    }

    Complex *data;
    int dataCount;
    double dataCountInverse;
    double twoPiOverDataCount;
    FFT1D *fft1D;
    static const double PI;
    static const Complex I;
};

const double FFT1DTest::PI = 3.141592653589793;
const Complex FFT1DTest::I = Complex(0.0, 1.0);

TEST_F(FFT1DTest, TestTransformDeltaFunction) {
    int index1 = 5;
    data[index1] = 3.0;
    int index2 = 3;
    Complex expect = data[5] * exp(-index1 * index2 * twoPiOverDataCount * I);
    fft1D->forward();
    ASSERT_NEAR(real(expect), real(data[index2]), 1e-14);
    ASSERT_NEAR(imag(expect), imag(data[index2]), 1e-14);
}

TEST_F(FFT1DTest, TestForwardAndReverseTransform) {
    int index = 5;
    double value = 3.0;
    data[index] = value;
    fft1D->forward();
    fft1D->reverse();
    ASSERT_NEAR(value * dataCount, real(data[index]), 1e-14);
}

}
