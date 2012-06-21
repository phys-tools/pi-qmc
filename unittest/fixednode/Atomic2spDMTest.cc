#include <gtest/gtest.h>
#include <cmath>
#include "fixednode/Atomic2spDM.h"

namespace {

class Atomic2spDMTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        Z = 1.0;
        weight = 1.0;
        pweight = 0.3;
        npart = 1;
        ifirst = 0;
        nfermion = 1;
        orbitalDM = 0;
        scale = 1.0;
        ipart = 0;
        ifermion = 0;
    }

    virtual void TearDown() {
        delete orbitalDM;
        orbitalDM = 0;
    }

    AtomicOrbitalDM* orbitalDM;
    AtomicOrbitalDM::VArray2* d1;
    AtomicOrbitalDM::VArray2* d2;
    AtomicOrbitalDM::Matrix result;

    double Z;
    double weight;
    double pweight;
    int npart;
    int ifirst;
    int nfermion;
    double scale;
    int ipart;
    int ifermion;
    static const double PI;

    void createOrbitalDM() {
        orbitalDM =
                new Atomic2spDM(Z, ifirst, npart, nfermion, pweight, weight);
        d1 = &orbitalDM->getD1Array();
        d2 = &orbitalDM->getD2Array();
        (*d1)(ifermion, ipart) = 1.0;
        (*d2)(ifermion, ipart) = 3.0;
        result.resize(nfermion, nfermion);
        result = 0.0;
    }


#if NDIM==3
    void calculateGradient(AtomicOrbitalDM::Vec& grad1,
            AtomicOrbitalDM::Vec& grad2) {
        grad1[0] = calculateNumericalDerivative((*d1)(ifermion,ipart)[0]);
        grad1[1] = calculateNumericalDerivative((*d1)(ifermion,ipart)[1]);
        grad1[2] = calculateNumericalDerivative((*d1)(ifermion,ipart)[2]);
        grad2[0] = calculateNumericalDerivative((*d2)(ifermion,ipart)[0]);
        grad2[1] = calculateNumericalDerivative((*d2)(ifermion,ipart)[1]);
        grad2[2] = calculateNumericalDerivative((*d2)(ifermion,ipart)[2]);
    }
#endif

    double getOrbitalValue() {
        result = 0.0;
        orbitalDM->evaluateValue(result, scale);
        double value = result(0,0);
        result = 0.0;
        return value;
    }

    double calculateNumericalDerivative(double& x) {
        const double DELTA = 1e-8;
        x += DELTA;
        double deltaY = getOrbitalValue();
        x -= 2*DELTA;
        deltaY -= getOrbitalValue();
        x += DELTA;
        return deltaY / (2*DELTA);
    }
};

const double Atomic2spDMTest::PI = 3.141592653589793;


TEST_F(Atomic2spDMTest, testValueAtOriginIsZero) {
    createOrbitalDM();
    (*d1)(ifermion, ipart) = 0.0;
    (*d2)(ifermion, ipart) = 0.0;
    orbitalDM->evaluateValue(result, scale);
    ASSERT_FLOAT_EQ(0.0, result(0,0));
}

TEST_F(Atomic2spDMTest, testThatCallingEvaluateTwiceDoublesValue) {
    createOrbitalDM();
    orbitalDM->evaluateValue(result, scale);
    double expect = 2*result(0,0);
    orbitalDM->evaluateValue(result, scale);
    ASSERT_FLOAT_EQ(expect, result(0,0));
}

TEST_F(Atomic2spDMTest, testScalingOfValue) {
    createOrbitalDM();
    orbitalDM->evaluateValue(result, scale);
    scale = 2.3;
    double expect = scale*result(0,0);
    result =0.0;
    orbitalDM->evaluateValue(result, scale);
    ASSERT_FLOAT_EQ(expect, result(0,0));
}

TEST_F(Atomic2spDMTest, testValueAwayFromOrigin) {
    createOrbitalDM();
    (*d1)(ifermion, ipart) = 0.0;
    (*d2)(ifermion, ipart) = 0.0;
    (*d1)(ifermion, ipart)[0] = 1.0;
    (*d2)(ifermion, ipart)[0] = 3.0;
    orbitalDM->evaluateValue(result, scale);
    double expect = 1.0*3.0*exp(-0.5*(1.0+3.0))/(32*PI)*(1.0+0.3);
    ASSERT_FLOAT_EQ(expect, result(0,0));
}

#if NDIM==3
TEST_F(Atomic2spDMTest, testAngularValueAtNinetyDegrees) {
    createOrbitalDM();
    (*d1)(ifermion, ipart) = 0.0;
    (*d2)(ifermion, ipart) = 0.0;
    (*d1)(ifermion, ipart)[0] = 1.0;
    (*d2)(ifermion, ipart)[2] = 3.0;
    orbitalDM->evaluateValue(result, scale);
    double expect = 1.0*3.0*exp(-0.5*(1.0+3.0))/(32*PI);
    ASSERT_FLOAT_EQ(expect, result(0,0));
}

TEST_F(Atomic2spDMTest, testGradient) {
    pweight = 1.0;
    createOrbitalDM();
    (*d1)(ifermion, ipart) = AtomicOrbitalDM::Vec(1.0, 2.0, 3.0);
    (*d2)(ifermion, ipart) = AtomicOrbitalDM::Vec(2.0, 1.0, 0.0);
    double r1 = sqrt(14.0);
    double r2 = sqrt(5.0);
    double r1DotR2 = 4.0;
    AtomicOrbitalDM::Vec expectedGrad1, expectedGrad2;
    calculateGradient(expectedGrad1, expectedGrad2);
    orbitalDM->evaluateValueAndGrad();
    AtomicOrbitalDM::ValAndGrad result =
            orbitalDM->operator ()(ifermion, ifermion);
    double expect = r1*r2*exp(-0.5*(r1+r2))/(32*PI)*(1.0+pweight*r1DotR2/(r1*r2));
    ASSERT_FLOAT_EQ(expect, result.val);
    ASSERT_FLOAT_EQ(expectedGrad1[0], result.grad1[0]);
    ASSERT_FLOAT_EQ(expectedGrad1[1], result.grad1[1]);
    ASSERT_FLOAT_EQ(expectedGrad1[2], result.grad1[2]);
    ASSERT_FLOAT_EQ(expectedGrad2[0], result.grad2[0]);
    ASSERT_FLOAT_EQ(expectedGrad2[1], result.grad2[1]);
    ASSERT_FLOAT_EQ(expectedGrad2[2], result.grad2[2]);
}

#endif

}
