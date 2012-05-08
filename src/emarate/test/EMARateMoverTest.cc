#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"
#include "SimulationInfo.h"
#include "SerialPaths.h"
#include "BeadFactory.h"

namespace {

class EMARateMoverTest: public testing::Test {
protected:

    virtual void SetUp() {
        separation = 5.0;
        simInfo.nslice = 128;
        simInfo.npart = 2;
        paths = new SerialPaths(simInfo.npart, simInfo.nslice, simInfo.tau,
                *simInfo.superCell, beadFactory);
        coefficient = 10.0;
    }

    virtual void TearDown() {
        delete paths;
    }

    void setDirectPaths() {
        Paths::Vec holePosition = Paths::Vec(0.0, 0.0, 0.0);
        Paths::Vec electronPosition = Paths::Vec(separation, separation, separation);
        for (int sliceIndex = 0; sliceIndex < simInfo.nslice; ++sliceIndex) {
            (*paths)(0, sliceIndex) = holePosition;
            (*paths)(1, sliceIndex) = electronPosition;
        }
    }

    void setExchangingPaths() {
        Paths::Vec afterPosition = Paths::Vec(0.0, 0.0, 0.0);
        Paths::Vec beforePosition = Paths::Vec(separation, separation, separation);
        (*paths)(0,0) = beforePosition;
        (*paths)(1,0) = afterPosition;

        double inverseSliceCount = 1.0 / paths->getNSlice();
        for (int sliceIndex = 1; sliceIndex < simInfo.nslice; ++sliceIndex) {
            double x = sliceIndex * inverseSliceCount;
            Paths::Vec position = (1.0 - x) * afterPosition + x * beforePosition;
            (*paths)(0, sliceIndex) = position;
            (*paths)(1, sliceIndex) = position;
        }
    }

    Species species1, species2;
    double coefficient;
    SimulationInfo simInfo;
    SerialPaths *paths;
    BeadFactory beadFactory;
    double separation;
};

TEST_F(EMARateMoverTest, testDirectPaths) {
    setDirectPaths();
}

}
