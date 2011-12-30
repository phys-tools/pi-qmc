#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"

#include "MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"


namespace {

class EMARateActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake();
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerInterface *sampler;
    Species species1, species2;
    SimulationInfo simInfo;
};

TEST_F(EMARateActionTest, testCreate) {
    EMARateAction action(simInfo, species1, species2, 1.0);
//    action.testableGetActionDifference(*sampler, 0);
}

}

