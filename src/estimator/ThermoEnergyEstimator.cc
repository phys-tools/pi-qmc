#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "ThermoEnergyEstimator.h"
#include "action/Action.h"
#include "action/DoubleAction.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "stats/ScalarAccumulator.h"
#include "util/shiny/Shiny.h"
#include <cstdlib>
#include <blitz/tinyvec.h>

ThermoEnergyEstimator::ThermoEnergyEstimator(const SimulationInfo& simInfo,
        const Action* action, const DoubleAction* doubleAction,
        const std::string& unitName, double scale, double shift,
        ScalarAccumulator *accumulator)
:   ScalarEstimator("thermo_energy", "scalar-energy/thermo-energy",
                unitName, scale, shift),
    action(action),
    doubleAction(doubleAction) {
    this->accumulator = accumulator;
}

ThermoEnergyEstimator::~ThermoEnergyEstimator() {
}

void ThermoEnergyEstimator::initCalc(const int nslice, const int firstSlice) {
    PROFILE_BEGIN(ThermoEnergyEstimator);
    accumulator->clearValue();
}

void ThermoEnergyEstimator::handleLink(const Vec& start, const Vec& end,
        const int ipart, const int islice, const Paths& paths) {
    PROFILE_BEGIN(HandleLink);
    if (action) {
        PROFILE_BEGIN(Action);
        double u(0), utau(0), ulambda(0);
        Vec fm, fp;
        action->getBeadAction(paths, ipart, islice, u, utau, ulambda, fm, fp);
        accumulator->addToValue(utau);
        PROFILE_END();
    }

    if (doubleAction) {
        PROFILE_BEGIN(DoubleAction);
        double u(0), utau(0), ulambda(0);
        Vec fm, fp;
        doubleAction->getBeadAction(paths, ipart, islice, u, utau, ulambda, fm,
                fp);
        accumulator->addToValue(utau);
        PROFILE_END();
    }
    PROFILE_END();
}

void ThermoEnergyEstimator::endCalc(const int lnslice) {
    accumulator->storeValue(lnslice);
    PROFILE_END();
}

double ThermoEnergyEstimator::calcValue() {
    return 0.0;
}

void ThermoEnergyEstimator::reset() {
    accumulator->reset();
}

void ThermoEnergyEstimator::evaluate(const Paths& paths) {
    paths.sumOverLinks(*this);
}


