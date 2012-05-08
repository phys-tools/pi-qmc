#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "EMARateEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Paths.h"
#include "util/SuperCell.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

EMARateEstimator::EMARateEstimator(const SimulationInfo& simInfo, double C)
: ScalarEstimator("ema_rate","scalar/ema-rate","",1.,0.),
  dtau(simInfo.getTau()),
  masse(1.0),
  massh(1.0),
  C(C),
  cell(*simInfo.getSuperCell()),
  lastSlice(simInfo.getNSlice()-1),
  sum(0),
  norm(0) {
}

void EMARateEstimator::initCalc(const int nslice, const int firstSlice) {
    actionDifference = 0.;
}

void EMARateEstimator::handleLink(const Vec& start, const Vec& end,
        const int ipart, const int islice, const Paths& paths) {
    if (islice==0 && ipart == 1) {
        evaluateElectronBeforeRecombination(start, end, paths);
    }
    if (islice==1 && ipart == 0) {
        evaluateHoleAfterRecombination(start, end, paths);
    }
}

void EMARateEstimator::evaluateElectronBeforeRecombination(
        const Vec &start, const Vec &end, const Paths &paths) {
    Vec delta = paths(0, 0) - start;
    cell.pbc(delta);
    actionDifference += 0.5 * masse * dot(delta, delta) / dtau;

    delta = end - start;
    cell.pbc(delta);
    actionDifference -= 0.5 * masse * dot(delta, delta) / dtau;
}

void EMARateEstimator::evaluateHoleAfterRecombination(
        const Vec &start, const Vec &end, const Paths &paths) {
    Vec delta = end - paths(1, 1, -1);
    cell.pbc(delta);
    actionDifference += 0.5 * massh * dot(delta, delta) / dtau;

    delta = end - start;
    cell.pbc(delta);
    actionDifference -= 0.5 * massh * dot(delta, delta) / dtau;
}

void EMARateEstimator::endCalc(const int nslice) {
    if ( C > 0.0 && (actionDifference - log(C)) > -40.0) {
        sum += 1.0 / (1.0 + C * exp(-actionDifference));
    } else {
        sum += 0.0;
    }
    norm += 1.0;
}
