#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "Measure.h"
#include "stats/Estimator.h"
#include "stats/PartitionWeight.h"
#include "util/shiny/Shiny.h"
#include <iostream>

Measure::Measure(Paths& paths, std::vector<Estimator*> estimator,
        PartitionWeight* weight)
:   paths(paths),
    estimator(estimator),
    weight(weight) {
}

void Measure::run() {
    PROFILE_BEGIN(Measure);
    if (weight) {
        weight->evaluate(&paths);
    }
    for (unsigned int i = 0; i < estimator.size(); ++i) {
        estimator[i]->evaluate(paths);
    }
    PROFILE_END();
}
