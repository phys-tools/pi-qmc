#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "Collect.h"
#include "stats/EstimatorManager.h"
#include "util/shiny/Shiny.h"

Collect::Collect(std::string& estName, EstimatorManager& estManager,  
          const int nstep) 
  : estManager(estManager), nstep(nstep), isFirstRun(true) {
}

void Collect::run() {
  PROFILE_BEGIN (Collect);
  if (isFirstRun) {
    estManager.startWritingGroup(nstep,"estimators");
    isFirstRun=false;
  }
  estManager.writeStep();
  PROFILE_END();
}
