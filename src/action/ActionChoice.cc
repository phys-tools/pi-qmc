#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ActionChoice.h"
#include "advancer/SectionChooser.h"
#include "base/EnumeratedModelState.h"
#include "base/Paths.h"

ActionChoice::ActionChoice(const int n)
  : CompositeAction(n) {
  enumModelState = new EnumeratedModelState(n);
  modelState = enumModelState;
}

ActionChoice::~ActionChoice() {
  delete modelState;
}

double ActionChoice::getActionDifference(
    const SectionSamplerInterface& sampler, int level) {
  int imodel = enumModelState->getModelState();
  double diff = actions[imodel]->getActionDifference(sampler,level);
  return diff;
}

double ActionChoice::getActionDifference(const Paths &paths, 
    const VArray &displacement, int nmoving, const IArray &movingIndex, 
    int iFirstSlice, int iLastSlice) {
  int imodel = enumModelState->getModelState();
  double diff = actions[imodel]->getActionDifference(paths,displacement,
      nmoving, movingIndex,iFirstSlice,iLastSlice);
  return diff;
}

double ActionChoice::getTotalAction(const Paths& paths, int level) const {
  int imodel = enumModelState->getModelState();
  double total=actions[imodel]->getTotalAction(paths,level);
  return total;
}

void ActionChoice::getBeadAction(const Paths& paths, int ipart, int islice,
       double& u, double& utau, double& ulambda, Action::Vec& fm, 
       Action::Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  int imodel = enumModelState->getModelState();
  actions[imodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
}

void ActionChoice::initialize(const SectionChooser& sectionChooser) {
  int imodel = enumModelState->getModelState();
  actions[imodel]->initialize(sectionChooser);
}

void ActionChoice::acceptLastMove() {
  int imodel = enumModelState->getModelState();
  actions[imodel]->acceptLastMove();
}

double ActionChoice::getActionChoiceDifference(const Paths &paths, int j) {
  jmodel = j;
  paths.sumOverLinks(*this);
  return actionDifference;
}

void ActionChoice::initCalc(const int nslice, const int firstSlice) {
  actionDifference = 0.;
}

void ActionChoice::handleLink(const LinkSummable::Vec &start,
    const LinkSummable::Vec &end, const int ipart, const int islice,
    const Paths &paths) {
  double u=0., utau=0, ulambda=0;
  int imodel = enumModelState->getModelState();
  LinkSummable::Vec fm=0.; LinkSummable::Vec fp=0.;
  actions[imodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
  actionDifference -= u;
  actions[jmodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
  actionDifference += u;
}

void ActionChoice::endCalc(const int nslice) {
}
