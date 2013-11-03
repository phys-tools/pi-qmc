#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "DoubleActionChoice.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/EnumeratedModelState.h"
#include "base/Paths.h"

DoubleActionChoice::DoubleActionChoice(const int n)
  : CompositeDoubleAction(n) {
  enumModelState = new EnumeratedModelState(n);
  modelState = enumModelState;
}

DoubleActionChoice::~DoubleActionChoice() {
  delete modelState;
}

double DoubleActionChoice::getActionDifference(
    const SectionSamplerInterface& sampler, const int level) {
  int imodel = enumModelState->getModelState();
  double diff = actions[imodel]->getActionDifference(sampler,level);
  return diff;
}

double DoubleActionChoice::getActionDifference(const Paths &paths, 
    const VArray &displacement, int nmoving, const IArray &movingIndex, 
    int iFirstSlice, int nslice) {
  int imodel = enumModelState->getModelState();
  double diff = actions[imodel]->getActionDifference(paths,displacement,nmoving,
                                           movingIndex,iFirstSlice,nslice);
  return diff;
}

double DoubleActionChoice::getTotalAction(const Paths& paths, int level) const {
  int imodel = enumModelState->getModelState();
  double total = actions[imodel]->getTotalAction(paths,level);
  return total;
}

void DoubleActionChoice::getBeadAction(const Paths& paths, 
       int ipart, int islice, double& u, double& utau, double& ulambda,
       Action::Vec& fm, Action::Vec& fp, bool check_for_node_crossing) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  int imodel = enumModelState->getModelState();
  actions[imodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp,
          check_for_node_crossing);
}

void DoubleActionChoice::initialize(const DoubleSectionChooser& 
                                             sectionChooser) {
  int imodel = enumModelState->getModelState();
  actions[imodel]->initialize(sectionChooser);
}

void DoubleActionChoice::acceptLastMove() {
  int imodel = enumModelState->getModelState();
  actions[imodel]->acceptLastMove();
}

double DoubleActionChoice::getActionChoiceDifference(const Paths &paths, int j) {
 jmodel = j;
 paths.sumOverLinks(*this);
 return actionDifference;
}

void DoubleActionChoice::initCalc(const int nslice, const int firstSlice) {
  actionDifference = 0.;
}

void DoubleActionChoice::handleLink(const LinkSummable::Vec &start,
    const LinkSummable::Vec &end, const int ipart, const int islice,
    const Paths &paths) {
  double u=0., utau=0, ulambda=0;
  LinkSummable::Vec fm=0.; LinkSummable::Vec fp=0.;
  int imodel = enumModelState->getModelState();
  actions[imodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp, false);
  actionDifference -= u;
  actions[jmodel]->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp, false);
  actionDifference += u;
}

void DoubleActionChoice::endCalc(const int nslice) {
}

