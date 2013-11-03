#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "BondLengthEstimator.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/Units.h"
#include "util/SuperCell.h"
#include <blitz/tinyvec-et.h>

BondLengthEstimator::BondLengthEstimator(const SimulationInfo& simInfo,
  const Species& species1, const Species& species2, const std::string& unitName)
  : ScalarEstimator("bond_length_"+species1.name+"-"+species2.name,
      unitName, "scalar-length/bondlength",
      1./simInfo.getUnits()->getLengthScaleIn(unitName),0.),
    length(0), norm1(0), norm2(0), avlength(0), names(simInfo.getNPart()),
    spec1(species1.name), spec2(species2.name), cell(*simInfo.getSuperCell()) {
  for (int i=0; i<names.size(); ++i)
    names(i)=simInfo.getPartSpecies(i).name;
  std::cout << "Bond Length Estimator "
            << spec1 << " " << spec2 << std::endl;
}

void BondLengthEstimator::initCalc(const int nslice, const int firstSlice) {
  length=norm1=0;
}

void BondLengthEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  if(names(ipart)==spec1) {
    for (int jpart=0; jpart<names.size(); ++jpart) {
      if(names(jpart)==spec2 && jpart!=ipart) {
        Vec delta=end-paths(jpart,islice);
        cell.pbc(delta);
        length+=sqrt(dot(delta,delta));
        norm1++;
      }
    }
  }
}

void BondLengthEstimator::endCalc(const int nslice) {
   avlength+=(length/norm1);
   norm2++;
}
