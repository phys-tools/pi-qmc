#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ArrayEstimator.h"

ArrayEstimator::ArrayEstimator(const std::string &name,
   const std::string &typeString, bool hasError)
 : Estimator(name,typeString,""), hasErrorFlag(hasError) {
}

ArrayEstimator::~ArrayEstimator() {}
