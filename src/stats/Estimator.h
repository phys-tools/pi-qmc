#ifndef _Estimator_h__
#define _Estimator_h__

#include <string>
#include <iostream>
class MPIManager;
class Paths;
class ReportWriters;
/** Base class for all estimators.

 Estimators should be subclassed to handle different
 data types.

 @author John Shumway */
class Estimator {
public:
    Estimator(const std::string& name) :
            name(name) {
    }
    Estimator(const std::string &name, const std::string &typeString,
            const std::string &unitName) :
            name(name), typeString(typeString), unitName(unitName) {
    }
    virtual ~Estimator() {
    }

    const std::string& getName() const {
        return name;
    }

    const std::string& getUnitName() const {
        return unitName;
    }

    const std::string& getTypeString() const {
        return typeString;
    }

    virtual void evaluate(const Paths&)=0;

    virtual void averageOverClones(const MPIManager*) {
    }

    virtual void startReport(ReportWriters* builder) = 0;
    virtual void reportStep(ReportWriters* builder) = 0;

private:
    /// The name of the estimator.
    const std::string name;
    /// The type string.
    const std::string typeString;
    /// The unit name.
    const std::string unitName;
};
#endif
