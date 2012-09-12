#ifndef __AccRejEstimator_h_
#define __AccRejEstimator_h_
#include "Estimator.h"
#include "ReportWriters.h"
#include <cstdlib>
#include <blitz/array.h>
class Paths;
/// Base class for scalar estimators.
/// @version $Revision$
/// @author John Shumway
/// @bug Need to add results for parallel MPI runs.
class AccRejEstimator: public Estimator {
public:
    typedef blitz::Array<long int, 1> IArray;
    /// Construct by giving the name and number of levels.
    AccRejEstimator(const std::string& name, const int nlevel);
    /// Virtual destructor.
    virtual ~AccRejEstimator() {
    }
    /// Get the number of levels.
    int getNLevel() const {
        return nlevel;
    }
    /// Get the number of accept.
    long int getNAccept(const int i) const {
        return naccept(i);
    }
    /// Get the number of trials.
    long int getNTrial(const int i) const {
        return ntrial(i);
    }
    /// Get the array for number of acceptances.
    IArray& getNAccept() {
        return naccept;
    }
    /// Get the number of trials.
    IArray& getNTrial() {
        return ntrial;
    }
    /// Get the array for number of acceptances.
    const IArray& getNAccept() const {
        return naccept;
    }
    /// Get the number of trials.
    const IArray& getNTrial() const {
        return ntrial;
    }
    /// Average over clones.
    virtual void averageOverClones(const MPIManager* mpi);
    /// Callback EstimatorReportBuilder method for integer arrays.
    virtual void startReport(ReportWriters *writers) {
        writers->startAccRejReport(this, 0);
    }
    /// Callback EstimatorReportBuilder method for integer arrays.
    virtual void reportStep(ReportWriters* writers) {
        writers->reportAccRejStep(this, 0);
    }
    virtual void evaluate(const Paths&) {
        ;
    }
    /// Record a trial move.
    void tryingMove(const int ilevel = 0) {
        ++ntrial(ilevel);
    }
    /// Record an accepted move.
    void moveAccepted(const int ilevel = 0) {
        ++naccept(ilevel);
    }
protected:
    /// Number of levels.
    const int nlevel;
    /// Arrays for
    IArray naccept;
    IArray ntrial;
    IArray sum;
};
#endif
