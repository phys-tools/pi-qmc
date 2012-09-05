#ifndef __SHOInteraction_h_
#define __SHOInteraction_h_
#include "action/Action.h"
class SimulationInfo;
class SectionSamplerInterface;
class Paths;
class Species;
#include <cstdlib>
#include <blitz/array.h>

/** Spring connecting two particles.
 * @bug Not fully implemented! Only used for debugging EMARateAction.
 * @author John Shumway
 * */

class SHOInteraction: public Action {
public:
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::Array<double, 1> Array;
    typedef blitz::Array<double, 2> Array2;
    typedef blitz::Array<double, 3> Array3;

    SHOInteraction(const SimulationInfo &simInfo, double omega,
            const Species *species1, const Species *species2, int maxlevel);
    virtual ~SHOInteraction();

    virtual double getActionDifference(const SectionSamplerInterface&,
            int level);

    virtual double getActionDifference(const Paths&, const VArray &displacement,
            int nmoving, const IArray &movingIndex, int iFirstSlice,
            int iLastSlice);

    virtual double getTotalAction(const Paths&, const int level) const;

    virtual void getBeadAction(const Paths&, const int ipart, const int islice,
            double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;

    double calculateAction(Vec delta1, Vec delta2, int level) const;
    double calculateTauDerivativeOfAction(Vec delta1, Vec delta2) const;
private:
    const double omega;
    const double mu;
    const double deltaTau;
    Array coshwt;
    Array sinhwt;
};
#endif
