#ifndef __EwaldCoulombEstimator_h_
#define __EwaldCoulombEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/ScalarEstimator.h"
#include "util/EwaldSum.h"
#include <cstdlib>
#include <blitz/array.h>
class Paths;
class Action;
class SimulationInfo;
class EwaldSum;
class SuperCell;
class ScalarAccumulator;
/** Coulomb energy estimator with Ewald sum. 
 *  @author John Shumway  */
class EwaldCoulombEstimator: public ScalarEstimator, public LinkSummable {
public:
    typedef blitz::Array<double, 1> Array;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::Array<Vec, 1> VArray;
    /// Constructor.trad ewald
    EwaldCoulombEstimator(const SimulationInfo& simInfo, const Action*,
            const double epsilon, const double rcut, const double kcut,
            const std::string& unitName, double scale, double shift,
            const double kappa, const int nImages, const bool testEwald,
            ScalarAccumulator*);
    /// Constructor. optewald
    EwaldCoulombEstimator(const SimulationInfo& simInfo, const Action*,
            const double epsilon, const double rcut, const double kcut,
            const std::string& unitName, double scale, double shift,
            ScalarAccumulator*);
    virtual ~EwaldCoulombEstimator();
    virtual void initCalc(const int nslice, const int firstSlice);
    virtual void handleLink(const Vec& start, const Vec& end, const int ipart,
            const int islice, const Paths&);
    virtual void endCalc(const int nslice);
    virtual double calcValue();
    virtual void reset();
    virtual void evaluate(const Paths& paths);
    void findBoxImageVectors(const SuperCell &a);
    void testEwaldTotalCharge(const Paths& paths);
private:
    const bool testEwald;
    double kappa;
    double kcut;
    // const double kappa;
    /// Ewald sum object.
    EwaldSum &ewaldSum;
    /// The supercell.
    const SuperCell &cell;
    /// Radial array for short range potential.
    Array vgrid;
    /// Number of points on radial array.
    int nradial;
    /// Spacing for radial array.
    double rcut, dr, drinv;
    /// Action.
    const Action* action;
    ///dielectric constant
    const double epsilon;
    /// The charges of the particles.
    Array q;
    /// Buffer to hold current particle postions.
    VArray r;
    std::vector<std::vector<double> > boxImageVecs;
    const int nImages;
    double sphereR;
};

#endif
