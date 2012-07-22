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
class MPIManager;
class EwaldSum;
class SuperCell;
/** Coulomb energy estimator with Ewald sum. 
 *  @version $Revision$
 *  @author John Shumway  */
class EwaldCoulombEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.trad ewald
    EwaldCoulombEstimator(const SimulationInfo& simInfo, const Action*, 
			  const double epsilon, const double rcut, const double kcut,
			  MPIManager *mpi,
			  const std::string& unitName, double scale, double shift, 
			  const double kappa, const int nImages, const bool testEwald);
 /// Constructor. optewald
    EwaldCoulombEstimator(const SimulationInfo& simInfo, const Action*, 
			  const double epsilon, const double rcut, const double kcut,
			  MPIManager *mpi,
			  const std::string& unitName, double scale, double shift);
    /// Virtual destructor.
  virtual ~EwaldCoulombEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of coulomb energy estimator.
  virtual double calcValue() {return value=etot/enorm;}
  /// Clear value of coulomb energy estimator.
  virtual void reset() {etot=enorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
  void findBoxImageVectors(  const SuperCell &a);
  void testEwaldTotalCharge( const Paths& paths);
private:
  const bool testEwald;
  double kappa;
  double kcut;
  // const double kappa;
  /// Ewald sum object.
  EwaldSum &ewaldSum;
  /// The supercell.
  const SuperCell &cell; 
  /// CoulombEnergy energy.
  double energy;
  /// The accumulated energy.
  double etot;
  /// The normalization of the accumualted energy.
  double enorm;
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
  /// The MPI manager.
  const MPIManager* mpi;
};

#endif
