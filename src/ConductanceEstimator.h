// $Id$
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifndef __ConductanceEstimator_h_
#define __ConductanceEstimator_h_
#include "stats/BlitzArrayBlkdEst.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <fftw3.h>
class Paths;
class Species;
class SimulationInfo;
class MPIManager;
/** Current-current correlation estimator for a homogeneous system.
 *  Calculates the current-current correlation function,
 *  @f[
 *  \chi_{J_\alpha^A  J_\beta^B} = -\frac{1}{\hbar}
 *  \langle J_\alpha^A J_\beta^B \rangle
 *  @f]
 *  where the @f$ J_\alpha^A @f$ are the particle (or charge) current
 *  in direction @f$\alpha@f$ for species @f$A@f$.
 *
 *  Results are stored in a rank 6 array,
 *  @f[ 
 *  \chi_{JJ}(i_{\text{order}},\alpha,\beta,A,B,i\omega_n).
 *  @f] 
 *  Flags control the behavior:
 *  - bool useCharge uses charge currents
 *  - int idim  selects a dimension, idim=-1 computes tensor.
 *  - species selects a species
 *  - bool speciesTensor computes a tensor of cross-species correlations
 *  - int nfreq number of frequencies to store.
 *  @version $Revision$
 *  @author John Shumway  */
class ConductanceEstimator : public BlitzArrayBlkdEst<6>, public LinkSummable {
public:
  typedef blitz::Array<std::complex<double>,6> CArray6;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  /// Constructor.
  //ConductanceEstimator(const SimulationInfo& simInfo, const int nfreq,
  //                     MPIManager *mpi);
  ConductanceEstimator(const SimulationInfo& simInfo, const int nfreq,
                       Species* species, bool useSpeciesTensor,
                       int idim, bool useCharge, MPIManager *mpi,
                       int norder);
  /// Virtual destructor.
  virtual ~ConductanceEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const blitz::TinyVector<double,NDIM>& start,
                          const blitz::TinyVector<double,NDIM>& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  // Get value of kinetic energy estimator.
  //virtual double calcValue() {return value=etot/enorm;}
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  ///
  const int npart, ifirst, nslice, nfreq;
  const double tau;
  const int idim;
  CArray6 temp, buff;
  int n; 
  //IArray ninbin,ninbinbuff;
  fftw_plan fwd, rev;
  MPIManager *mpi;
  /// The charge of each particle.
  Array charge;
  /// The species bin number of the particle.
  IArray ispecies;
  /// A reference to the supercell.
  const SuperCell &cell;
  /// The order of the polarizability (norder>1 for hyperpolarizability)
  const int norder;
};

#endif
