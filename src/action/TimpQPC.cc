#include "config.h"
#include "TimpQPC.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Beads.h"
#include "base/Paths.h"
#include "base/Species.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include <fstream>

TimpQPC::TimpQPC(const SuperCell& cell, const Species &species, const double tau, 
              const double width, const double length, const double vG,
              const double z, MPIManager *mpi)
  : tau(tau), width(width), length(length), vG(vG), z(z),
    ifirst(species.ifirst), npart(species.count),
    lx(cell[0]), ly(cell[1]) {
  if (mpi && mpi->isMain())
  {
    std::cout << "Timp QPC on "<<species.name
              <<" with w=" << width << ", l=" << length
              << ", vG=" << vG << ", z=" << z << std::endl;
    std::ofstream file("timp.dat");
    for (int j=-50; j<50; ++j) {
      for (int i=-50; i<50; ++i) {
        file << v(i*lx*0.01,j*ly*0.01) << std::endl;
      }
      file << std::endl;
    }
    file.close();
  }
}

double TimpQPC::getActionDifference(const SectionSamplerInterface& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const int nStride = 1 << level;
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  for (int islice=nStride; islice<nSlice; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      if (i>=ifirst && i<ifirst+npart) {
        // Add action for moving beads.
        Vec r=movingBeads(iMoving,islice); 
        //cell.pbc(r);
        deltaAction += tau*nStride*v(r[0],r[1]);
        // Subtract action for old beads.
        r=sectionBeads(i,islice);
        //cell.pbc(r);
        deltaAction -= tau*nStride*v(r[0],r[1]);
      }
    }
  }
  return deltaAction;
}

double TimpQPC::getTotalAction(const Paths& paths, const int level) const {
  return 0;
}

void TimpQPC::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  u=utau=0; fm=0.; fp=0.;
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec r=paths(ipart,islice);
    utau = v(r[0],r[1]);
    u = utau*tau;
  }
}

double TimpQPC::getActionDifference(const Paths &paths, 
       const VArray &displacement, int nmoving, const IArray &movingIndex,
       int iFirstSlice, int iLastSlice) {
  double deltaAction = 0;
  const SuperCell& cell = paths.getSuperCell();
  for (int i=0; i<nmoving; ++i) {
    int ipart = movingIndex(i);
    if (ipart>=ifirst && ipart<ifirst+npart) {
      for (int islice=iFirstSlice; islice<=iLastSlice; ++islice) {
        Vec r = paths(ipart,islice);
        deltaAction -= tau*v(r[0],r[1]);
        r += displacement(i);
	cell.pbc(r);
        deltaAction += tau*v(r[0],r[1]);
      }
    }
  }
  return deltaAction;
}


double TimpQPC::v(double x, double y) const {
  return f((2*x-length)*0.5/z,( 2*y+width)*0.5/z)
        -f((2*x+length)*0.5/z,( 2*y+width)*0.5/z)
        +f((2*x-length)*0.5/z,(-2*y+width)*0.5/z)
        -f((2*x+length)*0.5/z,(-2*y+width)*0.5/z);
}

double TimpQPC::f(const double u, const double v) const {
  return -vG/(2*PI)*(0.5*PI-atan(u)-atan(v)+atan(u*v/(sqrt(1+u*u+v*v))));
}

const double TimpQPC::PI = acos(-1.0);
