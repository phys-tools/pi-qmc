#ifdef __HAVE_CONFIG_H_
#include "config.h"
#endif

#include "GroundStateSNode.h"
#include "base/Species.h"
#include <cstdlib>
#include <blitz/tinyvec-et.h>

GroundStateSNode::GroundStateSNode(const Species& s, const int icenter, 
                                   const double tau) 
  : ifirst(s.ifirst), npart(s.count), icenter(icenter), 
    dunit(sqrt(2.0*s.mass/tau)) {
  std::cout << "GroundState S-Node with center = " << icenter << std::endl;
}

GroundStateSNode::~GroundStateSNode() {
}

NodeModel::DetWithFlag
GroundStateSNode::evaluate(const VArray &r1, const VArray &r2, 
                           const int islice, bool scaleMagnitude) {
  DetWithFlag result; result.err=false;
  double f=sqrt(dot(r1(ifirst+1)-r1(icenter),r1(ifirst+1)-r1(icenter))) 
          -sqrt(dot(r1(ifirst)-r1(icenter),r1(ifirst)-r1(icenter)));
  f*=sqrt(dot(r2(ifirst+1)-r2(icenter),r2(ifirst+1)-r2(icenter))) 
    -sqrt(dot(r2(ifirst)-r2(icenter),r2(ifirst)-r2(icenter)));
  result.det = f;
  return result;
}

void GroundStateSNode::evaluateDistance(const VArray &r1, const VArray &r2,
                             const int islice, Array& d1, Array& d2) {
  d1 = 100.;
  d2 = 100.;
  d1(ifirst) = d1(ifirst+1)
    = fabs(sqrt(dot(r1(ifirst)-r1(icenter),r1(ifirst)-r1(icenter)))
          -sqrt(dot(r1(ifirst+1)-r1(icenter),r1(ifirst+1)-r1(icenter))))
      *dunit;
  d2(ifirst) = d2(ifirst+1)
    = fabs(sqrt(dot(r2(ifirst)-r2(icenter),r2(ifirst)-r2(icenter)))
          -sqrt(dot(r2(ifirst+1)-r2(icenter),r2(ifirst+1)-r2(icenter))))
      *dunit;
}

void GroundStateSNode::evaluateGradLogDist(const VArray &r1, const VArray &r2,
          const int islice, VMatrix &gradd1, VMatrix& gradd2,
          const Array &dist1, const Array& dist2) {
/*  force=0;
  double x=fabs(sqrt(dot(r1(ifirst+1)-r1(icenter),r1(ifirst+1)-r1(icenter)))
               -sqrt(dot(r1(ifirst)-r1(icenter),r1(ifirst)-r1(icenter))));
  double y=fabs(sqrt(dot(r2(ifirst+1)-r2(icenter),r2(ifirst+1)-r2(icenter)))
               -sqrt(dot(r2(ifirst)-r2(icenter),r2(ifirst)-r2(icenter))));
  if (x<y) {
    force(ifirst)=r1(ifirst)-r1(icenter);
    force(ifirst+1)=r1(ifirst+1)-r1(icenter);
    double a=sqrt(dot(force(ifirst),force(ifirst)));
    double b=sqrt(dot(force(ifirst+1),force(ifirst+1)));
    force(ifirst)/=a*(a>b)?1:-1;
    force(ifirst+1)/=b*(b>a)?1:-1;
    force(icenter)=-(force(ifirst)+force(ifirst+1));
  } */
}
