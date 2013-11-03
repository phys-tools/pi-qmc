#ifndef __LinkSummable_h_
#define __LinkSummable_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <blitz/tinyvec.h>
class Paths;
/// Interface class for objects that can be summed over links.
/// @author John Shumway
class LinkSummable {
public:
    typedef blitz::TinyVector<double, NDIM> Vec;
    virtual ~LinkSummable() {
    }
    virtual void initCalc(const int nslice, const int firstSlice) {
    }
    virtual void handleLink(const Vec& start, const Vec& end, int ipart,
            int islice, const Paths&) = 0;
    virtual void endCalc(int nslice) {
    }
};
#endif
