#ifndef __BlitzArrayBlkdEst_h_
#define __BlitzArrayBlkdEst_h_
#include "ArrayEstimator.h"
#include "EstimatorReportBuilder.h"
#include "MPIManager.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
class Paths;

template<int N>
class BlitzArrayBlkdEst: public ArrayEstimator {
public:
    typedef blitz::Array<float, N> ArrayN;
    typedef blitz::TinyVector<int, N> IVecN;
    typedef blitz::TinyVector<double, N> VecN;

    BlitzArrayBlkdEst(const std::string &name, const std::string &typeString,
            const IVecN &n, bool hasError);
    virtual ~BlitzArrayBlkdEst() {
    }

    virtual void reset()=0;
    virtual void averageOverClones(const MPIManager* mpi);
    virtual int getNDim() const {
        return rank;
    }

    virtual int getExtent(const int idim) const {
        return n(idim);
    }
    virtual const void* getData() const {
        return accumvalue.data();
    }

    virtual const void* getError() const {
        return hasErrorFlag ? accumvalue2.data() : (float*) 0;
    }

    virtual void normalize() const {
        accumvalue /= accumnorm; //Temporarily normalize the data to write it out.
        if (hasErrorFlag) { //Temperorily convert accumvalue2 into rms error.
            accumvalue2 /= accumnorm;
            accumvalue2 -= accumvalue * accumvalue;
            accumvalue2 = sqrt(fabs(accumvalue2 / iblock));
        }
    }

    virtual void unnormalize() const {
        if (hasErrorFlag) {
            accumvalue2 *= accumvalue2 * iblock;
            accumvalue2 += accumvalue * accumvalue;
            accumvalue2 *= accumnorm;
        }
        accumvalue *= accumnorm;
    }

    virtual bool hasScale() const {
        return scale != 0;
    }

    virtual bool hasOrigin() const {
        return origin != 0;
    }

    virtual const void* getScale() const {
        return scale;
    }

    virtual const void* getOrigin() const {
        return origin;
    }

    virtual bool hasMin() const {
        return min != 0;
    }

    virtual bool hasMax() const {
        return max != 0;
    }

    virtual const void* getMin() const {
        return min;
    }

    virtual const void* getMax() const {
        return max;
    }

protected:
    static const int rank = N;
    const IVecN n;
    mutable ArrayN value, accumvalue, accumvalue2;
    double norm, accumnorm;
    bool hasError;
    int iblock;
    VecN* scale;
    VecN* origin;
    VecN* min;
    VecN* max;
};

template<int N>
BlitzArrayBlkdEst<N>::BlitzArrayBlkdEst(const std::string& name,
        const std::string& typeString, const IVecN& n, bool hasError) :
        ArrayEstimator(name, typeString, hasError), n(n), value(n), accumvalue(
                n), accumvalue2(hasError ? n : n * 0), norm(0), accumnorm(0), hasError(
                hasError), iblock(0), scale(0), origin(0), min(0), max(0) {
    value = 0;
    accumvalue = 0;
    accumvalue2 = 0;
}

template<int N>
void BlitzArrayBlkdEst<N>::averageOverClones(const MPIManager* mpi) {
#ifdef ENABLE_MPI
    if (mpi && mpi->isCloneMain()) {
        int rank = mpi->getCloneComm().Get_rank();
        int size = mpi->getCloneComm().Get_size();
        if (size>1) {
            reset();
            if (rank==0) {
#if MPI_VERSION==2
                mpi->getCloneComm().Reduce(MPI::IN_PLACE,&norm,1,MPI::DOUBLE,
                        MPI::SUM,0);
                mpi->getCloneComm().Reduce(MPI::IN_PLACE,value.data(),value.size(),
                        MPI::FLOAT,MPI::SUM,0);
#else
                double nbuff;
                ArrayN vbuff(n);
                mpi->getCloneComm().Reduce(&norm,&nbuff,1,MPI::DOUBLE,
                        MPI::SUM,0);
                mpi->getCloneComm().Reduce(value.data(),vbuff.data(),value.size(),
                        MPI::FLOAT,MPI::SUM,0);
                norm=nbuff;
                value=vbuff;
#endif
            } else {
                mpi->getCloneComm().Reduce(&norm,NULL,1,MPI::DOUBLE,MPI::SUM,0);
                mpi->getCloneComm().Reduce(value.data(),NULL,value.size(),
                        MPI::FLOAT,MPI::SUM,0);
            }
        }
    }
#endif
    // Next add value to accumvalue and accumvalue2.
    accumvalue += value / norm;
    if (hasErrorFlag)
        accumvalue2 += (value * value) / (norm * norm);
    accumnorm += 1.;
    value = 0.;
    norm = 0;
    ++iblock;
}

#endif
