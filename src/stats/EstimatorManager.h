#ifndef __EstimatorManager_h_
#define __EstimatorManager_h_
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <hdf5.h>
#include "Estimator.h"
class MPIManager;
class EstimatorReportBuilder;
class EstimatorIterator;
class PartitionWeight;
class ScalarAccumulator;

/** Class for managing estimators. 
 Stores all estimators and sets of estimators.

 We define an Estimator to be any of the mathematical expressions
 that are sampled during the simulation. These can be physical
 quantities, such as ThermalEnergyEstimator or ColoumbEnergyEstimator.
 There are also Estimator classes for algorithmic quantities, such
 accept/reject ratios.

 The EstimatorManager is responsible for coordinating all the data
 collection and recording it to disk. We have chosen an HDF5
 file format (http://hdf.ncsa.uiuc.edu/HDF5) so that we can compress
 the data and include meta data and other structure.

 @todo Need to make public interface for access by EstimatorReportBuilder.
 @author John Shumway */
class EstimatorManager {
public:
    class SimInfoWriter {
    public:
        virtual void writeH5(hid_t) const {}
        virtual void writeStdout(std::ostream&) const {}
        virtual void writeAscii(std::ostream&) const {}
    };
    EstimatorManager(const std::string& filename, MPIManager *mpi,
            const SimInfoWriter*);
    std::vector<Estimator*>& getEstimatorSet(const std::string &name);
    virtual ~EstimatorManager();
    void add(Estimator* e);

    void startWritingGroup(const int nstep, const std::string& name);
    void writeStep();
    void recordInputDocument(const std::string &filename);
    EstimatorIterator getEstimatorIterator();

    int getNStep() const;
    void setPartitionWeight(PartitionWeight *partitionWeight);
    void setIsSplitOverStates(bool);
    bool getIsSplitOverStates();
    ScalarAccumulator* createScalarAccumulator();
private:
    typedef std::list<Estimator*> EstimatorList;
    typedef EstimatorList::iterator EstimatorIter;
    typedef std::list<EstimatorReportBuilder*> BuilderList;
    typedef BuilderList::iterator BuilderIter;
    const std::string filename;
    EstimatorList estimator;
    std::map<std::string, std::vector<Estimator*> > estimatorSet;
    int nstep;
    int istep;
    MPIManager *mpi;
    const SimInfoWriter *simInfoWriter;
    BuilderList builders;
    PartitionWeight *partitionWeight;
    int isSplitOverStates;
    void createBuilders(const std::string& filename,
            const SimInfoWriter* simInfoWriter);
};
#endif
