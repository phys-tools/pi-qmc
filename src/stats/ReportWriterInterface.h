#ifndef REPORTWRITERINTERFACE_H_
#define REPORTWRITERINTERFACE_H_

template <class T>
class ReportWriterInterface {
public:
    ReportWriterInterface() {}
    virtual ~ReportWriterInterface() {}

    virtual void startReport(const T& estimator) {}
    virtual void startBlock(int istep) {}
    virtual void reportStep(const T& estimator) {}
};

#endif
