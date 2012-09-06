#ifndef REPORTWRITERINTERFACE_H_
#define REPORTWRITERINTERFACE_H_

template <class T>
class ReportWriterInterface {
public:
    ReportWriterInterface() {}
    virtual ~ReportWriterInterface() {}

    virtual void startScalarReport(const T& estimator) {}
    virtual void reportScalarStep(const T& estimator) {}
};

#endif
