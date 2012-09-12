#ifndef REPORTWRITERINTERFACE_H_
#define REPORTWRITERINTERFACE_H_

template <class T1, class T2>
class ReportWriterInterface {
public:
    ReportWriterInterface() {}
    virtual ~ReportWriterInterface() {}

    virtual void startReport(const T1* estimator, const T2* accumulator) {}
    virtual void startBlock(int istep) {}
    virtual void reportStep(const T1* estimator, const T2* accumulator) {}
};

#endif
