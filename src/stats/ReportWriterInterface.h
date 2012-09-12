#ifndef REPORTWRITERINTERFACE_H_
#define REPORTWRITERINTERFACE_H_

template <class T1, class T2>
class ReportWriterInterface {
public:
    ReportWriterInterface() {}
    virtual ~ReportWriterInterface() {}

    virtual void startReport(const T1* estimator, const T2* accumulator) = 0;
    virtual void startBlock(int istep) = 0;
    virtual void reportStep(const T1* estimator, const T2* accumulator) = 0;
};

#endif
