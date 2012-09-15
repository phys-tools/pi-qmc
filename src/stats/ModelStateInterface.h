#ifndef MODELSTATEINTERFACE_H_
#define MODELSTATEINTERFACE_H_

class ModelStateInterface {
public:
    virtual ~ModelStateInterface() {}
    virtual int getModelCount() const=0;
    virtual int getModelState() const=0;
};

#endif
