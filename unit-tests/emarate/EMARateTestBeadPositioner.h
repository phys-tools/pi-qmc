#ifndef EMARATETESTBEADPOSITIONER_H_
#define EMARATETESTBEADPOSITIONER_H_

class MultiLevelSamplerFake;

class EMARateTestBeadPositioner {
public:
    EMARateTestBeadPositioner(MultiLevelSamplerFake&);
    void setIdenticalPaths(double separation);
    void setRecombiningPaths(double separation);
private:
    MultiLevelSamplerFake& sampler;
};


#endif
