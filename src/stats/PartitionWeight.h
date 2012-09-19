#ifndef PARTITIONWEIGHT_H_
#define PARTITIONWEIGHT_H_

class PartitionWeight {
public:
    virtual ~PartitionWeight() {}
    virtual int getPartitionCount() const=0;
    virtual double getValue(int i) const=0;
};

#endif
