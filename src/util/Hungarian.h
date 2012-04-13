#ifndef HUNGARIAN_H_
#define HUNGARIAN_H_

class Hungarian {
public:
    Hungarian(int size);
    ~Hungarian();
    int operator[](int index) const;
    void solve(double *matrix);
    double getSum() const;
private:
    int size;
    static const int MODE;
    int *kindex;
    double sum;
};

#endif
