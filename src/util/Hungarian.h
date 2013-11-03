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
    int *work;
    double sum;
    void assndx(const int *mode, double *a, int *n,
        int *m, int *ida, int *k, double *sum, int *iw,
        int *idw);
};

#endif
