#ifndef MPILIFECYCLE_H_
#define MPILIFECYCLE_H_

class MPILifecycle {
public:
    static void initialize(int argc, char **argv);
    static void finalize();
    static int getRank() {return rank;}
    static bool isEnabled() {return enabled;}
private:
    static bool enabled;
    static int rank;
};

#endif
