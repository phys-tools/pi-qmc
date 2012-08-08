#include "Greeting.h"
#include "util/startup/MPILifecycle.h"
#include <config.h>
#include <iostream>

void Greeting::print() {
    int rank = MPILifecycle::getRank();;

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << "PIMC program: " << PACKAGE_STRING << "(";
        if (MPILifecycle::isEnabled()) {
            std::cout << "MPI enabled, ";
        } else {
            std::cout << "MPI disabled, ";
        }
#ifdef ENABLE_SPRNG
            std::cout << " SPRNG enabled, ";
#else
            std::cout << " SPRNG disabled, ";
#endif
            std::cout  << NDIM << "-d)" << std::endl;
        std::cout << std::endl;
    }
}
