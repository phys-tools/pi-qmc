#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "SectionChooser.h"
#include "action/Action.h"
#include "Paths.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "util/Permutation.h"
#include "util/RandomNumGenerator.h"
#include <cmath>
#include <iostream>
#include <gsl/gsl_qrng.h>

SectionChooser::SectionChooser(int nlevel, int npart, Paths &paths,
        Action &action, const BeadFactory &beadFactory) :
        CompositeAlgorithm(0), paths(&paths), action(&action), beads(
                beadFactory.getNewBeads(npart, (int) pow(2, nlevel) + 1)), permutation(
                new Permutation(npart)), nlevel(nlevel), qrng(
                gsl_qrng_alloc(gsl_qrng_sobol, 1)) {
}

SectionChooser::~SectionChooser() {
    delete beads;
    delete permutation;
    gsl_qrng_free(qrng);
}

void SectionChooser::run() {
    double x = RandomNumGenerator::getRand() * (1 - 1e-8);
    int ilo = paths->getLowestOwnedSlice(false) - 1;
    int ihi = paths->getHighestSampledSlice(beads->getNSlice() - 1, false);
    iFirstSlice = ilo + (int) ((ihi + 1 - ilo) * x);
    if (iFirstSlice > ihi)
        iFirstSlice = ihi;
    // Copy coordinates from allBeads to section Beads.
    paths->getBeads(iFirstSlice, *beads);
    permutation->reset();
    // Initialize the action.
    action->initialize(*this);
    // Run the sampling algorithm.
    CompositeAlgorithm::run();
    // Copy moved coordinates from sectionBeads to allBeads.
    paths->putBeads(iFirstSlice, *beads, *permutation);
    // Refresh the buffer slices.
    paths->setBuffers();
}
