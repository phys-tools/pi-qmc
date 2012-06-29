#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "NonZeroSectionChooser.h"
#include "action/Action.h"
#include "Paths.h"
#include "Beads.h"
#include "util/Permutation.h"
#include "util/RandomNumGenerator.h"
#include <cmath>
#include <iostream>

NonZeroSectionChooser::NonZeroSectionChooser(const int nlevel, Paths &paths,
        Action &action, const BeadFactory &beadFactory) :
        SectionChooser(nlevel, paths.getNPart(), paths, action, beadFactory) {
}

NonZeroSectionChooser::~NonZeroSectionChooser() {
}

void NonZeroSectionChooser::run() {
    double x = RandomNumGenerator::getRand() * (1 - 1e-8);
    int ilo = paths->getLowestOwnedSlice(false) - 1;
    if (ilo < 1)
        ilo = 1;
    int ihi = paths->getHighestSampledSlice(beads->getNSlice() - 1, false);
    if (ihi + beads->getNSlice() >= paths->getNSlice()) {
        ihi = paths->getNSlice() - beads->getNSlice();
    }
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
