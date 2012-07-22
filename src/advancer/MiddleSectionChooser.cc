#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MiddleSectionChooser.h"
#include "action/Action.h"
#include "base/Paths.h"
#include "base/Beads.h"
#include "base/BeadFactory.h"
#include "util/Permutation.h"
#include "util/RandomNumGenerator.h"
#include <cmath>
#include <iostream>

MiddleSectionChooser::MiddleSectionChooser(const int nlevel, Paths &paths,
        Action &action, const BeadFactory &beadFactory) :
        SectionChooser(nlevel, paths.getNPart(), paths, action, beadFactory) {
    iFirstSlice = paths.getNSlice() - beads->getNSlice() / 2;
    std::cout << "Middle section chooser with iFirstSlice = " << iFirstSlice
            << std::endl;
}

MiddleSectionChooser::~MiddleSectionChooser() {
}

void MiddleSectionChooser::run() {
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
