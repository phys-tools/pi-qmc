#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PermutationChooser.h"
#include "util/Permutation.h"

PermutationChooser::PermutationChooser(int size) :
        permutation(new Permutation(size)) {
}

PermutationChooser::~PermutationChooser() {
    delete permutation;
}

bool PermutationChooser::choosePermutation() {
    return true;
}
