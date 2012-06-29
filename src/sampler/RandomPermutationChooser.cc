#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "RandomPermutationChooser.h"
#include "util/RandomNumGenerator.h"
#include "util/Permutation.h"

RandomPermutationChooser::RandomPermutationChooser(const int nsize) :
        PermutationChooser(nsize), nsize(nsize) {
}

bool RandomPermutationChooser::choosePermutation() {
    for (int ifrom = 0; ifrom < nsize;) {
        int ito = (int) (nsize * RandomNumGenerator::getRand());
        if (ito == nsize)
            ito = nsize - 1;
        for (int jfrom = 0; jfrom <= ifrom; ++jfrom) {
            if (jfrom == ifrom)
                (*permutation)[ifrom++] = ito;
            if ((*permutation)[jfrom] == ito)
                break;
        }
    }
    return true;
}
