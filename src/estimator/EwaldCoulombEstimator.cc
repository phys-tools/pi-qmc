#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "EwaldCoulombEstimator.h"
#include "action/Action.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "base/Paths.h"
#include "stats/ScalarAccumulator.h"
#include "util/OptEwaldSum.h"
#include "util/SuperCell.h"
#include "util/TradEwaldSum.h"
#include <blitz/tinyvec-et.h>

// constructor trad ewald
EwaldCoulombEstimator::EwaldCoulombEstimator(const SimulationInfo& simInfo,
        const Action* action, const double epsilon, const double rcut,
        const double kcut, const std::string& unitName, double scale,
        double shift, const double kappa, const int nImages,
        const bool testEwald, ScalarAccumulator *accumulator)
    :   ScalarEstimator("coulomb_energy", "scalar-energy/coulomb-energy",
                unitName, scale, shift),
        testEwald(testEwald),
        kappa(kappa),
        kcut(kcut),
        ewaldSum(*new TradEwaldSum(*simInfo.getSuperCell(), simInfo.getNPart(),
                        rcut, kcut, kappa)),
        cell(*simInfo.getSuperCell()),
        vgrid(1001),
        nradial(1001),
        rcut(rcut),
        dr(rcut / 1000),
        drinv(1. / dr),
        action(action),
        epsilon(epsilon),
        q(simInfo.getNPart()),
        r(simInfo.getNPart()),
        nImages(nImages),
        sphereR(0.) {
    this->accumulator = accumulator;
    for (int i = 0; i < q.size(); ++i)
        q(i) = simInfo.getPartSpecies(i).charge;
    ewaldSum.getQArray() = q;
    ewaldSum.evalSelfEnergy();
    vgrid(0) = -ewaldSum.evalFR0() / epsilon;
    for (int i = 1; i < nradial; ++i) {
        vgrid(i) = -ewaldSum.evalFR(i * dr) / epsilon;
    }

    // set up for summing over images for trad ewald sum
    for (int i = 0; i < NDIM; i++) {
        sphereR = (cell[i] > sphereR) ? cell[i] : sphereR;
    }
    sphereR *= nImages;
    boxImageVecs.resize(0);
    findBoxImageVectors(cell);
    std::cout << "In Ewald CEstimator Using Trad Ewald with sphereR: "
            << sphereR << std::endl;
}

// constructor opt ewald
EwaldCoulombEstimator::EwaldCoulombEstimator(const SimulationInfo& simInfo,
        const Action* action, const double epsilon, const double rcut,
        const double kcut, const std::string& unitName,
        double scale, double shift, ScalarAccumulator *accumulator)
    :   ScalarEstimator("coulomb_energy", "scalar-energy/coulomb-energy",
                unitName, scale, shift),
        testEwald(false),
        ewaldSum(*new OptEwaldSum(*simInfo.getSuperCell(), simInfo.getNPart(),
                        rcut, kcut, 4 * kcut, 8)),
        cell(*simInfo.getSuperCell()),
        vgrid(1001),
        nradial(1001),
        rcut(rcut),
        dr(rcut / 1000),
        drinv(1. / dr),
        action(action),
        epsilon(epsilon),
        q(simInfo.getNPart()),
        r(simInfo.getNPart()),
        nImages(0) {
    this->accumulator = accumulator;
    for (int i = 0; i < q.size(); ++i)
        q(i) = simInfo.getPartSpecies(i).charge;
    ewaldSum.getQArray() = q;
    ewaldSum.evalSelfEnergy();
    vgrid(0) = -ewaldSum.evalFR0() / epsilon;
    for (int i = 1; i < nradial; ++i) {
        vgrid(i) = -ewaldSum.evalFR(i * dr) / epsilon;
    }
}

EwaldCoulombEstimator::~EwaldCoulombEstimator() {
    delete &ewaldSum;
}

void EwaldCoulombEstimator::initCalc(const int nslice, const int firstSlice) {
    accumulator->clearValue();
}

void EwaldCoulombEstimator::handleLink(const Vec& start, const Vec& end,
        const int ipart, const int islice, const Paths& paths) {
    double energy = 0.0;
    if (nImages > 1) {
        for (int jpart = 0; jpart < ipart; ++jpart) {
            for (unsigned int img = 0; img < boxImageVecs.size(); img++) {
                Vec boxImage;
                for (int l = 0; l < NDIM; l++)
                    boxImage[l] = boxImageVecs[img][l]; // eventually change data strucure to use tinyvecs.

                Vec delta = end - paths(jpart, islice);
                cell.pbc(delta);
                double r = sqrt(dot(delta + boxImage, delta + boxImage));
                energy += q(ipart) * q(jpart)
                        * (1. / (r * epsilon) - ewaldSum.evalFR(r) / epsilon); /// could use the vgrid later after testing
            }
        }
    } else {
        for (int jpart = 0; jpart < ipart; ++jpart) {
            Vec delta = end - paths(jpart, islice);
            cell.pbc(delta);
            double r = sqrt(dot(delta, delta));
            if (r < rcut) {
                int igrid = (int) (r * drinv);
                double x = r - igrid * dr;
                energy += q(ipart) * q(jpart)
                        * (1. / (r * epsilon) + (1 - x) * vgrid(igrid)
                                + x * vgrid(igrid + 1));
            }
        }
    }

    // Add long range contribution.
    if (ipart == 0) {
        paths.getSlice(islice, r);
        energy += ewaldSum.evalLongRange(r) / epsilon;
    }

    accumulator->addToValue(energy);
}

void EwaldCoulombEstimator::endCalc(const int lnslice) {
    accumulator->storeValue(lnslice);
}

double EwaldCoulombEstimator::calcValue() {
    return 0.0;
}

void EwaldCoulombEstimator::reset() {
    accumulator->reset();
}

void EwaldCoulombEstimator::evaluate(const Paths& paths) {
    paths.sumOverLinks(*this);
}

void EwaldCoulombEstimator::findBoxImageVectors(const SuperCell &a) {
    // 3D case first... to be extended to 2D and 1D soon...after testing 3D case
    std::vector<std::vector<double> > vertices(8);
    for (unsigned int i = 0; i < vertices.size(); i++) {
        vertices[i].resize(NDIM);
    }
    vertices[0][0] = -a[0] / 2;
    vertices[0][1] = -a[1] / 2;
    vertices[0][2] = -a[2] / 2;
    vertices[1][0] = a[0] / 2;
    vertices[1][1] = -a[1] / 2;
    vertices[1][2] = -a[2] / 2;
    vertices[2][0] = a[0] / 2;
    vertices[2][1] = a[1] / 2;
    vertices[2][2] = -a[2] / 2;
    vertices[3][0] = -a[0] / 2;
    vertices[3][1] = a[1] / 2;
    vertices[3][2] = -a[2] / 2;
    vertices[4][0] = -a[0] / 2;
    vertices[4][1] = -a[1] / 2;
    vertices[4][2] = a[2] / 2;
    vertices[5][0] = a[0] / 2;
    vertices[5][1] = -a[1] / 2;
    vertices[5][2] = a[2] / 2;
    vertices[6][0] = a[0] / 2;
    vertices[6][1] = a[1] / 2;
    vertices[6][2] = a[2] / 2;
    vertices[7][0] = -a[0] / 2;
    vertices[7][1] = a[1] / 2;
    vertices[7][2] = a[2] / 2;

    // Check that all the box images are inside the sphere of radius
    // nImages*max(cell side) and save the location of these boxes
    // (center of box).
    std::vector<double> L(NDIM);
    for (int nx = -nImages; nx <= nImages; nx++) {
        for (int ny = -nImages; ny <= nImages; ny++) {
            for (int nz = -nImages; nz <= nImages; nz++) {

                L[0] = nx * a[0];
                L[1] = ny * a[1];
                L[2] = nz * a[2];
                int flag = 1;
                for (int v = 0; v < 8; v++) {
                    double tmpR = 0;
                    for (int i = 0; i < NDIM; i++) {
                        tmpR += (L[i] + vertices[v][i])
                                * (L[i] + vertices[v][i]);
                    }
                    if (sqrt(tmpR) > sphereR) {
                        flag = 0;
                        break;
                    }
                }
                if (flag == 1) {
                    boxImageVecs.push_back(L);
                }

            }
        }
    }
}
