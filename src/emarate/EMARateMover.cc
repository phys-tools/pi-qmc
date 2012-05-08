#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "EMARateMover.h"
#include "Beads.h"
#include "sampler/MultiLevelSampler.h"
#include "util/RandomNumGenerator.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include "util/SuperCell.h"
#include "SimulationInfo.h"
#include "util/PeriodicGaussian.h"
#include <cmath>

EMARateMover::EMARateMover(const SimulationInfo& simInfo, const int maxlevel,
        double C)
: PermutationChooser(2), ParticleChooser(2),
  lambda(simInfo.getNPart()), tau(simInfo.getTau()), C(C) {
    for (int i=0; i<simInfo.getNPart(); ++i) {
        const Species* species=&simInfo.getPartSpecies(i);
        lambda(i)=0.5/species->mass;
    }
}

EMARateMover::~EMARateMover() {
}

double EMARateMover::makeMove(MultiLevelSampler& sampler, const int level) {
    const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
    Beads<NDIM>& movingBeads=sampler.getMovingBeads();
    const SuperCell& cell=sampler.getSuperCell();
    const int nStride=(int)pow(2,level);
    const int nSlice=sectionBeads.getNSlice();
    const blitz::Array<int,1>& index=sampler.getMovingIndex();
    const int nMoving=index.size();
    blitz::Array<Vec,1> gaussRand(nMoving); gaussRand=0.0;
    double toldOverTnew=0;
    forwardProb = 0;
    // Make radiating vs. direct choice at highest level.
    //std::cout << nStride << ", " << nSlice << std::endl;
    if (nStride+1 == nSlice) {
        earlierTransitions = 0.;
        //std::cout << "Choosing radiating vs. diagonal." << std::endl;
        int i = index(0);
        int j = index(1);
        Vec re1 = movingBeads(0,0);
        Vec re2 = movingBeads(0,nSlice-1);
        Vec rh1 = movingBeads(1,0);
        Vec rh2 = movingBeads(1,nSlice-1);
        Vec deltae = re2-re1; cell.pbc(deltae);
        Vec deltah = rh2-rh1; cell.pbc(deltah);
        Vec delta1 = re1-rh1; cell.pbc(delta1);
        Vec delta2 = re2-rh2; cell.pbc(delta2);
        double inv2sigma2e = 0.5/(lambda(i)*tau*nStride);
        double inv2sigma2h = 0.5/(lambda(j)*tau*nStride);
        double inv2sigma21 = 1./((lambda(i)+lambda(j))*tau*nStride);
        double te = dot(deltae,deltae)*inv2sigma2e;
        double th = dot(deltah,deltah)*inv2sigma2h;
        double t1 = dot(delta1,delta1)*inv2sigma21;
        double t2 = dot(delta2,delta2)*inv2sigma21;
        //std::cout << t1 << ", " << t2 << ", " << te << ", " << th << std::endl;
        double pRad = exp(-(t1+t2));
        double pDiag = exp(-(te+th));
        //std::cout << pRad << ", " << pDiag << "; "
        //          << pDiag/(pDiag+pRad) << std::endl;
        if (RandomNumGenerator::getRand()>pDiag/(pDiag+C*pRad)) {
            //std::cout << "Trying radiating move." << std::endl;
            isSamplingRadiating=true;
        } else {
            //std::cout << "Trying diagonal move." << std::endl;
            isSamplingRadiating=false;
        }
    } {

        for (int islice=nStride; islice<nSlice-nStride; islice+=2*nStride) {
            RandomNumGenerator::makeGaussRand(gaussRand);
            if (isSamplingRadiating) { // NEED TO CODE AND TEST THIS
                int iMoving1 = 0;
                int iMoving2 = 1;
                int i1 = index(0);
                int i2 = index(1);
                Vec mass1 = 0.5/lambda(i1);
                Vec mass2 = 0.5/lambda(i2);
                //std::cout << islice-nStride << ", "
                //          << islice+nStride << ", " << nSlice/2 << std::endl;
                Vec prev1 = (islice-nStride <= nSlice/2)
                          ? movingBeads(iMoving1, islice-nStride)
                                  : movingBeads(iMoving2, nSlice-1-(islice-nStride));
                Vec prev2 = (islice-nStride < nSlice/2)
                          ? movingBeads(iMoving1, nSlice-1-(islice-nStride))
                                  : movingBeads(iMoving2, islice-nStride);
                Vec next1 = (islice+nStride <= nSlice/2)
                          ? movingBeads(iMoving1, islice+nStride)
                                  : movingBeads(iMoving2, nSlice-1-(islice+nStride));
                Vec next2 = (islice+nStride < nSlice/2)
                          ? movingBeads(iMoving1, nSlice-1-(islice+nStride))
                                  : movingBeads(iMoving2, islice+nStride);
                //std::cout << prev1 << next1 << std::endl;
                //std::cout << prev2 << next2 << std::endl;
                Vec midpoint1, midpoint2;
                // Should be rewritten to use PBC.
                if ((islice-nStride < nSlice/2) && (islice+nStride > nSlice/2)) {
                    midpoint1 = (mass1*prev1 + mass2*next1)/(mass1+mass2);
                    midpoint2 = (mass1*prev2 + mass2*next2)/(mass1+mass2);
                } else {
                    midpoint1 = 0.5*(prev1 + next1);
                    midpoint2 = 0.5*(prev2 + next2);
                }
                Vec sigma1, sigma2;
                if (islice-nStride < nSlice/2) {
                    if (islice+nStride > nSlice/2) {
                        sigma1 = sqrt(0.5*(lambda(i1)+lambda(i2))*tau*nStride);
                        sigma2 = sigma1;
                    } else {
                        sigma1 = sqrt(lambda(i1)*tau*nStride);
                        sigma2 = sqrt(lambda(i1)*tau*nStride);
                    }
                } else {
                    sigma1 = sqrt(lambda(i2)*tau*nStride);
                    sigma2 = sqrt(lambda(i2)*tau*nStride);
                }
                Vec delta1 = gaussRand(iMoving1) * sigma1;
                Vec delta2 = gaussRand(iMoving2) * sigma2;
                if (islice < nSlice/2) {
                    movingBeads(iMoving1,islice)=midpoint1+delta1;
                    movingBeads(iMoving1,nSlice-1-islice)=midpoint2+delta2;
                } else if (islice==nSlice/2) {
                    movingBeads(iMoving1,islice)=midpoint1+delta1;
                    movingBeads(iMoving2,islice)=midpoint2+delta2;
                } else {
                    movingBeads(iMoving2,nSlice-1-islice)=midpoint1+delta1;
                    movingBeads(iMoving2,islice)=midpoint2+delta2;
                }
            } else { // THIS IS CORRECT
                for (int iMoving=0; iMoving<nMoving; ++iMoving) {
                    const int i=index(iMoving);
                    double sigma = sqrt(lambda(i)*tau*nStride);
                    //double inv2Sigma2 = 0.5/(sigma*sigma);
                    // Calculate the new position.
                    Vec midpoint=movingBeads.delta(iMoving,islice+nStride,-2*nStride);
                    cell.pbc(midpoint)*=0.5;
                    midpoint+=movingBeads(iMoving,islice-nStride);
                    Vec delta = gaussRand(iMoving); delta*=sigma;
                    (movingBeads(iMoving,islice)=midpoint)+=delta;
                    cell.pbc(movingBeads(iMoving,islice));
                }
            }

        }

        // Calculate 1/transition probability for move

        // For now, we'll assume that the only two radiating particles are being moved.
        int iMoving1 = index(0);
        int iMoving2 = index(1);

        double oldDiagAction = 0.;
        double oldRadAction = 0.;
        double newDiagAction = 0.;
        double newRadAction = 0.;

        Vec mass1 = 0.5/lambda(iMoving1);
        Vec mass2 = 0.5/lambda(iMoving2);

        const Vec inv2Sigma21 = 0.5*mass1/(tau*nStride);
        const Vec inv2Sigma22 = 0.5*mass2/(tau*nStride);
        Vec rePrev = movingBeads(0,0);
        Vec rhPrev = movingBeads(1,0);
        Vec reRadPrev = rePrev;
        Vec rhRadPrev = rhPrev;
        Vec rePrevOld = sectionBeads(iMoving1,0);
        Vec rhPrevOld = sectionBeads(iMoving2,0);
        Vec reRadPrevOld = rePrevOld;
        Vec rhRadPrevOld = rhPrevOld;
        for (int islice=nStride; islice<nSlice; islice+=nStride) {

            // Calculate action for moving beads.
            Vec re = movingBeads(0,islice);
            Vec rh = movingBeads(1,islice);

            Vec reRad = re;
            Vec rhRad = (islice==nSlice/2)?re:rh;

            Vec delta=re-rePrev; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                newDiagAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
            }
            delta=rh-rhPrev; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                newDiagAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
            }
            delta=reRad-reRadPrev; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                newRadAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
            }
            delta=rhRad-rhRadPrev; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                newRadAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
            }
            // Calculate action for old beads.
            Vec reOld = sectionBeads(iMoving1,islice);
            Vec rhOld = sectionBeads(iMoving2,islice);

            Vec reRadOld = re;
            Vec rhRadOld = (islice==nSlice/2)?re:rh;

            delta=reOld-rePrevOld; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                oldDiagAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
            }
            delta=rhOld-rhPrevOld; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                oldDiagAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
            }
            delta=reRadOld-reRadPrevOld; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                oldRadAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
            }
            delta=rhRadOld-rhRadPrevOld; cell.pbc(delta);
            for (int idim=0;idim<NDIM;++idim) {
                oldRadAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
            }

            // Set the previous positions.
            rePrev = re;
            rhPrev = rh;
            reRadPrev = (islice==nSlice/2)?rhPrev:rePrev;
            rhRadPrev = rhPrev;
            rePrevOld = reOld;
            rhPrevOld = rhOld;
            reRadPrevOld = (islice==nSlice/2)?rhPrevOld:rePrevOld;
            rhRadPrevOld = rhPrevOld;
        }

        double oldAction = oldDiagAction-log(1+C*exp(-oldRadAction+oldDiagAction));
        double newAction = newDiagAction-log(1+C*exp(-newRadAction+newDiagAction));
        if (C > 0.) {
            if (log(C)-oldRadAction+newDiagAction > 40) {
                oldAction = -log(C)+oldRadAction;
            }
            if (log(C)-newRadAction+newDiagAction > 40) {
                newAction = -log(C)+newRadAction;
            }
        }

        toldOverTnew = newAction-oldAction;

    }
    double returnValue = toldOverTnew - earlierTransitions;
    earlierTransitions = toldOverTnew;

    return returnValue;
    //return toldOverTnew; //Return the log of the probability.
}


// Delayed Rejection 
double EMARateMover::makeDelayedMove(MultiLevelSampler& sampler,
        const int level) {
    return 1.0;
}
