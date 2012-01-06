#include "CollectiveSectionMover.h"
#include <blitz/tinyvec-et.h>

CollectiveSectionMover::CollectiveSectionMover(double radius, Vec amplitude)
    :   radius(radius), amplitude(amplitude) {
}

CollectiveSectionMover::Vec CollectiveSectionMover::calcShift(const Vec & rin) const {
    double deltaR2 = dot(rin, rin);
    double g = 1 - deltaR2 / (radius*radius);
    if (g < 0) g = 0.;
    Vec rout = rin + g * amplitude;
    return rout;
}

CollectiveSectionMover::Vec CollectiveSectionMover::getAmplitude() const {
    return amplitude;
}

double CollectiveSectionMover::getRadius() const {
    return radius;
}

void CollectiveSectionMover::setAmplitude(Vec amplitude) {
    this->amplitude = amplitude;
}

void CollectiveSectionMover::setRadius(double radius) {
    this->radius = radius;
}



