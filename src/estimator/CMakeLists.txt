set (sources
    AngularMomentumEstimator.cc
    BondLengthEstimator.cc
    BoxEstimator.cc
    ConductanceEstimator.cc
    ConductivityEstimator.cc
    ConductivityEstimator2D.cc
    CoulombEnergyEstimator.cc
    CountCountEstimator.cc
    DensCountEstimator.cc
    DensDensEstimator.cc
    DensityCurrentEstimator.cc
    DensityEstimator.cc
    DiamagneticEstimator.cc
    DipoleMomentEstimator.cc
    DynamicPCFEstimator.cc
    EIndEstimator.cc
    EwaldCoulombEstimator.cc
    FreeEnergyEstimator.cc
    FrequencyEstimator.cc
    JEstimator.cc
    PermutationEstimator.cc
    PositionEstimator.cc
    SKOmegaEstimator.cc
    SpinChargeEstimator.cc
    ThermoEnergyEstimator.cc
    VIndEstimator.cc
    VirialEnergyEstimator.cc
    WindingEstimator.cc
)

add_library(estimator STATIC ${sources})
