set (sources
    AperiodicGaussian.cc
    EwaldSum.cc
    Hungarian.cc
    OptEwaldSum.cc
    PairDistance.cc
    PeriodicGaussian.cc
    Permutation.cc
    RandomNumGenerator.cc
    SuperCell.cc
    TradEwaldSum.cc
    WireEwald.cc
)

add_library(util STATIC ${sources})
