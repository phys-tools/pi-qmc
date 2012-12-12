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
    erf.cpp
    ipow.cpp
    fft/FFT1D.cpp
    math/VPolyFit.cpp
    propagator/KineticGrid.cpp
    propagator/PotentialGrid.cpp
    propagator/Propagator.cpp
    propagator/PropagatorGrid.cpp
    startup/CommandLineParser.cpp
    startup/Greeting.cpp
    startup/MPILifecycle.cpp
    startup/UsageMessage.cpp
)

add_library(util STATIC ${sources})
