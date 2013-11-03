set (sources
    FreeSpinMover.cc
    LatticeSpinPhase.cc
    MainSpinParser.cc
    SpinAction.cc
    SpinEstimator.cc
    SpinFixedPhaseAction.cc
    SpinMover.cc
    SpinPhase.cc
    SpinSetter.cc
)

add_library(spin STATIC ${sources})
