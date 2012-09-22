set (sources
    BeadFactory.cc
    Charges.cpp
    DoubleParallelPaths.cc
    EnumeratedModelState.cc
    MagneticFluxCalculator.cpp
    MagneticFluxWeight.cpp
    ParallelPaths.cc
    Paths.cc
    SerialPaths.cc
    SimInfoWriter.cc
    SimulationInfo.cc
    Species.cc
    SpinModelState.cc
    XMLWriter.cc
)

add_library(base STATIC ${sources})
