set (sources
    BeadFactory.cc
    DoubleParallelPaths.cc
    EnumeratedModelState.cc
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
