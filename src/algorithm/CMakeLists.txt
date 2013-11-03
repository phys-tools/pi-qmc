set (sources
    BinProbDensity.cc
    Collect.cc
    ConditionalDensityGrid.cc
    CubicLattice.cc
    Measure.cc
    PathReader.cc
    ProbDensityGrid.cc
    RingLattice.cc
    SeedRandom.cc
    StructReader.cc
    WorkerShifter.cc
    WritePaths.cc
    WriteProbDensity.cc
)

add_library(algorithm STATIC ${sources})
