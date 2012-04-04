add_subdirectory(demo)
add_subdirectory(action)
add_subdirectory(parser)
add_subdirectory(sampler)
add_subdirectory(emarate)
add_subdirectory(fixednode)
add_subdirectory(spin)
add_subdirectory(stats)

set (sources
    AngularMomentumEstimator.cc
    AperiodicGaussian.cc
    BeadFactory.cc
    BinProbDensity.cc
    BondLengthEstimator.cc
    BoxEstimator.cc
    Collect.cc
    ConditionalDensityGrid.cc
    ConductanceEstimator.cc
    ConductivityEstimator.cc
    ConductivityEstimator2D.cc
    CoulombEnergyEstimator.cc
    CountCountEstimator.cc
    CubicLattice.cc
    DensCountEstimator.cc
    DensDensEstimator.cc
    DensityCurrentEstimator.cc
    DensityEstimator.cc
    DiamagneticEstimator.cc
    DipoleMomentEstimator.cc
    DoubleParallelPaths.cc
    DynamicPCFEstimator.cc
    EIndEstimator.cc
    EnumeratedModelState.cc
    EwaldCoulombEstimator.cc
    EwaldSum.cc
    FreeEnergyEstimator.cc
    FrequencyEstimator.cc
    GridPotential.cc
    Hungarian.cc
    JEstimator.cc
    Measure.cc
    OptEwaldSum.cc
    PairDistance.cc
    PairIntegrator.cc
    PairPotential.cc
    ParallelPaths.cc
    PathReader.cc
    Paths.cc
    PeriodicGaussian.cc
    Permutation.cc
    PermutationEstimator.cc
    PositionEstimator.cc
    ProbDensityGrid.cc
    RandomNumGenerator.cc
    RingLattice.cc
    SHOPhase.cc
    SKOmegaEstimator.cc
    SeedRandom.cc
    SerialPaths.cc
    SimInfoWriter.cc
    SimulationInfo.cc
    Species.cc
    Spin4DPhase.cc
    SpinChargeEstimator.cc
    SpinModelState.cc
    StructReader.cc
    SuperCell.cc
    ThermoEnergyEstimator.cc
    TradEwaldSum.cc
    VIndEstimator.cc
    VirialEnergyEstimator.cc
    WindingEstimator.cc
    WireEwald.cc
    WorkerShifter.cc
    WritePaths.cc
    WriteProbDensity.cc
    XMLWriter.cc
)

add_executable(pi ${sources} main.cc)

TARGET_LINK_LIBRARIES(pi ${LIBXML2_LIBRARIES} )
target_link_libraries(pi ${BLAS_LIB})
target_link_libraries(pi ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
  target_link_libraries(pi ${LIBF2C_LIB})
endif(EXISTS ${LIBF2C_LIB})
target_link_libraries(pi ${FFTW3_LIB})
target_link_libraries(pi ${HDF5_LIB})
target_link_libraries(pi ${GSL_LIB})

set_target_properties(pi
    PROPERTIES
    LINKER_LANGUAGE CXX)

target_link_libraries(pi demo parser sampler emarate fixednode spin stats action)

#add_subdirectory(test)
INCLUDE_DIRECTORIES(${PI_BINARY_DIR} ${PI_SOURCE_DIR}/src)
