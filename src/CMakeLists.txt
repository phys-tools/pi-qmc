add_subdirectory(demo)
add_subdirectory(emarate)
add_subdirectory(fixednode)
add_subdirectory(parser)
add_subdirectory(sampler)
add_subdirectory(spin)
add_subdirectory(stats)

set (sources
    ActionChoice.cc
    AngularMomentumEstimator.cc
    AperiodicGaussian.cc
    BeadFactory.cc
    BinProbDensity.cc
    BondLengthEstimator.cc
    BoxEstimator.cc
    CaoBerneAction.cc
    Collect.cc
    CompositeAction.cc
    CompositeDoubleAction.cc
    ConditionalDensityGrid.cc
    ConductanceEstimator.cc
    ConductivityEstimator2D.cc
    ConductivityEstimator.cc
    CoulombAction.cc
    CoulombEnergyEstimator.cc
    CountCountEstimator.cc
    CubicLattice.cc
    DensCountEstimator.cc
    DensDensEstimator.cc
    DensityEstimator.cc
    DiamagneticEstimator.cc
    DipoleMomentEstimator.cc
    DotGeomAction.cc
    DoubleActionChoice.cc
    DoubleParallelPaths.cc
    DynamicPCFEstimator.cc
    EFieldAction.cc
    EIndEstimator.cc
    EnumeratedModelState.cc
    EwaldAction.cc
    EwaldCoulombEstimator.cc
    EwaldImagePairAction.cc
    EwaldSum.cc
    FixedPhaseAction.cc
    FreeEnergyEstimator.cc
    FrequencyEstimator.cc
    GateAction.cc
    GaussianAction.cc
    GaussianDotAction.cc
    GrapheneAction.cc
    GridPotential.cc
    Hungarian.cc
    HyperbolicAction.cc
    ImagePairAction.cc
    JelliumSlab.cc
    JEstimator.cc
    Measure.cc
    OptEwaldSum.cc
    OpticalLatticeAction.cc
    PairAction.cc
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
    PrimAnisSHOAction.cc
    PrimativePairAction.cc
    PrimColloidalAction.cc
    PrimCosineAction.cc
    PrimShellAction.cc
    PrimSHOAction.cc
    PrimTorusAction.cc
    ProbDensityGrid.cc
    QPCAction.cc
    RandomNumGenerator.cc
    RingGateAction.cc
    RingLattice.cc
    SeedRandom.cc
    SerialPaths.cc
    SHOAction.cc
    SHODotAction.cc
    SHOPhase.cc
    SimInfoWriter.cc
    SimulationInfo.cc
    SKOmegaEstimator.cc
    SmoothedGridPotential.cc
    Species.cc
    SphereAction.cc
    Spin4DPhase.cc
    SpinChargeEstimator.cc
    SpinModelState.cc
    SpringAction.cc
    SpringTensorAction.cc
    StillWebAction.cc
    StructReader.cc
    SuperCell.cc
    ThermoEnergyEstimator.cc
    TimpQPC.cc
    TradEwaldSum.cc
    TwoQDAction.cc
    VIndEstimator.cc
    VirialEnergyEstimator.cc
    WellImageAction.cc
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

target_link_libraries(pi demo emarate fixednode parser sampler spin stats)

#add_subdirectory(test)
INCLUDE_DIRECTORIES(${PI_BINARY_DIR} ${PI_SOURCE_DIR}/src)