    set (sources
    ActionChoice.cc
    ActionParser.cc
    AlternatingParticleChooser.cc
    AngularMomentumEstimator.cc
    AnisotropicNodes.cc
    AperiodicGaussian.cc
    AugmentedNodes.cc
    BeadFactory.cc
    BinProbDensity.cc
    BondLengthEstimator.cc
    BoxEstimator.cc
    CaoBerneAction.cc
    Collect.cc
    CollectiveMover.cc
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
    DampedFreeTensorMover.cc
    DensCountEstimator.cc
    DensDensEstimator.cc
    DensityEstimator.cc
    DiamagneticEstimator.cc
    DipoleMomentEstimator.cc
    DisplaceMoveSampler.cc
    DotGeomAction.cc
    DoubleActionChoice.cc
    DoubleDisplaceMoveSampler.cc
    DoubleMLSampler.cc
    DoubleParallelPaths.cc
    DoubleSectionChooser.cc
    DynamicPCFEstimator.cc
    EFieldAction.cc
    EIndEstimator.cc
    EMARateAction.cc
    EMARateEstimator.cc
    EMARateMover.cc
    EnumeratedModelState.cc
    EstimatorParser.cc
    EwaldAction.cc
    EwaldCoulombEstimator.cc
    EwaldImagePairAction.cc
    EwaldSum.cc
    ExchangeMover.cc
    ExcitonNodes.cc
    FixedNodeAction.cc
    FixedPhaseAction.cc
    FreeEnergyEstimator.cc
    FreeMover.cc
    FreeMoverPBC.cc
    FreeParticleNodes.cc
    FreePartNodesNoUpdate.cc
    FreeTensorMover.cc
    FrequencyEstimator.cc
    GateAction.cc
    GaussianAction.cc
    GaussianDotAction.cc
    GrapheneAction.cc
    GridPotential.cc
    GroundStateSNode.cc
    GroundStateWFNodes.cc
    Hungarian.cc
    HyperbolicAction.cc
    HyperbolicMover.cc
    ImagePairAction.cc
    JelliumSlab.cc
    JEstimator.cc
    MainParser.cc
    Measure.cc
    MiddleSectionChooser.cc
    ModelSampler.cc
    MultiLevelSampler.cc
    MultiSpeciesParticleChooser.cc
    NodeModel.cc
    NodeTester.cc
    NonZeroSectionChooser.cc
    OptEwaldSum.cc
    OpticalLatticeAction.cc
    PairAction.cc
    PairChooser.cc
    PairDistance.cc
    PairIntegrator.cc
    PairPotential.cc
    ParallelPaths.cc
    ParticleChooser.cc
    PathReader.cc
    Paths.cc
    PeriodicGaussian.cc
    Permutation.cc
    PermutationChooser.cc
    PermutationEstimator.cc
    PIMCParser.cc
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
    RandomPermutationChooser.cc
    RingGateAction.cc
    RingLattice.cc
    SectionChooser.cc
    SeedRandom.cc
    SerialPaths.cc
    SHOAction.cc
    SHODotAction.cc
    SHONodes.cc
    SHOPhase.cc
    SHORealNodes.cc
    SimInfoParser.cc
    SimInfoWriter.cc
    SimpleParticleChooser.cc
    SimulationInfo.cc
    SKOmegaEstimator.cc
    SmoothedGridPotential.cc
    Species.cc
    SpeciesParticleChooser.cc
    SphereAction.cc
    Spin4DPhase.cc
    SpinChargeEstimator.cc
    SpinChoiceFixedNodeAction.cc
    SpinModelSampler.cc
    SpinModelState.cc
    SpinStatePermutationChooser.cc
    SpringAction.cc
    SpringTensorAction.cc
    StillWebAction.cc
    StructReader.cc
    SuperCell.cc
    ThermoEnergyEstimator.cc
    TimpQPC.cc
    TradEwaldSum.cc
    TwoPairChooser.cc
    TwoQDAction.cc
    UniformMover.cc
    VIndEstimator.cc
    VirialEnergyEstimator.cc
    WalkingChooser.cc
    WellImageAction.cc
    WindingEstimator.cc
    WireEwald.cc
    WireNodes.cc
    WorkerShifter.cc
    WritePaths.cc
    WriteProbDensity.cc
    XMLParser.cc
    XMLUnitParser.cc
    XMLWriter.cc
    assndx.f
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

target_link_libraries(pi demo spin stats)

#add_subdirectory(test)
INCLUDE_DIRECTORIES(${PI_BINARY_DIR} ${PI_SOURCE_DIR}/src)