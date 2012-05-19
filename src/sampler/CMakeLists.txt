set (sources
    AlternatingParticleChooser.cc
    CollectiveMover.cc
    CollectiveSectionMover.cc
    CollectiveSectionSampler.cc
    DampedFreeTensorMover.cc
    DisplaceMoveSampler.cc
    DoubleCollectiveSectionSampler.cc
    DoubleDisplaceMoveSampler.cc
    DoubleMLSampler.cc
    DoubleSectionChooser.cc
    ExchangeMover.cc
    FreeMover.cc
    FreeMoverPBC.cc
    FreeTensorMover.cc
    HyperbolicMover.cc
    ModelSampler.cc
    MiddleSectionChooser.cc
    MultiLevelSampler.cc
    MultiSpeciesParticleChooser.cc
    NonZeroSectionChooser.cc
    PairChooser.cc
    ParticleChooser.cc
    PermutationChooser.cc
    RandomPermutationChooser.cc
    SectionChooser.cc
    SimpleParticleChooser.cc
    SpeciesParticleChooser.cc
    SpinModelSampler.cc
    SpinStatePermutationChooser.cc
    TwoPairChooser.cc
    UniformMover.cc
    WalkingChooser.cc
)

add_library(sampler STATIC ${sources})

add_subdirectory(test)
