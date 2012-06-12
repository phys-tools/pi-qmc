set (sources
    AnisotropicNodes.cc
    AtomicOrbitalDM.cc
    Atomic1sDM.cc
    Atomic2spDM.cc
    AugmentedNodes.cc
    ExcitonNodes.cc
    FixedNodeAction.cc
    FixedPhaseAction.cc
    FreeParticleNodes.cc
    FreePartNodesNoUpdate.cc
    GroundStateSNode.cc
    GroundStateWFNodes.cc
    NodeModel.cc
    NodeTester.cc
    SHONodes.cc
    SHOPhase.cc
    SHORealNodes.cc
    Spin4DPhase.cc
    SpinChoiceFixedNodeAction.cc
    WireNodes.cc
)

add_library(fixednode STATIC ${sources})

add_subdirectory(test)
