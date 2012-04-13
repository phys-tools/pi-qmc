set (sources
    AnisotropicNodes.cc
    AtomicOrbitalDM.cc
    Atomic1sDM.cc
    Atomic2spDM.cc
    AugmentedNodes.cc
    ExcitonNodes.cc
    FixedNodeAction.cc
    FreeParticleNodes.cc
    FreePartNodesNoUpdate.cc
    GroundStateSNode.cc
    GroundStateWFNodes.cc
    NodeModel.cc
    NodeTester.cc
    SHONodes.cc
    SHORealNodes.cc
    SpinChoiceFixedNodeAction.cc
    WireNodes.cc
)

add_library(fixednode STATIC ${sources})

add_subdirectory(test)
