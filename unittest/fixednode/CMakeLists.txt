set (sources
    AugmentedNodesTest.cc
    Atomic1sDMTest.cc
    Atomic2spDMTest.cc
)

add_library(fixednode_tests STATIC EXCLUDE_FROM_ALL ${sources})
