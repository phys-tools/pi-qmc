set (sources
    AperiodicGaussianTest.cc
    HungarianTest.cc
    PeriodicGaussianTest.cc
    PermutationTest.cc
    SuperCellTest.cc
)

add_library(util_tests STATIC EXCLUDE_FROM_ALL ${sources})
