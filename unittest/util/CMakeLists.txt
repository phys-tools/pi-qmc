set (dir util)
set_property(TARGET unittest APPEND PROPERTY sources
    ${dir}/AperiodicGaussianTest.cc
    ${dir}/HungarianTest.cc
    ${dir}/PeriodicGaussianTest.cc
    ${dir}/PermutationTest.cc
    ${dir}/SuperCellTest.cc
)
