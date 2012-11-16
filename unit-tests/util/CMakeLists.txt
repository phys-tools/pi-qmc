get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
    ${dir}/AperiodicGaussianTest.cc
    ${dir}/HungarianTest.cc
    ${dir}/PeriodicGaussianTest.cc
    ${dir}/PermutationTest.cc
    ${dir}/SuperCellTest.cc
    ${dir}/math/VPolyFitTest.cpp
    PARENT_SCOPE
)
