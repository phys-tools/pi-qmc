get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set_property(TARGET unittest APPEND PROPERTY sources
    ${dir}/AperiodicGaussianTest.cc
    ${dir}/HungarianTest.cc
    ${dir}/PeriodicGaussianTest.cc
    ${dir}/PermutationTest.cc
    ${dir}/SuperCellTest.cc
)
