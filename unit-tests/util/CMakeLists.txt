get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
    ${dir}/AperiodicGaussianTest.cc
    ${dir}/HungarianTest.cc
    ${dir}/PeriodicGaussianTest.cc
    ${dir}/PermutationTest.cc
    ${dir}/SuperCellTest.cc
    ${dir}/fft/FFT1DTest.cpp
    ${dir}/math/VPolyFitTest.cpp
    ${dir}/propagator/GridParametersTest.cpp
    ${dir}/propagator/KineticGridTest.cpp
    ${dir}/propagator/PotentialGridTest.cpp
    ${dir}/propagator/PropagatorGridTest.cpp
    ${dir}/propagator/Propagator1DTest.cpp
    PARENT_SCOPE
)
