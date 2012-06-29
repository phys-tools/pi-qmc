get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources} 
    ${dir}/DampedFreeTensorMover.cc
    ${dir}/FreeMover.cc
    ${dir}/FreeMoverPBC.cc
    ${dir}/FreeTensorMover.cc
    ${dir}/HyperbolicMover.cc
    PARENT_SCOPE
)
