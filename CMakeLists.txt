cmake_minimum_required (VERSION 2.6)

project (PI CXX)
set (PI_VERSION_MAJOR 1)
set (PI_VERSION_MINOR 4)

set(GMock_DIR ${PI_SOURCE_DIR}/CMakeModules)

configure_file (
  "${PI_SOURCE_DIR}/config.h.in"
  "${PI_BINARY_DIR}/config.h"
)

set (NDIM 3)

add_definitions(-DHAVE_CONFIG_H)
CONFIGURE_FILE(${PI_SOURCE_DIR}/CMakeConfig.h.in
               ${PI_BINARY_DIR}/config.h)

INCLUDE_DIRECTORIES(${PI_BINARY_DIR})


set(EXECUTABLE_OUTPUT_PATH ${PI_BINARY_DIR}/bin)

add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(src)
add_subdirectory(lib)
