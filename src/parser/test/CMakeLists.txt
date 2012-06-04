set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/test)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
include_directories(${SRC_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(parserTests EXCLUDE_FROM_ALL
    EstimatorParserTest.cc
)
target_link_libraries(parserTests gtest gtest_main -lpthread)
target_link_libraries(parserTests ${BLAS_LIB})
target_link_libraries(parserTests ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(parserTests ${LIBF2C_LIB})
endif()

add_dependencies(unittest parserTests)
