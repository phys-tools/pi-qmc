set (sources
    AccRejEstimator.cc
    ArrayEstimator.cc
    EstimatorIterator.cpp
    EstimatorReportBuilder.cpp
    EstimatorManager.cc
    MPIManager.cc
    ReportWriters.cpp
    SimpleScalarAccumulator.cpp
    ScalarEstimator.cc
    Units.cc
    ascii/AsciiReportBuilder.cc
    ascii/AsciiScalarReportWriter.cpp
    hdf5/H5ArrayReportWriter.cpp
    hdf5/H5Lib.cpp
    hdf5/H5ReportBuilder.cc
    hdf5/H5ScalarReportWriter.cpp
    hdf5/H5SplitScalarReportWriter.cpp
    stdout/StdoutAccRejReportWriter.cpp
    stdout/StdoutArrayReportWriter.cpp
    stdout/StdoutScalarReportWriter.cpp
    stdout/StdoutReportBuilder.cc
)

add_library(stats STATIC ${sources})
