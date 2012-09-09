set (sources
    AccRejEstimator.cc
    ArrayEstimator.cc
    AsciiReportBuilder.cc
    AsciiScalarReportWriter.cpp
    EstimatorIterator.cpp
    EstimatorReportBuilder.cpp
    EstimatorManager.cc
    H5ArrayReportWriter.cpp
    H5ReportBuilder.cc
    H5ScalarReportWriter.cpp
    MPIManager.cc
    ReportWriters.cpp
    ScalarEstimator.cc
    StdoutAccRejReportWriter.cpp
    StdoutArrayReportWriter.cpp
    StdoutScalarReportWriter.cpp
    StdoutReportBuilder.cc
    Units.cc
)

add_library(stats STATIC ${sources})
