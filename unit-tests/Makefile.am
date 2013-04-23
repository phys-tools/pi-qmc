bin_PROGRAMS = unit-test-pi-qmc


unit_test_pi_qmc_LDADD = \
    -L$(top_srcdir)/contrib/gtest-1.6.0 -lgtest -lpthread \
    -L$(top_srcdir)/src/parser -lparser \
    -L$(top_srcdir)/src/emarate -lemarate \
    -L$(top_srcdir)/src/action -laction \
    -L$(top_srcdir)/src/advancer -ladvancer \
    -L$(top_srcdir)/src/algorithm -lalgorithm \
    -L$(top_srcdir)/src/demo -ldemo \
    -L$(top_srcdir)/src/estimator -lestimator \
    -L$(top_srcdir)/src/fixednode -lfixednode \
    -L$(top_srcdir)/src/spin -lspin \
    -L$(top_srcdir)/src/stats -lstats \
    -L$(top_srcdir)/src/util -lutil
    -L$(top_srcdir)/src/base -lbase \

unit_test_pi_qmc_CXXFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src \
    -I$(top_srcdir)/contrib/blitz-0.9 \
    -I$(top_srcdir)/contrib/gtest-1.6.0/include

unit_test_pi_qmc_SOURCES = \
    unittest_main.cc \
    action/DotGeomActionTest.cc \
    action/EFieldActionTest.cc \
    action/GateActionTest.cc \
    action/GaussianActionTest.cc \
    action/GaussianDotActionTest.cc \
    action/GrapheneActionTest.cc \
    action/PrimAnisSHOActionTest.cc \
    action/PrimColloidalActionTest.cc \
    action/PrimCosineActionTest.cc \
    action/PrimSHOActionTest.cc \
    action/SHOActionTest.cc \
    action/coulomb/Coulomb1DLinkActionTest.cpp \
    action/coulomb/Coulomb3DLinkActionTest.cpp \
    action/coulomb/CoulombLinkActionTest.cpp \
    action/interaction/SHOInteractionTest.cpp \
    advancer/CollectiveSectionMoverTest.cc \
    advancer/CollectiveSectionSamplerTest.cc \
    advancer/MultiLevelSamplerFake.cc \
    advancer/MultiLevelSamplerTest.cc \
    base/FermionWeightTest.cpp \
    emarate/EMARateActionTest.cc \
    emarate/EMARateEstimatorTest.cc \
    emarate/EMARateMoverTest.cc \
    emarate/EMARateTestBeadPositioner.cc \
    fixednode/Atomic1sDMTest.cc \
    fixednode/Atomic2spDMTest.cc \
    fixednode/AugmentedNodesTest.cc \
    parser/EstimatorParserTest.cc \
    stats/ScalarEstimatorTest.cpp \
    stats/SimpleScalarAccumulatorTest.cpp \
    stats/UnitsTest.cpp \
    util/AperiodicGaussianTest.cc \
    util/HungarianTest.cc \
    util/PeriodicGaussianTest.cc \
    util/PermutationTest.cc \
    util/SuperCellTest.cc \
    util/fft/FFT1DTest.cpp \
    util/math/VPolyFitTest.cpp \
    util/propagator/GridParametersTest.cpp \
    util/propagator/KineticGridTest.cpp \
    util/propagator/PotentialGridTest.cpp \
    util/propagator/PropagatorGridTest.cpp \
    util/propagator/PropagatorTest.cpp

all-local: $(top_builddir)/bin unit-test-pi-qmc gtest
	cp unit-test-pi-qmc $(top_builddir)/bin/

.PHONY: gtest

gtest:
	(cd $(top_srcdir)/contrib/gtest-1.6.0; cmake .; make)

unit_test_pi_qmc_DEPENDENCIES = \
    gtest \
    $(top_srcdir)/src/action/libaction.la \
    $(top_srcdir)/src/advancer/libadvancer.la \
    $(top_srcdir)/src/base/libbase.la \
    $(top_srcdir)/src/emarate/libemarate.la \
    $(top_srcdir)/src/fixednode/libfixednode.la \
    $(top_srcdir)/src/parser/libparser.la \
    $(top_srcdir)/src/stats/libstats.la \
    $(top_srcdir)/src/util/libutil.la
