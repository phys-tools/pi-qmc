bin_PROGRAMS = unit-test-pi-qmc


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
    advancer/MultiLevelSamplerTest.cc \
    advancer/MultiLevelSamplerFake.cc \
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
    util/propagator/Propagator1DTest.cpp

unit_test_pi_qmc_LDADD = \
    -L$(top_builddir)/src/parser -lparser \
    -L$(top_builddir)/src/emarate -lemarate \
    -L$(top_builddir)/src/action -laction \
    -L$(top_builddir)/src/advancer -ladvancer \
    -L$(top_builddir)/src/algorithm -lalgorithm \
    -L$(top_builddir)/src/demo -ldemo \
    -L$(top_builddir)/src/estimator -lestimator \
    -L$(top_builddir)/src/fixednode -lfixednode \
    -L$(top_builddir)/src/spin -lspin \
    -L$(top_builddir)/src/stats -lstats \
    -L$(top_builddir)/src/base -lbase \
    -L$(top_builddir)/src/util -lutil \
    -L$(top_srcdir)/contrib/gtest-1.6.0 -lgtest -lpthread


all-local: $(top_builddir)/bin unit-test-pi-qmc gtest
	cp unit-test-pi-qmc $(top_builddir)/bin/

.PHONY: gtest

gtest:
	(cd $(top_srcdir)/contrib/gtest-1.6.0; cmake .; make)

unit_test_pi_qmc_DEPENDENCIES = \
    gtest \
    $(top_builddir)/src/action/libaction.la \
    $(top_builddir)/src/advancer/libadvancer.la \
    $(top_builddir)/src/base/libbase.la \
    $(top_builddir)/src/emarate/libemarate.la \
    $(top_builddir)/src/fixednode/libfixednode.la \
    $(top_builddir)/src/parser/libparser.la \
    $(top_builddir)/src/stats/libstats.la \
    $(top_builddir)/src/util/libutil.la
