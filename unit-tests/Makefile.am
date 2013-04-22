bin_PROGRAMS = unit-test-pi-qmc

unit_test_pi_qmc_CXXFLAGS = -I$(top_srcdir)/contrib/gtest-1.6.0/include

unit_test_pi_qmc_LDADD = -L$(top_srcdir)/contrib/gtest-1.6.0 -lgtest -lpthread

unit_test_pi_qmc_SOURCES = \
    unittest_main.cc

all-local: $(top_builddir)/bin unit-test-pi-qmc gtest
	cp unit-test-pi-qmc $(top_builddir)/bin/

.PHONY: gtest

gtest:
	(cd $(top_srcdir)/contrib/gtest-1.6.0; cmake .; make)

unit_test_pi_qmc_DEPENDENCIES = gtest
