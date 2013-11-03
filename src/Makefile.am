SUBDIRS = action advancer base demo emarate estimator fixednode parser \
	spin stats util algorithm

bin_PROGRAMS = pi-qmc

all-local: $(top_builddir)/bin pi-qmc
	cp pi-qmc $(top_builddir)/bin/

$(top_builddir)/bin:
	mkdir $(top_builddir)/bin

pi_qmc_CXXFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src \
    -I$(top_srcdir)/contrib/blitz-0.9

pi_qmc_DEPENDENCIES = \
    action/libaction.la \
    advancer/libadvancer.la \
    algorithm/libalgorithm.la \
    base/libbase.la \
    demo/libdemo.la \
    emarate/libemarate.la \
    estimator/libestimator.la \
    fixednode/libfixednode.la \
    parser/libparser.la \
    spin/libspin.la \
    stats/libstats.la \
    util/libutil.la
pi_qmc_LDADD = -Lparser -lparser -Lemarate -lemarate -Laction -laction \
    -Lalgorithm -lalgorithm -Ladvancer -ladvancer \
    -Lestimator -lestimator -Lfixednode -lfixednode \
    -Lspin -lspin -Lstats -lstats -Lbase -lbase -Lutil -lutil -Ldemo -ldemo
pi_qmc_SOURCES = \
    main.cc
