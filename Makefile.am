SUBDIRS = src doc unit-tests
ACLOCAL_AMFLAGS = -I m4

.PHONY: test unit-test system-test docs

docs:
	(cd doc; make docs)

test: unit-test system-test

unit-test:
	echo "Running unit tests..."
	(cd contrib/gtest-1.6.0; cmake .; make)
	(cd unit-tests; make)
	bin/unit-test-pi-qmc

system-test: 
	echo "Running system integration tests..."
	bin/system-test-pi-qmc
