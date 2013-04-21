SUBDIRS = src doc
ACLOCAL_AMFLAGS = -I m4
docs:
	(cd doc; make docs)
.PHONY: test

test: all
	(cd test/system; python run_tests.py)
