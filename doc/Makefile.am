docs:
	sed -e 's,top_srcdir,$(top_srcdir),' \
            -e 's/VERSION/$(VERSION)/' \
            < $(top_srcdir)/doc/etc/Doxyfile > Doxyfile
	doxygen Doxyfile
	@echo "Doxygen generated files are at doc/html/index.html"
.PHONY: docs
all:
	@echo "Run 'make docs' to generate Doxygen source code documentation."
clean:
	-rm -rf html
dist_noinst_DATA = etc/Doxyfile etc/doc.css etc/doxygen.sty etc/footer.html etc/header.html etc/header.tex etc/stylesheet.css
