SUBDIRS = RIsearch1/src RIsearch2/src
EXTRA_DIST = README tables
install-data-hook:
	@mkdir -p $(DESTDIR)$(pkgdatadir)/tables
	cp -r tables/* $(DESTDIR)$(pkgdatadir)/tables

