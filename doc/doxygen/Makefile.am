#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

SUBDIRS        = . html # latex

# The latex docuement is huge, and seems not to be too useful.

EXTRA_DIST  = doxyfied

noinst_DATA = doxyfied

# The dependece on configure.in is there to make sure that the
# documentation is recompiled if the release number in that file is
# changed.

toexamine = \
$(top_srcdir)/configure.in \
$(top_srcdir)/src/*.h \
$(top_srcdir)/src/*.cc

if MAINTAINER_MODE

doxyfied: Doxyfile $(toexamine)
	doxygen Doxyfile > /dev/null
	doxytag -t arts.tag html
	cp ../arts-splash.gif html/
	touch doxyfied

clean:
	rm -f doxyfied

else

doxyfied:
	@echo
	@echo "=========================================================="
	@echo
	@echo "The documentation is only re-generated in maintainer mode."
	@echo "To activate maintainer-mode call configure with the option"
	@echo "--enable-maintainer-mode"
	@echo
	@echo "=========================================================="
	@echo

endif
