#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

SUBDIRS = 

TESTDIRS = Abs AMSU ClearSky DISORT DOIT MonteCarlo OdinSMR

TESTFILES = .arts .ac .xml .xml.gz .xml.bin

check: \
	TestMonteCarloDataPrepare \
	TestMonteCarloGeneral \
	TestMonteCarloGeneralGaussian \
	TestRteCalcMC \
	TestOdinSMR \
	TestDOIT \
	TestClearSky \
	TestAbs \
	TestAMSUB \
	TestAMSUB_fast

# Broken tests:
#	TestMonteCarloSimple
#	TestMonteCarloWithIncomingLookup
#	TestDISORT

#
# Dependencies between tests
#
TestMonteCarloSimple: TestMonteCarloDataPrepare

TestMonteCarloGeneral: TestMonteCarloDataPrepare

TestMonteCarloGeneralGaussian: TestMonteCarloDataPrepare

TestMonteCarloWithIncomingLookup: TestMonteCarloDataPrepare

TestRteCalcMC: TestMonteCarloDataPrepare

TestAMSUB_fast: TestAMSUB

###############################################################

EXTRA_DIST = $(srcdir)/artsXML.py  $(srcdir)/testall.py $(foreach dir, $(TESTDIRS), $(foreach extension,$(TESTFILES),$(wildcard $(srcdir)/$(dir)/*$(extension))))

CLEANFILES = $(wildcard *.log) $(wildcard $(builddir)/*/*.generated.*.xml*) $(wildcard $(builddir)/*/*xml.generated*) $(wildcard $(builddir)/*/*.rep)

clean-local:
	@if [ "$(top_builddir)" != "$(top_srcdir)" ]; then \
		rm -f artsXML.py artsXML.pyc testall.py; \
		rm -rf $(TESTDIRS); \
	fi

Test%:
	@if [ "$(top_builddir)" != "$(top_srcdir)" ]; then \
		for i in $(TESTDIRS); do \
			test ! -e "$$i" && mkdir "$$i"; \
		done; \
		for i in  $(EXTRA_DIST); do \
			TESTDIR=$$(basename $$(dirname $$i)); \
			test "$$TESTDIR" = "tests" && TESTDIR=.; \
			cp -u "$$i" "$$TESTDIR/"; \
		done; \
	fi

	@if TOPSRCDIR=$(top_srcdir) python testall.py $@ > $@.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.log"; exit 1; fi


