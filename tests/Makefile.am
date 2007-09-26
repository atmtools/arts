#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

SUBDIRS = 

TESTDIRS = Abs ClearSky DISORT DOIT MonteCarlo

TESTFILES = .arts .ac .xml .xml.bin

check: \
	TestMonteCarloDataPrepare \
	TestMonteCarloGeneral \
	TestMonteCarloGeneralGaussian \
	TestDOIT \
	TestClearSky \
	TestAbs

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

###############################################################

Test%:
	@if python testall.py $@ > $@.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.log"; exit 1; fi

EXTRA_DIST = testall.py $(foreach dir, $(TESTDIRS), $(foreach extension,$(TESTFILES),$(wildcard $(dir)/*$(extension))))

CLEANFILES = *.log */*xml.generated* */*.rep

