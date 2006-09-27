#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

SUBDIRS = data 
check:
	python testall.py -v

check-parallel: checkall

checkall: MonteCarloTests MiscTests

MonteCarloTests: MonteCarloTest  MonteCarloTest2 \
		 MonteCarloTest3 MonteCarloTest4

MiscTests: DoitTest ClearSkyTest AbsTest

DoitTest:
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi

ClearSkyTest:
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi

AbsTest:
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi

MonteCarloTest: MCDataPrepareTest
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi

MonteCarloTest2: MCDataPrepareTest
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi

MonteCarloTest3: MCDataPrepareTest
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi

MonteCarloTest4: MCDataPrepareTest
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi

MCDataPrepareTest:
	@if python testall.py $@ > $@.test.log 2>&1; then \
		echo "$@ ok."; \
		else echo "$@ failed. Find details in $@.test.log"; exit 1; fi


EXTRA_DIST = *.arts testall.py

CLEANFILES = *.rep *.test.log

