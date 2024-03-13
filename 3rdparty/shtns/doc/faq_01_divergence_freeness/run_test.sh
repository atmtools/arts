#! /bin/bash

NOHEADER=""
for i in 1 2 4 6 8 16 32 64; do
	NLATS=$((32*i))
	NLONS=$((NLATS*2))
	./div_free_test.py $NLATS $NLONS $NOHEADER nooutput
	NOHEADER="noheader"
done
