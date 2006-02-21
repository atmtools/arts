#!/bin/bash

DISORT_SOURCES="\
        D1MACH.f \
        DISORT.f \
        DISOTEST.f \
        ErrPack.f \
        LINPAK.f \
        R1MACH.f"

if [ "`basename \`pwd\``" != "disort1.2" ]; then
    echo "ERROR: This script must be run from the disort1.2 directory"
    exit 1
fi

ERROR=0
for i in $DISORT_SOURCES
do
    f2c -r8 -C++ -C $i
    if [ $? -ne 0 ]; then ERROR=1; break; fi
done

if [ $ERROR -ne 0 ]; then
    echo "ERROR: f2c failed." > /dev/stderr
    exit 1
fi

echo
echo "==========================================="
echo "Moving disort C files to arts/src directory"
echo "==========================================="
echo
for i in *.c
do
    mv -v $i ../../src/disort_$i
done

