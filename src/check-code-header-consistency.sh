#!/bin/sh

LOGFILE=check-code-header-consistency.log
echo > $LOGFILE

trap "rm -f test_include.cc; exit" EXIT 1 2 3 9 15

FAILED=0
EXPECTED_FAIL=0
FAILED_HEADERS=
for i in $SRCDIR/*.h
do
    echo -n "Compiling header $i... "
    echo "Compiling header $i" >> $LOGFILE
    echo >> $LOGFILE
    echo "#include \"$i\"" > test_include.cc
    echo "int main (void) {return 0;}" >> test_include.cc

    $COMPILE -o check-includes.o $i >> $LOGFILE 2>&1

    if [ $? -eq 0 ]; then
        echo "OK"
    else
        echo >> $LOGFILE
        case "$i" in
            $SRCDIR/xml_io_instantiation.h)
            echo "FAILED (expected)"
            EXPECTED_FAIL=1
            ;;
            *)
            echo "FAILED"
            FAILED=1
            FAILED_HEADERS="$FAILED_HEADERS $i"
            ;;
        esac
    fi
done

if [ $FAILED -eq 0 -a $EXPECTED_FAIL -eq 0 ]; then
    echo "Compilation of all headers successful"
elif [ $FAILED -eq 0 -a $EXPECTED_FAIL -eq 1 ]; then
    echo "Compilation of all headers except those expected to fail successful"
else
    echo "Compilation of the following headers failed:"
    echo $FAILED_HEADERS
    echo
    echo "See $LOGFILE for details"
    exit 1
fi

exit 0

