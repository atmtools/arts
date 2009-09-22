#!/bin/sh

STATE=0
LIST=

for i in $SRCDIR/*.h $SRCDIR/*.cc; do
    if grep '^ *		* *.*$' $i > /dev/null; then
        STATE=1
        LIST="${LIST}$i "
    fi
done

if [ "$LIST" ]; then
    echo
    echo "The following source files use tabs for indentation:"
    echo $LIST
fi

exit $STATE

