#!/bin/sh

ARTS=$1

if [ ! -x "$ARTS" ]
then
	echo "Error: Unable to execute arts" 1>&2
	exit 1
fi

METHOD_LIST=`$ARTS -m all | tail +5 | sed 's/^-//' | sed 's/^.-*.$//'`

echo "%------------------------------------------------------------"
echo "% This file has been generated automatically by the script"
echo "% arts_methods_to_latex.sh, which calls arts with the -d flag."
echo "%"
echo "% Generation date: `date`"
echo "% ARTS version:    `$ARTS -v | head -n 1`"
echo "%"
echo "% DO NOT EDIT!"
echo "%------------------------------------------------------------"
echo "%"

TMPFILE=/tmp/arts_parse.$$

for i in $METHOD_LIST
do
    $ARTS -d $i > $TMPFILE

    echo "Processing $i..." >&2
    INSIDE_VERB=no
    while read j
    do
	if [ "`echo $j | grep '^\*--*\*$'`" ]
	then 
	    echo "%$j"
	elif [ "`echo $j | grep '^Workspace method = '`" ]
	then
	    echo "$j" | sed 's/^Workspace method = \(.*\)$/\\levelb{\1}/' | sed 's/_/\\_/g'
	elif [ "`echo $j | grep '^--*$'`" ]
	then
	    if [ x$INSIDE_VERB = xno ]
	    then
		echo '\footnotesize\begin{verbatim}'
		INSIDE_VERB=yes
	    else
		echo "$j"
#		echo '\end{verbatim}'
#		INSIDE_VERB=no
	    fi
	elif [ "`echo $j | grep '^Types ='`" ]
	then
	    echo "$j"
	    echo '\end{verbatim}'
	    INSIDE_VERB=no
	else
	    echo "$j"
	fi

	
    done < $TMPFILE

    rm -f $TMPFILE
done

