#!/bin/sh

ARTS=$1

if [ ! -x "$ARTS" ]
then
    echo "Error: Unable to execute arts" 1>&2
    exit 1
fi

METHOD_LIST=`$ARTS -m all | grep '^- ' | sed 's/^-//' | sed 's/^.-*.$//' | sort -f`

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
    INSIDE_DESC=no
    while read j
    do
        case "$j" in
         \*--*\*)
             echo "%$j"
             ;;
         "Workspace method = "*)
             echo "$j" \
              | sed "s/^Workspace method = \(.*\)$/\\\section{\1}/" \
              | sed "s/_/\\\_/g"
             ;;
         --------------------------------------------*)
             if [ x$INSIDE_DESC = xno ]
             then
                 cat << EOF
{\\footnotesize\\tt\\hbadness=10000
EOF
#\\footnotesize\\begin{verbatim}
                 INSIDE_DESC=yes
             else
                 echo "$j"
                 cat << EOF
\vspace{1em}
EOF
             fi
             ;;
         "Types ="*)
             echo "\\noindent $j" \
              | sed "s/_/\\\_/g"
             echo "}"
             #echo "\\end{verbatim}"
             INSIDE_DESC=no
             ;;
         *)
             if [ -z "$j" ]; then
                 cat << EOF
\vspace{1em}
EOF
             else
                 echo "\\noindent $j " \
                  | sed "s/_/\\\_/g" \
                  | sed "s/\#/\\\#/g" \
                  | sed "s/\^/\\\^/g"
                 echo
             fi
             ;;
        esac

    done < $TMPFILE

    rm -f $TMPFILE
done

# Local variables:
# mode: ksh
# ksh-indent: 4
# ksh-group-offset: -4
# ksh-brace-offset: -2
# ksh-case-item-offset: 1
# ksh-case-indent: 4
# ksh-tab-always-indent: nil
# indent-tabs-mode: nil
# End:
