#!/bin/sh

ARTS=$1

if [ ! -x "$ARTS" ]
then
    echo "Error: Unable to execute arts" 1>&2
    exit 1
fi

GROUP_LIST=`$ARTS -g | grep '^- ' | sed 's/^-//' | sed 's/^.-*.$//'`

echo "%------------------------------------------------------------"
echo "% This file has been generated automatically by the script"
echo "% arts_groups_to_latex.sh, which calls arts with the -g flag."
echo "%"
echo "% Generation date: `date`"
echo "% ARTS version:    `$ARTS -v | head -n 1`"
echo "%"
echo "% DO NOT EDIT!"
echo "%------------------------------------------------------------"
echo "%"

cat << EOF
\\begin{itemize}
EOF

for i in $GROUP_LIST
do
cat << EOF
\\item \\artsstyle{$i}
EOF
done

cat << EOF
\\end{itemize}
EOF

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
