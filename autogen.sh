#!/bin/sh
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# Remove the links:
rm -f INSTALL install-sh missing mkinstalldirs

rm -f config.cache
rm -f acconfig.h

echo "- aclocal."
aclocal -I m4

echo "- autoconf."
autoconf

# SAB 27.01.2000: This is part of the autotools package, so it 
#                 won´t work everywhere
#echo "- acconfig."
#acconfig

echo "- autoheader."
autoheader

echo "- automake."
automake --add-missing --copy --gnu


if test "`hostname | grep smiles`"; then
    ARTS_DATA_PATH="--with-arts-data=/pool/lookup2/arts-data"
else
    ARTS_DATA_PATH=
fi

./configure --enable-maintainer-mode $ARTS_DATA_PATH $@
exit
