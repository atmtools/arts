#!/usr/bin/env perl
#
# This script will look through all .h and .cc files for inclusion of plain
# C header files. These includes should be replaced by the corresponding
# C++ header.
#
# Author: Oliver Lemke <olemke@uni-bremen.de>
# Date:   2003-12-09

my @cheaders = (
    "assert.h",
    "ctype.h",
    "errno.h",
    "float.h",
    "limits.h",
    "locale.h",
    "math.h",
    "signal.h",
    "signal.h",
    "stdarg.h",
    "stddef.h",
    "stdio.h",
    "stdlib.h",
    "string.h",
    "time.h",
    "wchar.h",
    "wctype.h"
);

open (TPUT, "tput bold 2>/dev/null|");
while (<TPUT>) { $boldface = $_; }
close (TPUT);

open (TPUT, "tput sgr0 2>/dev/null|");
while (<TPUT>) { $normalface = $_; }
close (TPUT);

my @flist = <*.{cc,h}>;
my %headerstats;
foreach my $i (@flist) {
    my @includes;

    # Read include statements from sourcefile
    open (SOURCEFILE, "$i");
    while (<SOURCEFILE>) { if ($_ =~ "#include ") { push @includes, $_; } }
    close SOURCEFILE;

    # Loop over c header files and look for match in include list
    foreach my $j (@cheaders) {
        $headerstats->{$j}->{$i} = 0;
        foreach (@includes) {
            ($_ =~ "#include *<$j>") && $headerstats->{$j}->{$i}++;
        }
    }
}

my $cincludes = 0;
foreach my $header (keys %$headerstats) {
    my $first = 1;
    foreach my $source (keys %{$headerstats->{$header}}) {
        if ($headerstats->{$header}->{$source} > 0) {
            if ($first) {
                print "$boldface$header$normalface included by ";
                $first = 0;
            }
            print "$source ";
            ($headerstats->{$header}->{$source} > 1) && print "(multi) ";
            $cincludes += $headerstats->{$header}->{$source};
        }
    }
    (!$first) && print "\n";
}

if ($cincludes) {
    print "\n$cincludes plain C header files included.\n";
    print "These should be replaced with the C++ version.\n";
    die 1;
}

