/* D1MACH.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static logical c_true = TRUE_;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* RCS version control information: */
/* $Header: /srv/svn/cvs/cvsroot/arts/src/disort_D1MACH.c,v 1.4 2006/02/20 10:18:34 olemke Exp $ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
doublereal d1mach_(integer *i__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

/*  Double-precision machine constants (see R1MACH for documentation). */
/*  By default, returns values appropriate for a computer with IEEE */
/*  arithmetic.  This is an abbreviated version of a routine widely */
/*  used for 20+ years by numerical analysts.  Most of the values in */
/*  the original version pertain to computers which went to computer */
/*  heaven years ago and are of little if any interest. */

/*  If the values herein do not work for any reason, just look in */
/*  your Fortran manual for the correct values (usually in the part */
/*  discussing representations of numbers) and insert them. The exact */
/*  values are not that important; they can be a factor of 2-3 off */
/*  without causing any harm. */
/*  Only I = 1,2,4 is actually used by DISORT. */
/*  This routine is superseded in Fortran-90 by the intrinsic numeric */
/*  inquiry functions HUGE(1.D0), TINY(1.D0), and EPSILON(1.D0). */
/*  The original version can be found on NetLib (search by name): */
/*      http://www.netlib.org/ */
/* ==================================================================== */
    if (*i__ == 1) {
	ret_val = 2.3e-308;
/*        D1MACH = TINY(1.D0) */
    } else if (*i__ == 2) {
	ret_val = 1.7e308;
/*        D1MACH = HUGE(1.D0) */
    } else if (*i__ == 4) {
	ret_val = 2.3e-16;
/*        D1MACH = EPSILON(1.D0) */
    } else {
	errmsg_("D1MACH--argument incorrect", &c_true, (ftnlen)26);
    }
    return ret_val;
} /* d1mach_ */

#ifdef __cplusplus
	}
#endif
