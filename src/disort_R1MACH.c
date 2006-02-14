/* R1MACH.f -- translated by f2c (version 20050501).
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
/* $Header: /srv/svn/cvs/cvsroot/arts/src/disort_R1MACH.c,v 1.2 2006/02/14 15:41:17 olemke Exp $ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
doublereal r1mach_(integer *i__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

/*        Single-precision machine constants */
/*  Assume floating-point numbers are represented in the t-digit, */
/*  base-b form */
/*         sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) ) */
/*  where 0.le.x(i).lt.b  for  i = 1,...,t, */
/*  0.lt.x(1), and  emin.LE.e.LE.emax.  then */
/*  R1MACH(1) = b**(emin-1), the smallest positive magnitude */
/*              (use TINY(R) in Fortran 90, where R is a single */
/*              precision variable) */
/*  R1MACH(2) = b**emax*(1 - b**(-t)), the largest magnitude */
/*              (use HUGE(R) in Fortran 90, where R is a single */
/*              precision variable)) */
/*  R1MACH(3) = b**(-t), the smallest relative spacing. */
/*  R1MACH(4) = b**(1-t), the largest relative spacing.  i.e., */
/*              smallest positive eps such that  1+eps .ne. 1 */
/*              (use EPSILON(R) in Fortran 90, where R is a single */
/*              precision variable)) */
/*  R1MACH(5) = LOG10(b) */
/*  Reference: Fox P.A., Hall A.D., Schryer N.L.,'Framework For A */
/*               Portable Library', ACM Transactions On Mathematical */
/*               Software, Vol. 4, No. 2, June 1978, pp. 177-188. */
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
/*  inquiry functions HUGE(1.0), TINY(1.0), and EPSILON(1.0). */
/*  The original version can be found on NetLib (search by name): */
/*      http://www.netlib.org/ */
/* ==================================================================== */
    if (*i__ == 1) {
	ret_val = 1.2e-38;
/*        R1MACH = TINY(1.0) */
    } else if (*i__ == 2) {
	ret_val = 3.4e38;
/*        R1MACH = HUGE(1.0) */
    } else if (*i__ == 4) {
	ret_val = 1.2e-7;
/*        R1MACH = EPSILON(1.0) */
    } else {
	errmsg_("R1MACH--argument incorrect", &c_true, (ftnlen)26);
    }
    return ret_val;
} /* r1mach_ */

#ifdef __cplusplus
	}
#endif
