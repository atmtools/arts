/* ErrPack.f -- translated by f2c (version 20050501).
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

static integer c__1 = 1;
static logical c_true = TRUE_;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* RCS version control information: */
/* $Header: /srv/svn/cvs/cvsroot/arts/src/disort_ErrPack.c,v 1.2 2006/02/14 15:41:17 olemke Exp $ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* Subroutine */ int errmsg_(char *messag, logical *fatal, ftnlen messag_len)
{
    /* Initialized data */

    static integer nummsg = 0;
    static integer maxmsg = 100;
    static logical msglim = FALSE_;

    /* Format strings */
    static char fmt_99[] = "(//,\002 >>>>>>  TOO MANY WARNING MESSAGES -- \
 \002,\002They will no longer be printed  <<<<<<<\002,//)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, "(/,2A,/)", 0 };
    static cilist io___5 = { 0, 6, 0, "(/,2A,/)", 0 };
    static cilist io___6 = { 0, 6, 0, fmt_99, 0 };


/*        Print out a warning or error message;  abort if error */
    if (*fatal) {
	s_wsfe(&io___4);
	do_fio(&c__1, " ******* ERROR >>>>>>  ", (ftnlen)23);
	do_fio(&c__1, messag, messag_len);
	e_wsfe();
	s_stop("", (ftnlen)0);
    }
    ++nummsg;
    if (msglim) {
	return 0;
    }
    if (nummsg <= maxmsg) {
	s_wsfe(&io___5);
	do_fio(&c__1, " ******* WARNING >>>>>>  ", (ftnlen)25);
	do_fio(&c__1, messag, messag_len);
	e_wsfe();
    } else {
	s_wsfe(&io___6);
	e_wsfe();
	msglim = TRUE_;
    }
    return 0;
} /* errmsg_ */

logical wrtbad_(char *varnam, ftnlen varnam_len)
{
    /* Initialized data */

    static integer nummsg = 0;
    static integer maxmsg = 50;

    /* System generated locals */
    logical ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 6, 0, "(3A)", 0 };


/*          Write names of erroneous variables and return 'TRUE' */
/*      INPUT :   VarNam = Name of erroneous variable to be written */
/*                         ( CHARACTER, any length ) */
    ret_val = TRUE_;
    ++nummsg;
    s_wsfe(&io___9);
    do_fio(&c__1, " ****  Input variable  ", (ftnlen)23);
    do_fio(&c__1, varnam, varnam_len);
    do_fio(&c__1, "  in error  ****", (ftnlen)16);
    e_wsfe();
    if (nummsg == maxmsg) {
	errmsg_("Too many input errors.  Aborting...", &c_true, (ftnlen)35);
    }
    return ret_val;
} /* wrtbad_ */

logical wrtdim_(char *dimnam, integer *minval, ftnlen dimnam_len)
{
    /* System generated locals */
    logical ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, "(/,3A,I7)", 0 };


/*          Write name of too-small symbolic dimension and */
/*          the value it should be increased to;  return 'TRUE' */
/*      INPUT :  DimNam = Name of symbolic dimension which is too small */
/*                        ( CHARACTER, any length ) */
/*               Minval = Value to which that dimension should be */
/*                        increased (at least) */
    s_wsfe(&io___10);
    do_fio(&c__1, " ****  Symbolic dimension  ", (ftnlen)27);
    do_fio(&c__1, dimnam, dimnam_len);
    do_fio(&c__1, "  should be increased to at least ", (ftnlen)34);
    do_fio(&c__1, (char *)&(*minval), (ftnlen)sizeof(integer));
    e_wsfe();
    ret_val = TRUE_;
    return ret_val;
} /* wrtdim_ */

logical tstbad_(char *varnam, doublereal *relerr, ftnlen varnam_len)
{
    /* System generated locals */
    doublereal d__1;
    logical ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, "(/,3A,1P,E11.2,A)", 0 };


/*       Write name (VarNam) of variable failing self-test and its */
/*       percent error from the correct value;  return  'FALSE'. */
    ret_val = FALSE_;
    s_wsfe(&io___11);
    do_fio(&c__1, " Output variable ", (ftnlen)17);
    do_fio(&c__1, varnam, varnam_len);
    do_fio(&c__1, " differed by ", (ftnlen)13);
    d__1 = *relerr * 100.;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, " per cent from correct value.  Self-test failed.", (ftnlen)
	    48);
    e_wsfe();
    return ret_val;
} /* tstbad_ */

#ifdef __cplusplus
	}
#endif
