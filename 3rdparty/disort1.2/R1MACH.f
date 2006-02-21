c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: /srv/svn/cvs/cvsroot/arts/3rdparty/disort1.2/R1MACH.f,v 1.1 2006/02/21 16:23:28 olemke Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL FUNCTION R1MACH(I)

c        Single-precision machine constants

c  Assume floating-point numbers are represented in the t-digit,
c  base-b form

c         sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )

c  where 0.le.x(i).lt.b  for  i = 1,...,t,
c  0.lt.x(1), and  emin.LE.e.LE.emax.  then

c  R1MACH(1) = b**(emin-1), the smallest positive magnitude
c              (use TINY(R) in Fortran 90, where R is a single
c              precision variable)

c  R1MACH(2) = b**emax*(1 - b**(-t)), the largest magnitude
c              (use HUGE(R) in Fortran 90, where R is a single
c              precision variable))

c  R1MACH(3) = b**(-t), the smallest relative spacing.

c  R1MACH(4) = b**(1-t), the largest relative spacing.  i.e.,
c              smallest positive eps such that  1+eps .ne. 1
c              (use EPSILON(R) in Fortran 90, where R is a single
c              precision variable))

c  R1MACH(5) = LOG10(b)


c  Reference: Fox P.A., Hall A.D., Schryer N.L.,'Framework For A
c               Portable Library', ACM Transactions On Mathematical
c               Software, Vol. 4, No. 2, June 1978, pp. 177-188.


c  By default, returns values appropriate for a computer with IEEE 
c  arithmetic.  This is an abbreviated version of a routine widely
c  used for 20+ years by numerical analysts.  Most of the values in
c  the original version pertain to computers which went to computer
c  heaven years ago and are of little if any interest.
c 
c  If the values herein do not work for any reason, just look in
c  your Fortran manual for the correct values (usually in the part
c  discussing representations of numbers) and insert them. The exact
c  values are not that important; they can be a factor of 2-3 off
c  without causing any harm.

c  Only I = 1,2,4 is actually used by DISORT. 

c  This routine is superseded in Fortran-90 by the intrinsic numeric 
c  inquiry functions HUGE(1.0), TINY(1.0), and EPSILON(1.0).

c  The original version can be found on NetLib (search by name):
c      http://www.netlib.org/
c ====================================================================

      INTEGER I
      EXTERNAL  ERRMSG

      IF( I.EQ.1 )  THEN
         R1MACH = 1.2E-38
c        R1MACH = TINY(1.0)
      ELSE IF( I.EQ.2 )  THEN  
         R1MACH = 3.4E+38
c        R1MACH = HUGE(1.0)
      ELSE IF( I.EQ.4 )  THEN  
         R1MACH = 1.2E-07
c        R1MACH = EPSILON(1.0)
      ELSE
         CALL ERRMSG( 'R1MACH--argument incorrect', .TRUE.)
      END IF

      RETURN
      END

