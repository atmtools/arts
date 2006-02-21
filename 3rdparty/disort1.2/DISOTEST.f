
      PROGRAM  TESTDO

c    Runs test problems for DISORT and checks answers. These
c    problems test almost all logical branches in DISORT.

c    As distributed, runs all problems.  To run just a subset,
c    change DATA statement for DOPROB;  e.g. to run just problem 4,
c    set DOPROB(4) TRUE and all others FALSE.

c    It is HIGHLY recommended that you use the code below as a template
c    for creating your own CALLs to DISORT, rather than starting from
c    scratch.  This will prevent a lot of mistakes and ensure that every
c    input argument gets a value.  Note in particular how GETMOM is
c    sometimes called to fill an array section of PMOM (for one layer);
c    several people have done this incorrectly in attempting to write it
c    ab initio (passing array sections for arrays that do not start at 
c    element 1 is tricky).

c    Note that the ratio to the 'correct answer' may occasionally be
c    significantly different from unity -- even so different that
c    the ratio just prints as ****** rather than a number.  However,
c    this mostly occurs for values of flux or intensity that are very
c    small compared to the forcing functions (that is, small compared
c    to internal thermal emission and/or radiation incident at the 
c    boundaries).  The printed number 'SERIOUSLY NON-UNIT RATIOS'
c    attempts to count just the cases where there is a real disagreement
c    and not those where quantitites are down at their noise level 
c    (defined as 10^(-6) times their maximum value).

c    Further documentation can be found in the file DISOTEST.doc.


c  Routines called :

c    DISORT:   The discrete ordinates radiative transfer program

c    GETMOM:   Sets phase function Legendre coefficients

c    PRTFIN:   Prints fluxes and intensities and their ratios to
c              the correct values

c    CHEKDO:   Data block containing correct fluxes and intensities

c    RATIO :   Ratio of calculated to correct value with underflow
c              and overflow protection (kept in file DISORT.f)

c+---------------------------------------------------------------------+

C     IMPLICIT  NONE

c                 ** DISORT I/O specifications **

      INTEGER  MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU
      PARAMETER ( MAXCLY = 6, MAXCMU = 48, MAXPHI = 3, MAXULV = 5,
     &            MAXUMU = 10 )
      CHARACTER  HEADER*127
      LOGICAL  DELTAM, LAMBER, PLANK, ONLYFL, PRNT(7), USRANG, USRTAU
      INTEGER  IBCND, NLYR, NUMU, NSTR, NPHI, NTAU
      REAL     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
     &         HL( 0:MAXCMU ), PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ),
     &         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     &         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV )

      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     &         DFDT( MAXULV ), UAVG( MAXULV ), U0U( MAXUMU, MAXULV ),
     &         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     &         TRNMED( MAXUMU )

c+---------------------------------------------------------------------+

c                     ** Correct answers **

      INTEGER  MXPROB, MXCASE, MAXTAU, MAXMU, MAXAZ
      PARAMETER  ( MXPROB = 9, MXCASE = 8, MAXTAU = 5, MAXMU = 10,
     &             MAXAZ = 3 )
      REAL  TSTFIR( MAXTAU, MXCASE, MXPROB ),
     &      TSTFDN( MAXTAU, MXCASE, MXPROB ),
     &      TSTFUP( MAXTAU, MXCASE, MXPROB ),
     &      TSTDFD( MAXTAU, MXCASE, MXPROB ),
     &      TSTUU ( MAXTAU, MAXMU, MAXAZ, MXCASE, MXPROB )
      COMMON / DOCHEK / TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU

c+---------------------------------------------------------------------+

      INTEGER  MXTAU, MXMU, MXPHI
      PARAMETER     ( MXTAU = 5, MXMU = 10, MXPHI = 3 )
      CHARACTER*1   ABC(18)*1, TITLE*100, BLANKS*3
      LOGICAL       AZMAVG, DOPROB( 13 )
      INTEGER       ICAS, IOD, ISS, IU, J, K, LC, LENTIT, LU, NPROB
      REAL          CMPFIR( MXTAU ), CMPFDN( MXTAU ), CMPFUP( MXTAU ),
     &              CMPDFD( MXTAU ), CMPUU ( MXTAU, MXMU, MXPHI ),
     &              PI

c     .. External Subroutines ..

      EXTERNAL  DISORT, ERRMSG, GETMOM, PRTFIN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, FLOAT, INDEX
c     ..
      DATA  ABC / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
     &            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r' /,
     &      PRNT / .TRUE., 5*.FALSE., .TRUE. /, ACCUR / 0.0 /,
     &      BLANKS / '   ' /,  DOPROB / 13*.TRUE. /


      PI = 2.* ASIN( 1.0 )

      IF( DOPROB(1) )  THEN

c **********************************************************************
c ****  Test Problem 1:  Isotropic Scattering                       ****
c ****  (Compare to Ref. VH1, Table 12)                             ****
c **********************************************************************

      NSTR = 16
      NLYR = 1
      CALL  GETMOM( 1, 0.0, NSTR, PMOM )
      USRTAU    = .TRUE.
      NTAU      = 2
      UTAU( 1 ) = 0.0
      USRANG    = .TRUE.
      NUMU      = 6
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.5
      UMU( 3 )  = -0.1
      UMU( 4 )  =  0.1
      UMU( 5 )  =  0.5
      UMU( 6 )  =  1.0
      NPHI      = 1
      PHI( 1 )  = 0.0
      IBCND     = 0
      UMU0      = 0.1
      PHI0      = 0.0
      LAMBER    = .TRUE.
      ALBEDO    = 0.0
      DELTAM    = .FALSE.
      PLANK     = .FALSE.
      ONLYFL    = .FALSE.

      DO 10  ICAS = 1, 6

         IF ( ICAS.EQ.1 ) THEN

            UTAU( 2 )  = 0.03125
            SSALB( 1 ) = 0.2
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.2 ) THEN

            UTAU( 2 )  = 0.03125
            SSALB( 1 ) = 1.0
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.3 ) THEN

            UTAU( 2 )  = 0.03125
            SSALB( 1 ) = 0.99
            FBEAM      = 0.0
            FISOT      = 1.0

         ELSE IF ( ICAS.EQ.4 ) THEN

            UTAU( 2 )  = 32.0
            SSALB( 1 ) = 0.2
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.5 ) THEN

            UTAU( 2 )  = 32.0
            SSALB( 1 ) = 1.0
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.6 ) THEN

            UTAU( 2 )  = 32.0
            SSALB( 1 ) = 0.99
            FBEAM      = 0.0
            FISOT      = 1.0

         END IF

         DTAUC( 1 ) = UTAU( 2 )

         WRITE( HEADER,'(3A,F9.5,A,F5.2)') 'Test Case No. 1',ABC(ICAS),
     &          ':  Isotropic Scattering, Ref. VH1, Table 12:  b =',
     &          UTAU(2), ', a =', SSALB(1)

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         AZMAVG = .TRUE.
         NPROB = 1
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &       ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
  10  CONTINUE

      ENDIF


      IF( DOPROB(2) )  THEN

c **********************************************************************
c ****  Test Problem 2:  Rayleigh Scattering, Beam Source           ****
c ****  (Compare To Ref. SW, Table 1)                               ****
c **********************************************************************

      NSTR = 16
      NLYR = 1
      CALL  GETMOM( 2, 0.0, NSTR, PMOM )
      USRTAU    = .TRUE.
      NTAU      = 2
      UTAU( 1 ) = 0.0
      USRANG    = .TRUE.
      NUMU      = 6
      UMU( 1 )  = -0.981986
      UMU( 2 )  = -0.538263
      UMU( 3 )  = -0.018014
      UMU( 4 )  =  0.018014
      UMU( 5 )  =  0.538263
      UMU( 6 )  =  0.981986
      NPHI      = 1
      PHI( 1 )  = 0.0
      IBCND     = 0
      FBEAM     = PI
      UMU0      = 0.080442
      PHI0      = 0.0
      FISOT     = 0.0
      LAMBER    = .TRUE.
      ALBEDO    = 0.0
      DELTAM    = .FALSE.
      PLANK     = .FALSE.
      ONLYFL    = .FALSE.

      ICAS = 0
      DO 20  IOD = 1, 2

         IF ( IOD.EQ.1 )  UTAU( 2 ) = 0.2
         IF ( IOD.EQ.2 )  UTAU( 2 ) = 5.0

         DTAUC( 1 ) = UTAU( 2 )

         DO 20  ISS = 1, 2

            IF( ISS.EQ.1 )  SSALB( 1 ) = 0.5
            IF( ISS.EQ.2 )  SSALB( 1 ) = 1.0
            ICAS = ICAS + 1

            WRITE( HEADER, '(3A,F5.2,A,F9.6,A,F4.2)')
     &             'Test Case No. 2', ABC(ICAS),
     &             ', Rayleigh Scattering, Ref. SW, Table 1:  tau =',
     &             UTAU(2), ', mu0 =', UMU0, ', ss-albedo =', SSALB(1)

            CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                    WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                    NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                    PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                    TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                    HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                    MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                    UU, U0U, ALBMED, TRNMED )

            AZMAVG = .FALSE.
            NPROB = 2
            IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                    FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                    TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                    TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                    MAXTAU, MAXMU, MAXAZ )
   20 CONTINUE

      ENDIF


      IF( DOPROB(3) )  THEN

c **********************************************************************
c ****  Test Problem 3:  Henyey-Greenstein Scattering               ****
c ****  (Compare To Ref. VH2, Table 35)                             ****
c **********************************************************************

      NSTR = 16
      NLYR = 1
      CALL  GETMOM( 3, 0.75, NSTR, PMOM )
      USRTAU    = .TRUE.
      NTAU      = 2
      UTAU( 1 ) = 0.0
      USRANG    = .TRUE.
      NUMU      = 6
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.5
      UMU( 3 )  = -0.1
      UMU( 4 )  =  0.1
      UMU( 5 )  =  0.5
      UMU( 6 )  =  1.0
      NPHI      = 1
      PHI( 1 )  = 0.0
      IBCND     = 0
      UMU0      = 0.1
      PHI0      = 0.0
      DELTAM    = .TRUE.
      LAMBER    = .TRUE.
      ONLYFL    = .FALSE.
      ALBEDO    = 0.0
      PLANK     = .FALSE.

      DO 30  ICAS = 1, 6

         IF ( ICAS.EQ.1 ) THEN

            UTAU( 2 )  = 1.0
            SSALB( 1 ) = 0.2
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.2 ) THEN

            UTAU( 2 )  = 1.0
            SSALB( 1 ) = 1.0
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.3 ) THEN

            UTAU( 2 )  = 1.0
            SSALB( 1 ) = 0.99
            FBEAM      = 0.0
            FISOT      = 1.0

         ELSE IF ( ICAS.EQ.4 ) THEN

            UTAU( 2 )  = 8.0
            SSALB( 1 ) = 0.2
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.5 ) THEN

            UTAU( 2 )  = 8.0
            SSALB( 1 ) = 1.0
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.6 ) THEN

            UTAU( 2 )  = 8.0
            SSALB( 1 ) = 0.99
            FBEAM      = 0.0
            FISOT      = 1.0

         END IF

         DTAUC( 1 ) = UTAU( 2 )

         WRITE( HEADER, '(3A,F9.5,A,F5.2)') 'Test Case No. 3',ABC(ICAS),
     &         ', Henyey-Greenstein Scattering, Ref. VH2, Table 35,'
     &         //' g = 0.75, b =', UTAU(2), ', a =', SSALB(1)

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         AZMAVG = .TRUE.
         NPROB = 3
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
30    CONTINUE

      ENDIF


      IF( DOPROB(4) )  THEN

c **********************************************************************
c ****  Test Problem 4:  Haze-L Scattering, Beam Source             ****
c ****  (Compare to Ref. GS, Tables 12-16)                          ****
c **********************************************************************

      NSTR = 32
      NLYR = 1
      CALL  GETMOM( 4, 0.0, NSTR, PMOM )
      DTAUC( 1 ) = 1.0
      USRTAU     = .TRUE.
      NTAU       = 3
      UTAU( 1 )  = 0.0
      UTAU( 2 )  = 0.5
      UTAU( 3 )  = 1.0
      USRANG     = .TRUE.
      NUMU       = 6
      UMU( 1 )   = -1.0
      UMU( 2 )   = -0.5
      UMU( 3 )   = -0.1
      UMU( 4 )   =  0.1
      UMU( 5 )   =  0.5
      UMU( 6 )   =  1.0
      IBCND      = 0
      FBEAM      = PI
      PHI0       = 0.0
      FISOT      = 0.0
      LAMBER     = .TRUE.
      ALBEDO     = 0.0
      DELTAM     = .TRUE.
      PLANK      = .FALSE.
      ONLYFL     = .FALSE.

      DO 40 ICAS = 1, 3

         WRITE( TITLE, '(3A)' ) 'Test Case No. 4', ABC(ICAS),
     &          ', Haze-L Scattering, Ref. GS, Table '
         LENTIT = INDEX( TITLE,BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            SSALB( 1 ) = 1.0
            NPHI       = 1
            PHI( 1 )   = 0.0
            UMU0       = 1.0
            HEADER = TITLE(1:LENTIT) // ' 12'

         ELSE IF ( ICAS.EQ.2 ) THEN

            SSALB( 1 ) = 0.9
            NPHI       = 1
            PHI( 1 )   = 0.0
            UMU0       = 1.0
            HEADER = TITLE(1:LENTIT) // ' 13'

         ELSE IF ( ICAS.EQ.3 ) THEN

            SSALB( 1 ) = 0.9
            NPHI       = 3
            PHI( 1 )   = 0.0
            PHI( 2 )   = 90.0
            PHI( 3 )   = 180.0
            UMU0       = 0.5
            HEADER = TITLE(1:LENTIT) // ' 14-16'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         AZMAVG = .FALSE.
         NPROB = 4
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
40    CONTINUE

      ENDIF


      IF( DOPROB(5) )  THEN

c **********************************************************************
c ****  Test Problem 5:  Cloud C.1 Scattering, Beam Source          ****
c ****  (Compare to Ref. GS, Tables 19-20)                          ****
c **********************************************************************

      NSTR = 48
      NLYR = 1
      CALL  GETMOM( 5, 0.0, NSTR, PMOM )
      DTAUC( 1 ) = 64.0
      USRTAU     = .TRUE.
      NTAU       = 3
      USRANG     = .TRUE.
      NUMU       = 6
      UMU( 1 )   = -1.0
      UMU( 2 )   = -0.5
      UMU( 3 )   = -0.1
      UMU( 4 )   =  0.1
      UMU( 5 )   =  0.5
      UMU( 6 )   =  1.0
      NPHI       = 1
      PHI( 1 )   = 0.0
      IBCND      = 0
      FBEAM      = PI
      UMU0       = 1.0
      PHI0       = 0.0
      FISOT      = 0.0
      LAMBER     = .TRUE.
      ALBEDO     = 0.0
      DELTAM     = .TRUE.
      PLANK      = .FALSE.
      ONLYFL     = .FALSE.

      DO 50 ICAS = 1, 2

         WRITE( TITLE, '(3A)' ) 'Test Case No. 5', ABC(ICAS),
     &          ', Cloud C.1 Scattering, Ref. GS, Table '
         LENTIT = INDEX( TITLE,BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 32.0
            UTAU( 3 )  = 64.0
            SSALB( 1 ) = 1.0
            HEADER = TITLE(1:LENTIT) // ' 19'

         END IF

         IF ( ICAS.EQ.2 ) THEN

            UTAU( 1 )  =  3.2
            UTAU( 2 )  = 12.8
            UTAU( 3 )  = 48.0
            SSALB( 1 ) = 0.9
            HEADER = TITLE(1:LENTIT) // ' 20'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         AZMAVG = .FALSE.
         NPROB = 5
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
50    CONTINUE

      ENDIF


      IF( DOPROB(6) )  THEN

c **********************************************************************
c ****  Test Problem 6:  No Scattering, Increasingly Complex Sources****
c **********************************************************************

      NSTR       = 16
      NLYR       = 1
      SSALB( 1 ) = 0.0
      WVNMLO     = 0.0
      WVNMHI     = 50000.
      USRTAU     = .TRUE.
      USRANG     = .TRUE.
      NUMU       = 4
      UMU( 1 )   = -1.0
      UMU( 2 )   = -0.1
      UMU( 3 )   =  0.1
      UMU( 4 )   =  1.0
      NPHI       = 1
      PHI( 1 )   = 90.0
      IBCND      = 0
      FBEAM      = 200.0
      UMU0       = 0.5
      PHI0       = 0.0
      FISOT      = 0.0
      TEMIS      = 1.0
      ONLYFL     = .FALSE.
      DELTAM     = .FALSE.

      DO 60  ICAS = 1, 8

         WRITE( TITLE, '(3A)' ) 'Test Case No. 6', ABC(ICAS),
     &          ': No Scattering; Source = Beam'
         LENTIT = INDEX( TITLE, BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            NTAU = 2
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 0.0

         ELSE IF ( ICAS.GT.1 ) THEN

            NTAU = 3
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 0.5
            UTAU( 3 ) = 1.0

         END IF

         IF ( ICAS.EQ.1 ) THEN
c                                    ** Transparent medium, beam source
            DTAUC( 1 ) = 0.0
            LAMBER = .TRUE.
            ALBEDO = 0.0
            PLANK = .FALSE.
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = 0'

         ELSE IF ( ICAS.EQ.2 ) THEN
c                                    ** Add some optical depth
            DTAUC( 1 ) = 1.0
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = 0'

         ELSE IF ( ICAS.EQ.3 ) THEN
c                                   ** Add some isotropic reflection
            LAMBER = .TRUE.
            ALBEDO = 0.50
            PLANK = .FALSE.
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo=0.5 Lambert'

         ELSE IF ( ICAS.EQ.4 ) THEN
c                                   ** Use non-isotropic reflection
            DTAUC( 1 ) = 1.0
            LAMBER = .FALSE.
            DO  59  K = 0, NSTR
               HL( K ) = 0.7**K
59          CONTINUE
            PLANK = .FALSE.
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = Non-Lambert'

         ELSE IF ( ICAS.EQ.5 ) THEN
c                                   ** Add some bottom-boundary emission
            DTAUC( 1 ) = 1.0
            TEMPER( 0 ) = 0.0
            TEMPER( 1 ) = 0.0
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 0.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &               ', Bottom Emission; Bott Alb = Non-Lambert'

         ELSE IF ( ICAS.EQ.6 ) THEN
c                                   ** Add some top-boundary diffuse
c                                      incidence (prescribed + emitted)
            DTAUC( 1 ) = 1.0
            TEMPER( 0 ) = 0.0
            TEMPER( 1 ) = 0.0
            FISOT  = 100.0 / PI
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 250.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &               ', Bottom+Top Emission; Bott Alb = Non-Lambert'

         ELSE IF ( ICAS.EQ.7 ) THEN
c                                   ** Add some internal emission
            DTAUC( 1 ) = 1.0
            TEMPER( 0 ) = 250.0
            TEMPER( 1 ) = 300.0
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 250.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &         ', Bottom+Top+Internal Emission; Bott Alb = Non-Lambert'

         ELSE IF ( ICAS.EQ.8 ) THEN
c                                   ** Increase the optical depth
            DTAUC( 1 ) = 10.0
            TEMPER( 0 ) = 250.0
            TEMPER( 1 ) = 300.0
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 1.0
            UTAU( 3 ) = 10.0
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 250.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &         ', Bottom+Top+Internal Emission; Bott Alb = Non-Lambert'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         AZMAVG = .FALSE.
         NPROB = 6
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
60    CONTINUE

      ENDIF


      IF( DOPROB(7) )  THEN

c **********************************************************************
c ****  Test Problem 7:  Absorption + Scattering + All Possible     ****
c ****  Sources, Various Surface Reflectivities ( One Layer )       ****
c **********************************************************************

      NSTR = 12
      NLYR = 1
      CALL  GETMOM( 3, 0.8, NSTR, PMOM )
      DTAUC( 1 )  = 1.0
      SSALB( 1 )  = 0.5
      TEMPER( 0 ) = 300.0
      TEMPER( 1 ) = 200.0
      WVNMLO      = 0.0
      WVNMHI      = 50000.
      USRTAU      = .TRUE.
      NTAU        = 3
      UTAU( 1 )   = 0.0
      UTAU( 2 )   = 0.5
      UTAU( 3 )   = 1.0
      USRANG      = .TRUE.
      NUMU        = 4
      UMU( 1 )    = -1.0
      UMU( 2 )    = -0.1
      UMU( 3 )    =  0.1
      UMU( 4 )    =  1.0
      NPHI        = 2
      PHI( 1 )    = 0.0
      PHI( 2 )    = 90.0
      IBCND       = 0
      FBEAM       = 200.0
      UMU0        = 0.5
      PHI0        = 0.0
      FISOT       = 100.0
      BTEMP       = 320.0
      TTEMP       = 100.0
      TEMIS       = 1.0
      DELTAM      = .TRUE.
      PLANK       = .TRUE.
      ONLYFL      = .FALSE.

      DO 70  ICAS = 1, 3

         WRITE( TITLE, '(3A)' ) 'Test Case No. 7', ABC(ICAS),
     &       ': Absorption + Henyey-Greenstein Scattering, All Sources'

         LENTIT = INDEX( TITLE, BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            LAMBER = .TRUE.
            ALBEDO = 0.0
            HEADER = TITLE(1:LENTIT) // ', Bottom Albedo = 0'

         ELSE IF ( ICAS.EQ.2 ) THEN

            LAMBER = .TRUE.
            ALBEDO = 1.0
            HEADER = TITLE(1:LENTIT) // ', Bottom Albedo = 1'

         ELSE IF ( ICAS.EQ.3 ) THEN

            LAMBER = .FALSE.
            DO 69 K = 0, NSTR
               HL( K ) = 0.7**K
69          CONTINUE
            HEADER = TITLE(1:LENTIT) // ', Bottom Albedo = BDR Function'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         AZMAVG = .FALSE.
         NPROB = 7
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
70    CONTINUE

      ENDIF


      IF( DOPROB(8) )  THEN

c **********************************************************************
c ****  Test Problem 8:  Absorbing/Isotropic-Scattering Medium      ****
c ****  With Two Computational Layers                               ****
c **** (Compare Fluxes To Ref. OS, Table 1)                         ****
c **********************************************************************

      NSTR = 8
      NLYR = 2
      CALL  GETMOM( 1, 0.0, NSTR, PMOM(0,1) )
      CALL  GETMOM( 1, 0.0, NSTR, PMOM(0,2) )
      USRTAU   = .TRUE.
      USRANG   = .TRUE.
      NUMU     = 4
      UMU( 1 ) = -1.0
      UMU( 2 ) = -0.2
      UMU( 3 ) =  0.2
      UMU( 4 ) =  1.0
      NPHI     = 1
      PHI( 1 ) = 60.0
      IBCND    = 0
      FBEAM    = 0.0
      FISOT    = 1.0 / PI
      LAMBER   = .TRUE.
      ALBEDO   = 0.0
      PLANK    = .FALSE.
      DELTAM   = .FALSE.
      ONLYFL   = .FALSE.

      DO 80  ICAS = 1, 3

         IF ( ICAS.EQ.1 ) THEN

            DTAUC( 1 ) = 0.25
            DTAUC( 2 ) = 0.25
            SSALB( 1 ) = 0.5
            SSALB( 2 ) = 0.3
            NTAU       = 3
            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 0.25
            UTAU( 3 )  = 0.5
            HEADER = 'Test Case No. 8A:  Ref. OS, Table 1,'
     &               // ' Line 4 (Two Inhomogeneous Layers)'

         ELSE IF ( ICAS.EQ.2 ) THEN

            DTAUC( 1 ) = 0.25
            DTAUC( 2 ) = 0.25
            SSALB( 1 ) = 0.8
            SSALB( 2 ) = 0.95
            NTAU       = 3
            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 0.25
            UTAU( 3 )  = 0.5
            HEADER = 'Test Case No. 8b:  Ref. OS, Table 1,'
     &               // ' Line 1 (Two Inhomogeneous Layers)'

         ELSE IF ( ICAS.EQ.3 ) THEN

            DTAUC( 1 ) = 1.0
            DTAUC( 2 ) = 2.0
            SSALB( 1 ) = 0.8
            SSALB( 2 ) = 0.95
            NTAU       = 3
            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 1.0
            UTAU( 3 )  = 3.0
            HEADER = 'Test Case No. 8c:  Ref. OS, Table 1,'
     &               // ' Line 13 (Two Inhomogeneous Layers)'
         ENDIF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         NPROB = 8
         AZMAVG = .FALSE.
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
  80  CONTINUE

      ENDIF


      IF( DOPROB(9) )  THEN

c **********************************************************************
c ****  Test Problem 9:  General Emitting/Absorbing/Scattering      ****
c ****  Medium with Every Computational Layer Different.            ****
c **** (Compare 9a,b Fluxes to Ref. DGIS, Tables VI-VII, beta = 0)  ****
c **********************************************************************

      NSTR = 8
      NLYR = 6
      DO 86  LC = 1, NLYR
         DTAUC( LC ) = LC
         SSALB( LC ) = 0.6 + LC*0.05
 86   CONTINUE
      USRTAU    = .TRUE.
      NTAU      = 5
      UTAU( 1 ) = 0.0
      UTAU( 2 ) = 1.05
      UTAU( 3 ) = 2.1
      UTAU( 4 ) = 6.0
      UTAU( 5 ) = 21.0
      USRANG    = .TRUE.
      NUMU      = 4
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.2
      UMU( 3 )  =  0.2
      UMU( 4 )  =  1.0
      NPHI      = 1
      PHI( 1 )  = 60.0
      IBCND     = 0
      FBEAM     = 0.0
      FISOT     = 1.0 / PI
      LAMBER    = .TRUE.
      DELTAM    = .TRUE.
      ONLYFL    = .FALSE.

      DO 90  ICAS = 1, 3

         IF ( ICAS.EQ.1 ) THEN

            DO 87  LC = 1, NLYR
               CALL  GETMOM( 1, 0.0, NSTR, PMOM(0,LC) )
 87         CONTINUE
            ALBEDO = 0.0
            PLANK = .FALSE.
            HEADER = 'Test Case No. 9a:  Ref. DGIS, Tables VI-VII,'
     &               // ' beta=l=0 (multiple inhomogeneous layers)'

         ELSE IF ( ICAS.EQ.2 ) THEN

            PMOM(0,1) = 1.0
            PMOM(1,1) = 2.00916/3.
            PMOM(2,1) = 1.56339/5.
            PMOM(3,1) = 0.67407/7.
            PMOM(4,1) = 0.22215/9.
            PMOM(5,1) = 0.04725/11.
            PMOM(6,1) = 0.00671/13.
            PMOM(7,1) = 0.00068/15.
            PMOM(8,1) = 0.00005/17.
            DO 88  LC = 2, NLYR
               DO 88  K = 0, 8
                  PMOM(K,LC) = PMOM(K,1)
 88         CONTINUE
            HEADER = 'Test Case No. 9b:  Ref. DGIS, Tables VI-VII,'
     &               // ' beta=0,l=8 (multiple inhomogeneous layers)'

         ELSE IF ( ICAS.EQ.3 ) THEN

            TEMPER( 0 ) = 600.0
            DO 89  LC = 1, NLYR
               CALL  GETMOM( 3, FLOAT(LC)/7.0, NSTR, PMOM(0,LC) )
               TEMPER( LC ) = 600.0 + LC*10.0
 89         CONTINUE
            NPHI = 3
            PHI( 1 ) = 60.0
            PHI( 2 ) = 120.0
            PHI( 3 ) = 180.0
            FBEAM    = PI
            UMU0     = 0.5
            PHI0     = 0.0
            FISOT    = 1.0
            ALBEDO   =  0.5
            PLANK    = .TRUE.
            WVNMLO   =  999.0
            WVNMHI   = 1000.0
            BTEMP    = 700.0
            TTEMP    = 550.0
            TEMIS    = 1.0
            HEADER = 'Test Case No. 9c:  Generalization of 9A '//
     &               'to include all possible complexity'
         ENDIF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         NPROB = 9
         AZMAVG = .FALSE.
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                 FLUP, DFDT, U0U, UU, TSTFIR(1,ICAS,NPROB),
     &                 TSTFDN(1,ICAS,NPROB), TSTFUP(1,ICAS,NPROB),
     &                 TSTDFD(1,ICAS,NPROB), TSTUU(1,1,1,ICAS,NPROB),
     &                 MAXTAU, MAXMU, MAXAZ )
  90  CONTINUE

      ENDIF


      IF( DOPROB(10) )  THEN

c **********************************************************************
c ****  Test Problem 10: Compare USRANG = True With USRANG = False  ****
c ****  take Problem 9c (our most general case) but only 4 Streams  ****
c **********************************************************************

      NSTR = 4
      NLYR = 6
      TEMPER( 0 ) = 600.0
      DO 97  LC = 1, NLYR
         DTAUC( LC ) = LC
         SSALB( LC ) = 0.6 + LC*0.05
         CALL  GETMOM( 3, FLOAT(LC)/(NLYR+1), NSTR, PMOM(0,LC) )
         TEMPER( LC ) = 600.0 + LC*10.0
  97  CONTINUE
      USRTAU    = .TRUE.
      NTAU      = 3
      UTAU( 1 ) = 0.0
      UTAU( 2 ) = 2.1
      UTAU( 3 ) = 21.0
      NPHI      = 2
      PHI( 1 )  = 60.0
      PHI( 2 )  = 120.0
      IBCND     = 0
      FBEAM     = PI
      UMU0      = 0.5
      PHI0      = 0.0
      FISOT     = 1.0
      LAMBER    = .TRUE.
      DELTAM    = .TRUE.
      ALBEDO    =  0.5
      PLANK     = .TRUE.
      WVNMLO    =  999.0
      WVNMHI    = 1000.0
      BTEMP     = 700.0
      TTEMP     = 550.0
      TEMIS     = 1.0
      ONLYFL    = .FALSE.

      DO 100  ICAS = 1, 2

         IF ( ICAS.EQ.1 ) THEN

            USRANG = .TRUE.
            NUMU      = 4
            UMU( 1 )  = - 0.788675129
            UMU( 2 )  = - 0.211324871
            UMU( 3 )  =   0.211324871
            UMU( 4 )  =   0.788675129
            PRNT( 2 ) = .TRUE.
            PRNT( 5 ) = .TRUE.
            HEADER = 'Test Case No. 10a:  like 9c, USRANG = True'

         ELSE IF ( ICAS.EQ.2 ) THEN

            USRANG    = .FALSE.
            NUMU      = 0
            PRNT( 2 ) = .FALSE.
            PRNT( 5 ) = .FALSE.
            HEADER = 'Test Case No. 10b:  like 9C, USRANG = False'

         ENDIF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         IF ( ICAS.EQ.1 ) THEN
c                               ** Save results to compare to case 2
            DO 98  LU = 1, NTAU
               CMPFIR( LU ) = RFLDIR( LU )
               CMPFDN( LU ) = RFLDN( LU )
               CMPFUP( LU ) = FLUP( LU )
               CMPDFD( LU ) = DFDT( LU )
               DO 98  IU = 1, NUMU
                  DO 98  J = 1, NPHI
                     CMPUU( LU, IU, J ) = UU( IU, LU, J )
   98       CONTINUE

         ELSE IF ( ICAS.EQ.2 ) THEN

            AZMAVG = .FALSE.
            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                    FLUP, DFDT, U0U, UU, CMPFIR, CMPFDN, CMPFUP,
     &                    CMPDFD, CMPUU, MXTAU, MXMU, MXPHI )
         END IF

 100  CONTINUE

      ENDIF


      IF( DOPROB(11) )  THEN

c **********************************************************************
c ****  Test Problem 11: Single-Layer vs. Multiple Layers           ****
c ****  11a: Results at user levels for one computational layer     ****
c ****  11b: Single layer of 11a subdivided into multiple           ****
c ****       computational layers at the 11a user levels            ****
c **********************************************************************

      NSTR      = 16
      USRANG    = .TRUE.
      NUMU      = 4
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.1
      UMU( 3 )  =  0.1
      UMU( 4 )  =  1.0
      NPHI      = 2
      PHI( 1 )  = 0.0
      PHI( 2 )  = 90.0
      IBCND     = 0
      FBEAM     = 1.0
      UMU0      = 0.5
      PHI0      = 0.0
      FISOT     = 0.5 / PI
      LAMBER    = .TRUE.
      ALBEDO    = 0.5
      DELTAM    = .FALSE.
      PLANK     = .FALSE.
      ONLYFL    = .FALSE.

      DO 110  ICAS = 1, 2

         IF ( ICAS.EQ.1 ) THEN

            NLYR = 1
            DTAUC( 1 ) = 1.0
            SSALB( 1 ) = 0.9
            CALL  GETMOM( 1, 0.0, NSTR, PMOM )
            USRTAU    = .TRUE.
            NTAU      = 4
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 0.05
            UTAU( 3 ) = 0.5
            UTAU( 4 ) = 1.0
            PRNT( 2 ) = .TRUE.
            PRNT( 5 ) = .TRUE.
            HEADER = 'Test Case No. 11a: One Isotropic-Scattering Layer'

         ELSE IF ( ICAS.EQ.2 ) THEN

            NLYR = NTAU - 1
            DO 107 LC = 1, NLYR
               DTAUC( LC ) = UTAU(LC+1) - UTAU(LC)
               SSALB( LC ) = 0.9
               CALL  GETMOM( 1, 0.0, NSTR, PMOM(0,LC) )
107         CONTINUE
            USRTAU    = .FALSE.
            PRNT( 2 ) = .FALSE.
            PRNT( 5 ) = .FALSE.
            HEADER = 'Test Case No. 11b: Same as 11a but treated as' //
     &               ' multiple layers'
         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         IF ( ICAS.EQ.1 ) THEN
c                               ** Save results to compare to case 2
            DO 108  LU = 1, NTAU
               CMPFIR( LU ) = RFLDIR( LU )
               CMPFDN( LU ) = RFLDN( LU )
               CMPFUP( LU ) = FLUP( LU )
               CMPDFD( LU ) = DFDT( LU )
               DO 108  IU = 1, NUMU
                  DO 108  J = 1, NPHI
                     CMPUU( LU, IU, J ) = UU( IU, LU, J )
  108       CONTINUE

         ELSE IF ( ICAS.EQ.2 ) THEN

            AZMAVG = .FALSE.
            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                    FLUP, DFDT, U0U, UU, CMPFIR, CMPFDN, CMPFUP,
     &                    CMPDFD, CMPUU, MXTAU, MXMU, MXPHI )
         END IF

110   CONTINUE

      ENDIF


      IF( DOPROB(12) )  THEN

c **********************************************************************
c ****  Test Problem 12: Test Absorption-Optical-Depth Shortcut     ****
c ****  compares cases where the DISORT shortcut for absorption     ****
c ****  optical depth .GT. 10 is not used (12a), then is used (12b) ****
c ****  (this shortcut is only employed when  PLANK = False.)       ****
c **********************************************************************

      NSTR    = 20
      USRANG  = .TRUE.
      NUMU     = 4
      UMU( 1 ) = -1.0
      UMU( 2 ) = -0.1
      UMU( 3 ) =  0.1
      UMU( 4 ) =  1.0
      NPHI     = 1
      PHI( 1 ) = 0.0
      IBCND    = 0
      FBEAM    = 1.0
      UMU0     = 1.0
      PHI0     = 0.0
      FISOT    = 0.0
      LAMBER   = .TRUE.
      ALBEDO   = 1.0
      DELTAM   = .TRUE.
      PLANK    = .FALSE.
      ONLYFL   = .FALSE.

      DO 120  ICAS = 1, 2

         IF ( ICAS.EQ.1 ) THEN

            NLYR = 1
            DTAUC( 1 ) = 20.1
            SSALB( 1 ) = 0.5
            CALL  GETMOM( 3, 0.9, NSTR, PMOM )
            USRTAU    = .TRUE.
            NTAU      = 4
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 10.0
            UTAU( 3 ) = 19.9
            UTAU( 4 ) = 20.1
            PRNT( 2 ) = .TRUE.
            PRNT( 5 ) = .TRUE.
            HEADER = 'Test Case No. 12a:  Overhead Beam Striking '//
     &               'Absorbing/Scattering Medium'

         ELSE IF ( ICAS.EQ.2 ) THEN

            NLYR = NTAU - 1
            DO 117 LC = 1, NLYR
               DTAUC( LC ) = UTAU(LC+1) - UTAU(LC)
               SSALB( LC ) = 0.5
               CALL  GETMOM( 3, 0.9, NSTR, PMOM(0,LC) )
117         CONTINUE
            USRTAU    = .FALSE.
            PRNT( 2 ) = .FALSE.
            PRNT( 5 ) = .FALSE.
            HEADER = 'Test Case No. 12b: Same as 12a but uses shortcut'
     &               // ' for absorption optical depth .GT. 10'
         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )

         IF ( ICAS.EQ.1 ) THEN
c                               ** Save results to compare to case 2
            DO 118  LU = 1, NTAU
               CMPFIR( LU ) = RFLDIR( LU )
               CMPFDN( LU ) = RFLDN( LU )
               CMPFUP( LU ) = FLUP( LU )
               CMPDFD( LU ) = DFDT( LU )
               DO 118  IU = 1, NUMU
                  DO 118  J = 1, NPHI
                     CMPUU( LU, IU, J ) = UU( IU, LU, J )
  118       CONTINUE

         ELSE IF ( ICAS.EQ.2 ) THEN

            AZMAVG = .FALSE.
            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN,
     &                    FLUP, DFDT, U0U, UU, CMPFIR, CMPFDN, CMPFUP,
     &                    CMPDFD, CMPUU, MXTAU, MXMU, MXPHI )
         END IF

120   CONTINUE

      ENDIF


      IF( DOPROB(13) )  THEN

c **********************************************************************
c ****  Test Problem 13: Test shortcut for flux albedo, transmission ***
c **** ( shortcut gives flux albedo, transmission of entire medium  ****
c ****   as a function of sun angle )                               ****
c ****  13a,c = Shortcut;  13b,d = Brute Force Method               ****
c **********************************************************************

      NSTR   = 16
      NPHI   = 0
      PHI0   = 0.0
      ALBEDO = 0.5
      DELTAM = .TRUE.

      DO 130  ICAS = 1, 4

         IF ( ICAS.EQ.1 ) THEN

            IBCND      = 1
            NLYR       = 1
            DTAUC( 1 ) = 1.0
            SSALB( 1 ) = 0.99
            CALL  GETMOM( 3, 0.8, NSTR, PMOM )
            PRNT( 6 )  = .TRUE.
            PRNT( 2 )  = .FALSE.
            USRANG     = .TRUE.
            NUMU       = 1
            UMU( 1 )   =  0.5
            HEADER = 'Test Case No. 13a:  Albedo and Transmissivity'//
     &               ' from Shortcut, Single Layer'

         ELSE IF ( ICAS.EQ.2 ) THEN

            IBCND     = 0
            USRTAU    = .TRUE.
            NTAU      = 2
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 1.0
            UMU0      = 0.5
            FBEAM     = 1.0 / UMU0
            FISOT     = 0.0
            LAMBER    = .TRUE.
            PLANK     = .FALSE.
            ONLYFL    = .TRUE.
            PRNT( 6 ) = .FALSE.
            PRNT( 2 ) = .TRUE.
            HEADER = 'Test Case No. 13b:  Albedo and Transmissivity'//
     &               ' by Regular Method, Single Layer'

         ELSE IF ( ICAS.EQ.3 ) THEN

            IBCND     = 1
            PRNT( 6 ) = .TRUE.
            PRNT( 2 ) = .FALSE.
            NLYR      = 2
            DO 125  LC = 1, NLYR
               DTAUC( LC ) = 1.0 / NLYR
               CALL  GETMOM( 3, 0.8, NSTR, PMOM(0,LC) )
  125       CONTINUE
            SSALB( 1 ) = 0.99
            SSALB( 2 ) = 0.50
            HEADER = 'Test Case No. 13c:  Albedo and Transmissivity'//
     &               ' from Shortcut, Multiple Layer'

         ELSE IF ( ICAS.EQ.4 ) THEN

            IBCND = 0
            PRNT( 6 ) = .FALSE.
            PRNT( 2 ) = .TRUE.
            HEADER = 'Test Case No. 13d:  Albedo and Transmissivity'//
     &               ' by Regular Method, Multiple Layer'
         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                 WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG,
     &                 NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0,
     &                 PHI0, FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP,
     &                 TEMIS, DELTAM, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU,
     &                 MAXPHI, RFLDIR, RFLDN, FLUP, DFDT, UAVG,
     &                 UU, U0U, ALBMED, TRNMED )
130   CONTINUE

      ENDIF

      STOP
      END

      SUBROUTINE  GETMOM( IPHAS, GG, NMOM, PMOM )

c        Calculate phase function Legendre expansion coefficients
c        in various special cases


c       INPUT: IPHAS   Phase function options
c                      1 : Isotropic
c                      2 : Rayleigh
c                      3 : Henyey-Greenstein with asymmetry factor GG
c                      4 : Haze L as specified by Garcia/Siewert
c                      5 : Cloud C.1 as specified by Garcia/Siewert

c              GG      Asymmetry factor for Henyey-Greenstein case

c              NMOM    Index of highest Legendre coefficient needed
c                        ( = number of streams 'NSTR'  chosen
c                         for the discrete ordinate method)

c      OUTPUT: PMOM(K)  Legendre expansion coefficients (K=0 to NMOM)
c                         (be sure to dimension '0:maxval' in calling
c                          program)

c      Reference:  Garcia, R. and C. Siewert, 1985: Benchmark Results
c                     in Radiative Transfer, Transp. Theory and Stat.
c                     Physics 14, 437-484, Tables 10 And 17
c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   IPHAS, NMOM
      REAL      GG
c     ..
c     .. Array Arguments ..

      REAL      PMOM( 0:NMOM )
c     ..
c     .. Local Scalars ..

      INTEGER   K
c     ..
c     .. Local Arrays ..

      REAL      CLDMOM( 299 ), HAZELM( 82 )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..

      DATA HAZELM /  2.41260, 3.23047, 3.37296, 3.23150, 2.89350,
     A               2.49594, 2.11361, 1.74812, 1.44692, 1.17714,
     B               0.96643, 0.78237, 0.64114, 0.51966, 0.42563,
     C               0.34688, 0.28351, 0.23317, 0.18963, 0.15788,
     D               0.12739, 0.10762, 0.08597, 0.07381, 0.05828,
     E               0.05089, 0.03971, 0.03524, 0.02720, 0.02451,
     F               0.01874, 0.01711, 0.01298, 0.01198, 0.00904,
     G               0.00841, 0.00634, 0.00592, 0.00446, 0.00418,
     H               0.00316, 0.00296, 0.00225, 0.00210, 0.00160,
     I               0.00150, 0.00115, 0.00107, 0.00082, 0.00077,
     J               0.00059, 0.00055, 0.00043, 0.00040, 0.00031,
     K               0.00029, 0.00023, 0.00021, 0.00017, 0.00015,
     L               0.00012, 0.00011, 0.00009, 0.00008, 0.00006,
     M               0.00006, 0.00005, 0.00004, 0.00004, 0.00003,
     N               0.00003, 3*0.00002, 8*0.00001 /

      DATA  ( CLDMOM(K), K = 1, 159 ) /
     A  2.544,  3.883,  4.568,  5.235,  5.887,  6.457,  7.177,  7.859,
     B  8.494,  9.286,  9.856, 10.615, 11.229, 11.851, 12.503, 13.058,
     C 13.626, 14.209, 14.660, 15.231, 15.641, 16.126, 16.539, 16.934,
     D 17.325, 17.673, 17.999, 18.329, 18.588, 18.885, 19.103, 19.345,
     E 19.537, 19.721, 19.884, 20.024, 20.145, 20.251, 20.330, 20.401,
     F 20.444, 20.477, 20.489, 20.483, 20.467, 20.427, 20.382, 20.310,
     G 20.236, 20.136, 20.036, 19.909, 19.785, 19.632, 19.486, 19.311,
     H 19.145, 18.949, 18.764, 18.551, 18.348, 18.119, 17.901, 17.659,
     I 17.428, 17.174, 16.931, 16.668, 16.415, 16.144, 15.883, 15.606,
     J 15.338, 15.058, 14.784, 14.501, 14.225, 13.941, 13.662, 13.378,
     K 13.098, 12.816, 12.536, 12.257, 11.978, 11.703, 11.427, 11.156,
     L 10.884, 10.618, 10.350, 10.090,  9.827,  9.574,  9.318,  9.072,
     M  8.822, 8.584, 8.340, 8.110, 7.874, 7.652, 7.424, 7.211, 6.990,
     N  6.785, 6.573, 6.377, 6.173, 5.986, 5.790, 5.612, 5.424, 5.255,
     O  5.075, 4.915, 4.744, 4.592, 4.429, 4.285, 4.130, 3.994, 3.847,
     P  3.719, 3.580, 3.459, 3.327, 3.214, 3.090, 2.983, 2.866, 2.766,
     Q  2.656, 2.562, 2.459, 2.372, 2.274, 2.193, 2.102, 2.025, 1.940,
     R  1.869, 1.790, 1.723, 1.649, 1.588, 1.518, 1.461, 1.397, 1.344,
     S  1.284, 1.235, 1.179, 1.134, 1.082, 1.040, 0.992, 0.954, 0.909 /
      DATA  ( CLDMOM(K), K = 160, 299 ) /
     T  0.873, 0.832, 0.799, 0.762, 0.731, 0.696, 0.668, 0.636, 0.610,
     U  0.581, 0.557, 0.530, 0.508, 0.483, 0.463, 0.440, 0.422, 0.401,
     V  0.384, 0.364, 0.349, 0.331, 0.317, 0.301, 0.288, 0.273, 0.262,
     W  0.248, 0.238, 0.225, 0.215, 0.204, 0.195, 0.185, 0.177, 0.167,
     X  0.160, 0.151, 0.145, 0.137, 0.131, 0.124, 0.118, 0.112, 0.107,
     Y  0.101, 0.097, 0.091, 0.087, 0.082, 0.079, 0.074, 0.071, 0.067,
     Z  0.064, 0.060, 0.057, 0.054, 0.052, 0.049, 0.047, 0.044, 0.042,
     A  0.039, 0.038, 0.035, 0.034, 0.032, 0.030, 0.029, 0.027, 0.026,
     B  0.024, 0.023, 0.022, 0.021, 0.020, 0.018, 0.018, 0.017, 0.016,
     C  0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.010, 0.009,
     D  0.009, 3*0.008, 2*0.007, 3*0.006, 4*0.005, 4*0.004, 6*0.003,
     E  9*0.002, 18*0.001 /


      IF ( IPHAS.LT.1 .OR. IPHAS.GT.5 )
     &     CALL ERRMSG( 'GETMOM--bad input variable IPHAS',.TRUE.)

      IF ( IPHAS.EQ.3 .AND. (GG.LE.-1.0 .OR. GG.GE.1.0) )
     &     CALL ERRMSG( 'GETMOM--bad input variable GG',.TRUE.)

      IF ( NMOM.LT.2 )
     &     CALL ERRMSG( 'GETMOM--bad input variable NMOM',.TRUE.)


      PMOM(0) = 1.0
      DO  10  K = 1, NMOM
         PMOM(K) = 0.0
   10 CONTINUE


      IF ( IPHAS.EQ.2 )  THEN
c                                       ** Rayleigh phase function
         PMOM(2) = 0.1

      ELSE IF ( IPHAS.EQ.3 ) THEN
c                                       ** Henyey-Greenstein phase fcn
         DO  20  K = 1, NMOM
            PMOM(K) = GG**K
   20    CONTINUE

      ELSE IF ( IPHAS.EQ.4 ) THEN
c                                        ** Haze-L phase function
         DO  30  K = 1, MIN(82,NMOM)
            PMOM(K) = HAZELM(K) / ( 2*K+1 )
   30    CONTINUE

      ELSE IF ( IPHAS.EQ.5 ) THEN
c                                        ** Cloud C.1 phase function
         DO  40  K = 1, MIN(298,NMOM)
            PMOM(K) = CLDMOM(K) / ( 2*K+1 )
40       CONTINUE

      END IF

      END

      SUBROUTINE  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, AZMAVG, RFLDIR, RFLDN, FLUP,
     &                    DFDT, U0U, UU, TSTFIR, TSTFDN, TSTFUP,
     &                    TSTDFD, TSTUU, MAXTAU, MAXMU, MAXAZ )

c        Print DISORT results and, directly beneath them, their
c        ratios to the correct answers;  print number of non-unit
c        ratios that occur but try to count just the cases where 
c        there is a real disagreement and not those where flux or
c        intensity are down at their noise level (defined as 10^(-6)
c        times their maximum value).  d(flux)/d(tau) is treated the
c        same as fluxes in this noise estimation even though it
c        is a different type of quantity (although with flux units).

c     INPUT :   TSTFIR  correct direct flux
c               TSTFDN  correct diffuse down flux
c               TSTFUP  correct diffuse up flux
c               TSTDFD  correct d(flux)/d(optical depth)
c               TSTUU   correct intensity or azim-avg intensity
c               AZMAVG  TRUE, check azimuthally-averaged intensity only
c                       FALSE, check full intensity
c               (remaining input = DISORT I/O variables)

c --------------------------------------------------------------------

c     .. Parameters ..

      INTEGER   MAXRAT
      PARAMETER ( MAXRAT = 100 )
c     ..
c     .. Scalar Arguments ..

      LOGICAL   AZMAVG, ONLYFL
      INTEGER   MAXAZ, MAXMU, MAXTAU, MAXULV, MAXUMU, NPHI, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      DFDT( * ), FLUP( * ), PHI( * ), RFLDIR( * ), RFLDN( * ),
     &          TSTDFD( * ), TSTFDN( * ), TSTFIR( * ), TSTFUP( * ),
     &          TSTUU( MAXTAU, MAXMU, MAXAZ ), U0U( MAXUMU, * ),
     &          UMU( * ), UTAU( * ), UU( MAXUMU, MAXULV, * )
c     ..
c     .. Local Scalars ..

      INTEGER  IU, J, LU, NUMBAD
      REAL     FLXMAX, FNOISE, RAT, RAT1, RAT2, RAT3, RAT4, UMAX, UNOISE
c     ..
c     .. Local Arrays ..

      REAL      RATV( MAXRAT )
c     ..
c     .. External Functions ..

      REAL      RATIO
      EXTERNAL  RATIO
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Statement Functions ..

      LOGICAL   BADRAT
c     ..
c     .. Statement Function definitions ..

      BADRAT( RAT ) = (RAT.LT.0.99) .OR. (RAT.GT.1.01)
c     ..


      IF ( NTAU.GT.MAXTAU .OR. NUMU.GT.MAXMU .OR. NPHI.GT.MAXAZ )  CALL
     &   ERRMSG( 'PRTFIN--out of bounds in comparator arrays', .TRUE.)

      FLXMAX = 0.0
      DO 5  LU = 1, NTAU
         FLXMAX = MAX( FLXMAX, TSTFIR(LU), TSTFDN(LU), TSTFUP(LU) )
 5    CONTINUE
      FNOISE = 1.E-6 * FLXMAX
      IF( FLXMAX.LE.0.0 ) 
     &    CALL ERRMSG( 'PRTFIN--all fluxes zero or negative', .FALSE.)
      IF( FNOISE.LE.0.0 ) 
     &    CALL ERRMSG( 'PRTFIN--all fluxes near underflowing', .FALSE.)

      NUMBAD = 0

      WRITE(*,'(//,A,/,A,/,A)')
     &  '                  <-------------- FLUXES -------------->',
     &  '    Optical       Downward       Downward         Upward'//
     &  '    d(Net Flux)',
     &  '      Depth         Direct        Diffuse        Diffuse'//
     &  '    / d(Op Dep)'

      DO 10  LU = 1, NTAU

         WRITE( *,'(0P,F11.4,1P,4E15.4)')  UTAU(LU), RFLDIR(LU), 
     &          RFLDN(LU), FLUP(LU), DFDT(LU)

         RAT1 = RATIO( RFLDIR(LU), TSTFIR(LU) )
         RAT2 = RATIO( RFLDN(LU),  TSTFDN(LU) )
         RAT3 = RATIO(  FLUP(LU),  TSTFUP(LU) )
         RAT4 = RATIO(  DFDT(LU),  TSTDFD(LU) )
         WRITE( *,'(11X,4( ''    ('',F9.4,'')''))')  
     &          RAT1, RAT2, RAT3, RAT4

         IF( BADRAT(RAT1) .AND. ABS(RFLDIR(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT2) .AND. ABS(RFLDN(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT3) .AND. ABS(FLUP(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT4) .AND. ABS(DFDT(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1

10    CONTINUE

      IF ( ONLYFL )  GO TO 100

      IF ( NUMU.GT.MAXRAT .OR. NPHI.GT.MAXRAT )
     &     CALL ERRMSG( 'PRTFIN--increase parameter MAXRAT', .TRUE.)


      IF ( AZMAVG ) THEN
c                                       ** Print az-avg intensities

         IF ( NUMU.GT.8 ) CALL ERRMSG
     &      ( 'PRTFIN--az-avg-intensity FORMATs inadequate',.FALSE.)

         UMAX = 0.0
         DO 16  LU = 1, NTAU
            DO 15  IU = 1, NUMU
               UMAX = MAX( UMAX, TSTUU(LU,IU,1) )
 15         CONTINUE
 16      CONTINUE
         UNOISE = 1.E-6 * UMAX
         IF( UMAX.LE.0.0 )  CALL ERRMSG
     &     ( 'PRTFIN--all az-avg intensities zero or negative',.FALSE.)
         IF( UNOISE.LE.0.0 ) CALL ERRMSG
     &     ( 'PRTFIN--all az-avg intensities near underflowing',.FALSE.)

         WRITE (*,'(//,A,//,A,/,A,8F14.5)')
     &  ' ********* AZIMUTHALLY AVERAGED INTENSITIES  *********',
     &  '   Optical   Polar Angle Cosines',
     &  '     Depth', ( UMU(IU), IU = 1, NUMU )

         DO 30  LU = 1, NTAU

            WRITE( *,'(0P,F10.3,1P,8E14.4)') UTAU(LU), 
     &             ( U0U(IU,LU), IU = 1, NUMU )

            DO 20  IU = 1, NUMU
               RATV(IU) = RATIO( U0U(IU,LU), TSTUU(LU,IU,1) )
               IF( BADRAT(RATV(IU)) .AND. ABS(U0U(IU,LU)).GT.UNOISE )
     &               NUMBAD = NUMBAD + 1
   20       CONTINUE

            WRITE( *,'(10X, 8(:,''   ('',F9.4,'')''))')  
     &           ( RATV(IU), IU = 1, NUMU )

   30    CONTINUE


      ELSE
c                                       ** Print intensities

         IF ( NPHI.GT.8 ) CALL ERRMSG
     &      ( 'PRTFIN--intensity FORMATs inadequate',.FALSE.)

         UMAX = 0.0
         DO 36  LU = 1, NTAU
            DO 35  IU = 1, NUMU
               DO 34  J = 1, NPHI
                  UMAX = MAX( UMAX, TSTUU(LU,IU,J) )
 34            CONTINUE
 35         CONTINUE
 36      CONTINUE
         UNOISE = 1.E-6 * UMAX
         IF( UMAX.LE.0.0 )  CALL ERRMSG
     &     ( 'PRTFIN--all intensities zero or negative',.FALSE.)
         IF( UNOISE.LE.0.0 ) CALL ERRMSG
     &     ( 'PRTFIN--all intensities near underflowing',.FALSE.)

         WRITE( *,'(//,A,//,A,/,A,/,A,8(F10.1,4X))' )
     &            ' ********  I N T E N S I T I E S  *********',
     &            '             Polar   Azimuthal Angles (Degrees)',
     &            '   Optical   Angle',
     &            '     Depth  Cosine', ( PHI(J), J = 1, NPHI )

         DO 60  LU = 1, NTAU

            DO 50  IU = 1, NUMU

               IF( IU.EQ.1 ) WRITE( *,'(/,0P,F10.3,F8.3,1P,8E14.4)')
     &             UTAU(LU), UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

               IF( IU.GT.1 ) WRITE( *,'(10X,0P,F8.3, 1P,8E14.4)') 
     &                       UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

               DO 40  J = 1, NPHI
                  RATV(J) = RATIO( UU(IU,LU,J), TSTUU(LU,IU,J) )
                  IF( BADRAT(RATV(J)) .AND. ABS(UU(IU,LU,J)).GT.UNOISE )
     &                  NUMBAD = NUMBAD + 1
   40          CONTINUE

               WRITE( *,'(18X, 8(:,''   ('',F9.4,'')''))') 
     &              ( RATV(J), J = 1, NPHI )

   50       CONTINUE
   60    CONTINUE

      END IF

  100 CONTINUE
      IF( NUMBAD.GT.0 )  WRITE( *,300)  ' ====  ', NUMBAD,
     &    '  SERIOUSLY NON-UNIT RATIOS    ===='

      RETURN

300   FORMAT( //,1X,45('='),/,A,I4,A,/,1X,45('=') )
      END

      BLOCK DATA  CHEKDO

c       Correct answers to test problems ( as produced by DISORT
c       running in 14-digit precision on a Cray computer )

c     .. Parameters ..

      INTEGER   MXPROB, MXCASE, MAXTAU, MAXMU, MAXAZ
      PARAMETER ( MXPROB = 9, MXCASE = 8, MAXTAU = 5, MAXMU = 10,
     &            MAXAZ = 3 )
c     ..
c     .. Local Scalars ..

      INTEGER   I, J, K
c     ..
c     .. Common blocks ..

      COMMON    / DOCHEK / TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU

      REAL      TSTDFD( MAXTAU, MXCASE, MXPROB ),
     &          TSTFDN( MAXTAU, MXCASE, MXPROB ),
     &          TSTFIR( MAXTAU, MXCASE, MXPROB ),
     &          TSTFUP( MAXTAU, MXCASE, MXPROB ),
     &          TSTUU( MAXTAU, MAXMU, MAXAZ, MXCASE, MXPROB )
c     ..

c ********************* Test Case 1A *********************************

      DATA (TSTFIR(I,1,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,1,1), I = 1, 2) / 0.0, 7.94108E-02 /
      DATA (TSTFUP(I,1,1), I = 1, 2) / 7.99451E-02, 0.0 /
      DATA (TSTDFD(I,1,1), I = 1, 2) / 2.54067E+01, 1.86531E+01 /
      DATA ((TSTUU(I,J,1,1,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.17771E-01, 2.64170E-02, 1.34041E-02,
     &    1.33826E-02, 2.63324E-02, 1.15898E-01, 3*0.0 /

c ********************* Test Case 1B *********************************

      DATA (TSTFIR(I,2,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,2,1), I = 1, 2) / 0.0, 4.20233E-01 /
      DATA (TSTFUP(I,2,1), I = 1, 2) / 4.22922E-01, 0.0 /
      DATA (TSTDFD(I,2,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 6.22884E-01, 1.39763E-01, 7.09192E-02,
     &    7.08109E-02, 1.39337E-01, 6.13458E-01, 3*0.0 /

c ********************* Test Case 1C *********************************

      DATA (TSTFIR(I,3,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,3,1), I = 1, 2) / 3.14159E+00, 3.04897E+00 /
      DATA (TSTFUP(I,3,1), I = 1, 2) / 9.06556E-02, 0.0 /
      DATA (TSTDFD(I,3,1), I = 1, 2) / 6.66870E-02, 5.88936E-02 /
      DATA ((TSTUU(I,J,1,3,1), J = 1, 6), I = 1, 2)
     &  / 3*1.0, 1.33177E-01, 2.99879E-02, 1.52233E-02,
     &    9.84447E-01, 9.69363E-01, 8.63946E-01, 3*0.0 /

c ********************* Test Case 1D *********************************

      DATA (TSTFIR(I,4,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,4,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,4,1), I = 1, 2) / 2.59686E-01, 0.0 /
      DATA (TSTDFD(I,4,1), I = 1, 2) / 2.57766E+01, 0.0 /
      DATA ((TSTUU(I,J,1,4,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 2.62972E-01, 9.06967E-02, 5.02853E-02,
     &    1.22980E-15, 1.30698E-17, 6.88841E-18, 3*0.0 /

c ********************* Test Case 1E *********************************

      DATA (TSTFIR(I,5,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,5,1), I = 1, 2) / 0.0, 6.76954E-02 /
      DATA (TSTFUP(I,5,1), I = 1, 2) / 3.07390E+00, 0.0 /
      DATA (TSTDFD(I,5,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,5,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.93321E+00, 1.02732E+00, 7.97199E-01,
     &    2.71316E-02, 1.87805E-02, 1.16385E-02, 3*0.0 /

c ********************* Test Case 1F *********************************

      DATA (TSTFIR(I,6,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,6,1), I = 1, 2) / 3.14159E+00, 4.60048E-03 /
      DATA (TSTFUP(I,6,1), I = 1, 2) / 2.49618E+00, 0.0 /
      DATA (TSTDFD(I,6,1), I = 1, 2) / 1.14239E-01, 7.93633E-05 /
      DATA ((TSTUU(I,J,1,6,1), J = 1, 6), I = 1, 2)
     &  / 3*1.0, 8.77510E-01, 8.15136E-01, 7.52715E-01,
     &    1.86840E-03, 1.26492E-03, 7.79281E-04, 3*0.0 /


c ********************* Test Case 2A *********************************

      DATA (TSTFIR(I,1,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,1,2), I = 1, 2) / 0.0, 4.41791E-02 /
      DATA (TSTFUP(I,1,2), I = 1, 2) / 5.35063E-02, 0.0 /
      DATA (TSTDFD(I,1,2), I = 1, 2) / 1.66570E+00, 1.89848E-01 /
      DATA ((TSTUU(I,J,1,1,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.61796E-01, 2.11501E-02, 7.86713E-03,
     &    7.71897E-03, 2.00778E-02, 2.57685E-02, 3*0.0 /

c ********************* Test Case 2B *********************************

      DATA (TSTFIR(I,2,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,2,2), I = 1, 2) / 0.0, 1.06123E-01 /
      DATA (TSTFUP(I,2,2), I = 1, 2) / 1.25561E-01, 0.0 /
      DATA (TSTDFD(I,2,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 3.47678E-01, 4.87120E-02, 1.89387E-02,
     &    1.86027E-02, 4.64061E-02, 6.77603E-02, 3*0.0 /

c ********************* Test Case 2C *********************************

      DATA (TSTFIR(I,3,2), I = 1, 2) / 2.52716E-01, 2.56077E-28 /
      DATA (TSTFDN(I,3,2), I = 1, 2) / 0.0, 2.51683E-04 /
      DATA (TSTFUP(I,3,2), I = 1, 2) / 6.24730E-02, 0.0 /
      DATA (TSTDFD(I,3,2), I = 1, 2) / 1.67462E+00, 1.75464E-04 /
      DATA ((TSTUU(I,J,1,3,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.62566E-01, 2.45786E-02, 1.01498E-02,
     &    1.70004E-04, 3.97168E-05, 1.32472E-05, 3*0.0 /

c ********************* Test Case 2D *********************************

      DATA (TSTFIR(I,4,2), I = 1, 2) / 2.52716E-01, 2.56077E-28 /
      DATA (TSTFDN(I,4,2), I = 1, 2) / 0.0, 2.68008E-02 /
      DATA (TSTFUP(I,4,2), I = 1, 2) / 2.25915E-01, 0.0 /
      DATA (TSTDFD(I,4,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,4,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 3.64010E-01, 8.26993E-02, 4.92370E-02,
     &    1.05950E-02, 7.69337E-03, 3.79276E-03, 3*0.0 /


c ********************* Test Case 3A *********************************

      DATA (TSTFIR(I,1,3), I = 1, 2) / 3.14159E+00, 1.42628E-04 /
      DATA (TSTFDN(I,1,3), I = 1, 2) / 0.0, 4.70303E-02 /
      DATA (TSTFUP(I,1,3), I = 1, 2) / 1.69951E-01, 0.0 /
      DATA (TSTDFD(I,1,3), I = 1, 2) / 2.59221E+01, 7.43638E-02 /
      DATA ((TSTUU(I,J,1,1,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 4.86923E-01, 5.08802E-02, 1.03351E-02,
     &    7.88518E-03, 2.22081E-02, 2.90003E-03, 3*0.0 /

c ********************* Test Case 3B *********************************

      DATA (TSTFIR(I,2,3), I = 1, 2) / 3.14159E+00, 1.42628E-04 /
      DATA (TSTFDN(I,2,3), I = 1, 2) / 0.0, 1.31455E+00 /
      DATA (TSTFUP(I,2,3), I = 1, 2) / 1.82690E+00, 0.0 /
      DATA (TSTDFD(I,2,3), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 3.72886E+00, 6.32069E-01, 1.51211E-01,
     &    2.13998E-01, 5.84010E-01, 3.77162E-01, 3*0.0 /

c ********************* Test Case 3C *********************************

      DATA (TSTFIR(I,3,3), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,3,3), I = 1, 2) / 3.14159E+00, 2.49660E+00 /
      DATA (TSTFUP(I,3,3), I = 1, 2) / 5.83122E-01, 0.0 /
      DATA (TSTDFD(I,3,3), I = 1, 2) / 8.06945E-02, 4.35243E-02 /
      DATA ((TSTUU(I,J,1,3,3), J = 1, 6), I = 1, 2) /
     &    3*1.0, 5.66618E-01, 2.32213E-01, 7.60274E-02,
     &    9.11175E-01, 7.44031E-01, 4.03607E-01, 3*0.0 /

c ********************* Test Case 3D *********************************

      DATA (TSTFIR(I,4,3), I = 1, 2) / 3.14159E+00, 5.67011E-35 /
      DATA (TSTFDN(I,4,3), I = 1, 2) / 0.0, 1.99813E-05 /
      DATA (TSTFUP(I,4,3), I = 1, 2) / 1.70062E-01, 0.0 /
      DATA (TSTDFD(I,4,3), I = 1, 2) / 2.59222E+01, 1.91941E-05 /
      DATA ((TSTUU(I,J,1,4,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 4.86926E-01, 5.09190E-02, 1.03691E-02,
     &    2.07944E-05, 9.76964E-07, 1.88004E-07, 3*0.0 /

c ********************* Test Case 3E *********************************

      DATA (TSTFIR(I,5,3), I = 1, 2) / 3.14159E+00, 5.67011E-35 /
      DATA (TSTFDN(I,5,3), I = 1, 2) / 0.0, 5.91973E-01 /
      DATA (TSTFUP(I,5,3), I = 1, 2) / 2.54962E+00, 0.0 /
      DATA (TSTDFD(I,5,3), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,5,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 3.85776E+00, 8.68015E-01, 3.79793E-01,
     &    2.36362E-01, 1.64357E-01, 9.15016E-02, 3*0.0 /

c ********************* Test Case 3F *********************************

      DATA (TSTFIR(I,6,3), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,6,3), I = 1, 2) / 3.14159E+00, 1.01167E+00 /
      DATA (TSTFUP(I,6,3), I = 1, 2) / 1.68750E+00, 0.0 /
      DATA (TSTDFD(I,6,3), I = 1, 2) / 1.00459E-01, 1.71385E-02 /
      DATA ((TSTUU(I,J,1,6,3), J = 1, 6), I = 1, 2)
     &  / 3*1.0, 7.57316E-01, 5.87262E-01, 4.33051E-01,
     &    4.16334E-01, 2.75956E-01, 1.51785E-01, 3*0.0 /


c ********************* Test Case 4A *********************************

      DATA (TSTFIR(I,1,4), I = 1, 3)
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,1,4), I = 1, 3)
     &  / 0.0, 1.17401E+00, 1.81264E+00 /
      DATA (TSTFUP(I,1,4), I = 1, 3)
     &  / 1.73223E-01, 1.11113E-01, 0.0 /
      DATA (TSTDFD(I,1,4), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,4), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 9.29865E-02,
     &    6.61234E-02, 3.55365E-02, 2.49435E+00, 1.19070E-01,
     &    1.35138E-01, 1.24070E-01, 4.02884E-02, 1.73581E-02,
     &    3.34666E+00, 2.19624E-01, 1.57002E-01, 3*0.0 /

c ********************* Test Case 4B *********************************

      DATA (TSTFIR(I,2,4), I = 1, 3)
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,2,4), I = 1, 3)
     &  / 0.0, 1.01517E+00, 1.51554E+00 /
      DATA (TSTFUP(I,2,4), I = 1, 3)
     &  / 1.23666E-01, 7.88691E-02, 0.0 /
      DATA (TSTDFD(I,2,4), I = 1, 3)
     &  / 3.43724E-01, 3.52390E-01, 3.19450E-01 /
      DATA ((TSTUU(I,J,1,2,4), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 6.55782E-02,
     &    4.56643E-02, 2.74243E-02, 2.22302E+00, 9.64101E-02,
     &    9.62927E-02, 8.44926E-02, 2.80216E-02, 1.35087E-02,
     &    2.94684E+00, 1.67508E-01, 1.08212E-01, 3*0.0 /

c ********************* Test Case 4C *********************************

      DATA (TSTFIR(I,3,4), I = 1, 3)
     &  / 1.57080E+00, 5.77864E-01, 2.12584E-01 /
      DATA (TSTFDN(I,3,4), I = 1, 3)
     &  / 0.0, 7.02764E-01, 8.03294E-01 /
      DATA (TSTFUP(I,3,4), I = 1, 3)
     &  / 2.25487E-01, 1.23848E-01, 0.0 /
      DATA (TSTDFD(I,3,4), I = 1, 3)
     &  / 3.85003E-01, 3.37317E-01, 2.16403E-01 /
      DATA (((TSTUU(I,J,K,3,4), J = 1, 6), I = 1, 3), K = 1, 3)
     &  / 3*0.0, 8.70016E-01,
     &    2.24760E-01, 2.28322E-02, 4.76042E-02, 3.00258E+00,
     &    1.41140E+00, 6.97400E-01, 1.09065E-01, 9.35116E-03,
     &    8.37538E-02, 2.68792E+00, 8.76317E-01, 6*0.0,
     &    8.86109E-02, 5.76913E-02, 2.28322E-02,
     &    4.76042E-02, 5.82549E-02, 1.04636E-01, 9.15334E-02,
     &    2.95681E-02, 9.35116E-03, 8.37538E-02, 9.43348E-02,
     &    8.95961E-02, 6*0.0, 6.95291E-02,
     &    4.93283E-02, 2.28322E-02, 4.76042E-02, 2.59416E-02,
     &    6.24961E-02, 5.90188E-02, 2.44592E-02, 9.35116E-03,
     &    8.37538E-02, 4.00024E-02, 4.66785E-02, 3*0.0 /


c ********************* Test Case 5A *********************************

      DATA (TSTFIR(I,1,5), I = 1, 3)
     &  / 3.14159E+00, 3.97856E-14, 5.03852E-28 /
      DATA (TSTFDN(I,1,5), I = 1, 3)
     &  / 0.0, 2.24768E+00, 4.79851E-01 /
      DATA (TSTFUP(I,1,5), I = 1, 3)
     &  / 2.66174E+00, 1.76783E+00, 0.0 /
      DATA (TSTDFD(I,1,5), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,5), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 4.57660E-01,
     &    7.68942E-01, 1.03122E+00, 7.53662E-01, 6.96362E-01,
     &    6.50541E-01, 6.27631E-01, 5.81809E-01, 5.24532E-01,
     &    1.95230E-01, 1.31990E-01, 7.20655E-02, 3*0.0 /

c ********************* Test Case 5B *********************************

      DATA (TSTFIR(I,2,5), I = 1, 3)
     &  / 1.28058E-01, 8.67322E-06, 4.47729E-21 /
      DATA (TSTFDN(I,2,5), I = 1, 3)
     &  / 1.74767E+00, 2.33975E-01, 6.38347E-05 /
      DATA (TSTFUP(I,2,5), I = 1, 3)
     &  / 2.70485E-01, 3.74253E-02, 1.02904E-05 /
      DATA (TSTDFD(I,2,5), I = 1, 3)
     &  / 3.10129E-01, 4.52671E-02, 1.25022E-05 /
      DATA ((TSTUU(I,J,1,2,5), J = 1, 6), I = 1, 3)
     &  / 1.50205E+01, 2.18778E-01, 1.36520E-01, 1.14001E-01,
     &    8.71231E-02, 8.55016E-02, 1.85090E-01, 4.92726E-02,
     &    2.65448E-02, 2.02154E-02, 1.29660E-02, 9.51226E-03,
     &    3.41287E-05, 1.39916E-05, 7.47042E-06, 5.65604E-06,
     &    3.58246E-06, 2.57859E-06 /


c ********************* Test Case 6A *********************************

      DATA (TSTFIR(I,1,6), I = 1, 2) / 2*100.0 /
      DATA (TSTFDN(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTDFD(I,1,6), I = 1, 2) / 2*200.0 /
      DATA ((TSTUU(I,J,1,1,6), J = 1, 4), I = 1, 2) / 8*0.0 /

c ********************* Test Case 6B *********************************

      DATA (TSTFIR(I,2,6), I = 1, 3)
     &  / 1.000000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTDFD(I,2,6), I = 1, 3)
     &  / 2.00000E+02, 7.35759E+01, 2.70671E+01 /
      DATA ((TSTUU(I,J,1,2,6), J = 1, 4), I = 1, 3) / 12*0.0 /

c ********************* Test Case 6C *********************************

      DATA (TSTFIR(I,3,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,3,6), I = 1, 3)
     &  / 1.48450E+00, 2.99914E+00, 6.76676E+00 /
      DATA (TSTDFD(I,3,6), I = 1, 3)
     &  / 2.02010E+02, 7.79962E+01, 4.06006E+01 /
      DATA ((TSTUU(I,J,1,3,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 9.77882E-05, 7.92386E-01,
     &    2*0.0, 1.45131E-02, 1.30642E+00,
     &    2*0.0, 2.15393E+00, 2.15393E+00 /

c ********************* Test Case 6D *********************************

      DATA (TSTFIR(I,4,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,4,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,4,6), I = 1, 3)
     &  / 6.19105E-01, 1.33245E+00, 3.47790E+00 /
      DATA (TSTDFD(I,4,6), I = 1, 3)
     &  / 2.00899E+02, 7.57612E+01, 3.65186E+01 /
      DATA ((TSTUU(I,J,1,4,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 5.07347E-05, 2.72396E-01,
     &    2*0.0, 7.52970E-03, 4.49105E-01,
     &    2*0.0, 1.11751E+00, 7.40449E-01 /

c ********************* Test Case 6E *********************************

      DATA (TSTFIR(I,5,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,5,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,5,6), I = 1, 3)
     &  / 8.28836E+01, 1.65215E+02, 3.59895E+02 /
      DATA (TSTDFD(I,5,6), I = 1, 3)
     &  / 3.10552E+02, 3.11018E+02, 6.79533E+02 /
      DATA ((TSTUU(I,J,1,5,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 3.20501E-03, 4.64836E+01,
     &    2*0.0, 4.75666E-01, 7.66386E+01,
     &    2*0.0, 7.05951E+01, 1.26356E+02 /

c ********************* Test Case 6F *********************************

      DATA (TSTFIR(I,6,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,6,6), I = 1, 3)
     &  / 3.21497E+02, 1.42493E+02, 7.05306E+01 /
      DATA (TSTFUP(I,6,6), I = 1, 3)
     &  / 8.52994E+01, 1.70343E+02, 3.72842E+02 /
      DATA (TSTDFD(I,6,6), I = 1, 3)
     &  / 9.57001E+02, 5.29253E+02, 8.07923E+02 /
      DATA ((TSTUU(I,J,1,6,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 3.58884E-03, 4.75119E+01,
     &    6.20697E+01, 6.89532E-01, 5.32631E-01, 7.83338E+01,
     &    3.76472E+01, 4.64603E-03, 7.90494E+01, 1.29151E+02 /

c ********************* Test Case 6G *********************************

      DATA (TSTFIR(I,7,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,7,6), I = 1, 3)
     &  / 3.21497E+02, 3.04775E+02, 3.63632E+02 /
      DATA (TSTFUP(I,7,6), I = 1, 3)
     &  / 3.36115E+02, 4.13977E+02, 4.43050E+02 /
      DATA (TSTDFD(I,7,6), I = 1, 3)
     &  / 5.81342E+02, 1.28621E+02, -1.71119E+02 /
      DATA ((TSTUU(I,J,1,7,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 7.80731E+01, 1.17126E+02,
     &    9.78748E+01, 1.01048E+02, 1.15783E+02, 1.36113E+02,
     &    1.10061E+02, 1.38631E+02, 1.33373E+02, 1.42865E+02 /

c ********************* Test Case 6H *********************************

      DATA (TSTFIR(I,8,6), I = 1, 3)
     &  / 1.00000E+02, 1.35335E+01, 2.06115E-07 /
      DATA (TSTFDN(I,8,6), I = 1, 3)
     &  / 3.21497E+02, 2.55455E+02, 4.43444E+02 /
      DATA (TSTFUP(I,8,6), I = 1, 3)
     &  / 2.37350E+02, 2.61130E+02, 4.56205E+02 /
      DATA (TSTDFD(I,8,6), I = 1, 3)
     &  / 4.23780E+02, 6.19828E+01, -3.17719E+01 /
      DATA ((TSTUU(I,J,1,8,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 7.12616E+01, 7.80737E+01,
     &    8.49992E+01, 7.73186E+01, 7.88310E+01, 8.56424E+01,
     &    1.38631E+02, 1.45441E+02, 1.44088E+02, 1.45547E+02 /


c ********************* Test Case 7A *********************************

      DATA (TSTFIR(I,1,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,1,7), I = 1, 3)
     &  / 3.19830E+02, 3.54099E+02, 3.01334E+02 /
      DATA (TSTFUP(I,1,7), I = 1, 3)
     &  / 4.29572E+02, 4.47018E+02, 5.94576E+02 /
      DATA (TSTDFD(I,1,7), I = 1, 3)
     &  /-8.04270E+01, 2.51589E+02, 7.15964E+02 /
      DATA (((TSTUU(I,J,K,1,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.41637E+02, 1.48416E+02,
     &    1.06346E+02, 1.34218E+02, 1.02514E+02, 1.58864E+02,
     &    9.64160E+01, 8.87832E+01, 1.89259E+02, 1.89259E+02,
     &    1.01805E+02, 1.01805E+02, 1.27778E+02, 1.48416E+02,
     &    1.06346E+02, 1.06238E+02, 9.41347E+01, 1.58864E+02,
     &    9.64160E+01, 7.48649E+01, 1.89259E+02, 1.89259E+02 /

c ********************* Test Case 7B *********************************

      DATA (TSTFIR(I,2,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,2,7), I = 1, 3)
     &  / 3.19830E+02, 3.50555E+02, 2.92063E+02 /
      DATA (TSTFUP(I,2,7), I = 1, 3)
     &  / 3.12563E+02, 2.68126E+02, 3.05596E+02 /
      DATA (TSTDFD(I,2,7), I = 1, 3)
     &  /-1.68356E+02, 1.01251E+02, 4.09326E+02 /
      DATA (((TSTUU(I,J,K,2,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.35840E+02, 9.56588E+01,
     &    1.05967E+02, 1.28779E+02, 9.00044E+01, 8.87624E+01,
     &    9.53651E+01, 7.47553E+01, 9.72743E+01, 9.72743E+01,
     &    1.01805E+02, 1.01805E+02, 1.21980E+02, 9.56588E+01,
     &    1.05967E+02, 1.00799E+02, 8.16247E+01, 8.87624E+01,
     &    9.53651E+01, 6.08370E+01, 9.72743E+01, 9.72743E+01 /

c ********************* Test Case 7C *********************************

      DATA (TSTFIR(I,3,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,7), I = 1, 3)
     &  / 3.19830E+02, 3.53250E+02, 2.98583E+02 /
      DATA (TSTFUP(I,3,7), I = 1, 3)
     &  / 4.07040E+02, 4.11058E+02, 5.30504E+02 /
      DATA (TSTDFD(I,3,7), I = 1, 3)
     &  /-9.85693E+01, 2.17724E+02, 6.23936E+02 /
      DATA (((TSTUU(I,J,K,3,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.40347E+02, 1.40534E+02,
     &    1.06259E+02, 1.32945E+02, 9.92793E+01, 1.48621E+02,
     &    9.61392E+01, 8.45884E+01, 1.54901E+02, 1.76268E+02,
     &    1.01805E+02, 1.01805E+02, 1.26343E+02, 1.40534E+02,
     &    1.06259E+02, 1.04823E+02, 9.02686E+01, 1.48621E+02,
     &    9.61392E+01, 6.97291E+01, 1.37958E+02, 1.76268E+02 /


c ********************* Test Case 8A *********************************

      DATA (TSTFIR(I,1,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,1,8), I = 1, 3) /
     &  1.00000E+00, 7.22235E-01, 5.13132E-01 /
      DATA (TSTFUP(I,1,8), I = 1, 3) / 9.29634E-02, 2.78952E-02, 0.0 /
      DATA (TSTDFD(I,1,8), I = 1, 3) /
     &  1.12474E+00, 6.51821E-01, 5.63361E-01 /
      DATA ((TSTUU(I,J,1,1,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 5.62566E-02, 1.94423E-02,
     &  2.62711E-01, 1.36952E-01, 1.84909E-02, 5.52188E-03,
     &  2.10014E-01, 5.60376E-02, 2*0.0 /

c ********************* Test Case 8B *********************************

      DATA (TSTFIR(I,2,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,2,8), I = 1, 3) /
     &  1.00000E+00, 7.95332E-01, 6.50417E-01 /
      DATA (TSTFUP(I,2,8), I = 1, 3) / 2.25136E-01, 1.26349E-01, 0.0 /
      DATA (TSTDFD(I,2,8), I = 1, 3) /
     &  5.12692E-01, 3.56655E-01, 5.68095E-02 /
      DATA ((TSTUU(I,J,1,2,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 1.23687E-01, 4.95581E-02,
     &  2.77499E-01, 1.83950E-01, 8.35695E-02, 2.50575E-02,
     &  2.40731E-01, 1.29291E-01, 2*0.0 /

c ********************* Test Case 8C *********************************

      DATA (TSTFIR(I,3,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,3,8), I = 1, 3) /
     &  1.00000E+00, 4.86157E-01, 1.59984E-01 /
      DATA (TSTFUP(I,3,8), I = 1, 3) / 3.78578E-01, 2.43397E-01, 0.0 /
      DATA (TSTDFD(I,3,8), I = 1, 3) /
     &  5.65095E-01, 2.76697E-01, 1.35679E-02 /
      DATA ((TSTUU(I,J,1,3,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 1.49335E-01, 1.04766E-01,
     &  1.89020E-01, 9.88158E-02, 9.65192E-02, 6.54445E-02,
     &  6.84762E-02, 2.96698E-02, 2*0.0 /


c ********************* Test Case 9A *********************************

      DATA (TSTFIR(I,1,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,1,9), I = 1, 5) /
     &  1.00000E+00, 3.55151E-01, 1.44265E-01, 6.71445E-03, 6.16968E-07/
      DATA (TSTFUP(I,1,9), I = 1, 5) /
     &  2.27973E-01, 8.75098E-02, 3.61819E-02, 2.19291E-03, 0.0 /
      DATA (TSTDFD(I,1,9), I = 1, 5) /
     &  8.82116E-01, 2.32366E-01, 9.33443E-02, 3.92782E-03, 1.02500E-07/
      DATA ((TSTUU(I,J,1,1,9), J = 1, 4), I = 1, 5) /
     &  2*3.18310E-01, 9.98915E-02, 5.91345E-02,
     &  1.53507E-01, 5.09531E-02, 3.67006E-02, 2.31903E-02,
     &  7.06614E-02, 2.09119E-02, 1.48545E-02, 9.72307E-03,
     &  3.72784E-03, 1.08815E-03, 8.83317E-04, 5.94743E-04,
     &  2.87656E-07, 1.05921E-07, 2*0.0 /

c ********************* Test Case 9B *********************************

      DATA (TSTFIR(I,2,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,2,9), I = 1, 5) /
     &  1.00000E+00, 4.52357E-01, 2.36473E-01, 2.76475E-02, 7.41854E-05/
      DATA (TSTFUP(I,2,9), I = 1, 5) /
     &  1.00079E-01, 4.52015E-02, 2.41941E-02, 4.16017E-03, 0.0 /
      DATA (TSTDFD(I,2,9), I = 1, 5) /
     &  8.04577E-01, 2.55330E-01, 1.30976E-01, 1.36227E-02, 1.22022E-05/
      DATA ((TSTUU(I,J,1,2,9), J = 1, 4), I = 1, 5) /
     &  2*3.18310E-01, 7.39198E-02, 1.32768E-02,
     &  1.96609E-01, 5.92369E-02, 3.00230E-02, 7.05566E-03,
     &  1.15478E-01, 3.01809E-02, 1.52672E-02, 4.06932E-03,
     &  1.46177E-02, 3.85590E-03, 2.38301E-03, 7.77891E-04,
     &  3.37743E-05, 1.20858E-05, 2*0.0 /

c ********************* Test Case 9C *********************************

      DATA (TSTFIR(I,3,9), I = 1, 5) /
     &  1.57080E+00, 1.92354E-01, 2.35550E-02, 9.65131E-06, 9.03133E-19/
      DATA (TSTFDN(I,3,9), I = 1, 5) /
     &  6.09217E+00, 4.97279E+00, 4.46616E+00, 4.22731E+00, 4.73767E+00/
      DATA (TSTFUP(I,3,9), I = 1, 5) /
     &  4.68414E+00, 4.24381E+00, 4.16941E+00, 4.30667E+00, 5.11524E+00/
      DATA (TSTDFD(I,3,9), I = 1, 5) /
     &  3.49563E+00, 8.81206E-01, 3.50053E-01, 1.93471E-02, 7.15349E-02/
      DATA (((TSTUU(I,J,K,3,9), J = 1, 4), I = 1, 5), K = 1, 3)
     & / 2*1.93920E+00, 1.61855E+00, 1.43873E+00,
     &   1.66764E+00, 1.44453E+00, 1.38339E+00, 1.33890E+00,
     &   1.48511E+00, 1.35009E+00, 1.33079E+00, 1.32794E+00,
     &   1.34514E+00, 1.35131E+00, 1.35980E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00,
     &   1.93920E+00, 1.93920E+00, 1.57895E+00, 1.43873E+00,
     &   1.66764E+00, 1.42925E+00, 1.37317E+00, 1.33890E+00,
     &   1.48511E+00, 1.34587E+00, 1.32921E+00, 1.32794E+00,
     &   1.34514E+00, 1.35129E+00, 1.35979E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00,
     &   1.93920E+00, 1.93920E+00, 1.56559E+00, 1.43873E+00,
     &   1.66764E+00, 1.42444E+00, 1.37034E+00, 1.33890E+00,
     &   1.48511E+00, 1.34469E+00, 1.32873E+00, 1.32794E+00,
     &   1.34514E+00, 1.35128E+00, 1.35979E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 2*1.62823E+00 /

      END


