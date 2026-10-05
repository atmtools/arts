      PROGRAM  SCATCNV
C       Converts a scattering file (e.g. from MIESCAT) to DDA radiative
C       transfer format file.  Scattering files have the six unique elements
C       of the Stokes scattering matrix for randomly oriented particles in
C       the form of a Legendre series in cosine of the scattering angle.
C       DDA radiative transfer files have the scattering matrix, extinction,
C       matrix, and emission vector listed and ready for use in a discrete
C       angle radiative transfer code (i.e. RT4).  The Stokes 4x4
C       scattering matrix is listed for the incident and outgoing quadrature
C       angles and for the Fourier modes in azimuth.  The extinction matrix
C       and emission vector are listed for the quadrature incident angles.
C       The guts of the code are taken from RADSCAT3.FOR.
C
      INTEGER  NSTOKES, NUMMU, AZIORDER, NLEGEN
      REAL*8   MU_VALUES(32), QUAD_WEIGHTS(32)
      REAL*8   COEF(6,100),  EXTINCT, ALBEDO, CONST
      CHARACTER*64  SCAT_FILE, OUT_FILE
      CHARACTER*8   QUADTYPE


      NSTOKES = 4

      WRITE (*,'(1X,A)') 'Scattering file name : '
      READ (*,'(A)') SCAT_FILE

      WRITE (*,'(1X,A)') 'Output file name : '
      READ (*,'(A)') OUT_FILE

      WRITE (*,'(1X,A)') 'Number of quadrature angles : '
      READ (*,*) NUMMU

      WRITE (*,'(1X,A)')
     .  'Type of quadrature scheme (Gaussian, Lobatto) : '
      READ (*,'(A)') QUADTYPE

      WRITE (*,'(1X,A)') 'Azimuth order : '
      READ (*,*) AZIORDER



      CALL READ_SCAT_FILE (SCAT_FILE, NLEGEN, COEF, EXTINCT, ALBEDO)

      IF (QUADTYPE(1:1) .EQ. 'L') THEN
        QUADTYPE = 'LOBATTO'
        CALL LOBATTO_QUADRATURE (NUMMU, MU_VALUES, QUAD_WEIGHTS)
      ELSE
        QUADTYPE = 'GAUSSIAN'
        CALL GAUSS_LEGENDRE_QUADRATURE (NUMMU, MU_VALUES, QUAD_WEIGHTS)
      ENDIF


      OPEN (UNIT=2, FILE=OUT_FILE, STATUS='UNKNOWN')
      WRITE (2,11) NUMMU, AZIORDER, QUADTYPE
11    FORMAT (1X,I3,2X,I3,3X,1H',A8,1H')

      CONST = EXTINCT*ALBEDO/ (4.0D0*3.1415926535897932384D0)
      CALL SCATTERING (NUMMU, AZIORDER, NSTOKES, CONST,
     .                     MU_VALUES, NLEGEN, COEF)

      CALL OUTPUT_EXTINCT_EMIS (NUMMU, MU_VALUES, EXTINCT, ALBEDO)

      CLOSE (2)

      END





      SUBROUTINE OUTPUT_EXTINCT_EMIS (NUMMU, MU_VALUES,
     .                                EXTINCT, ALBEDO)
      INTEGER  NUMMU
      REAL*8   MU_VALUES(NUMMU), EXTINCT, ALBEDO
      INTEGER  J, L
      REAL*8   MU, ABSORB

      WRITE (2,*) 'C  EXTINCTION MATRIX'
      DO 100 L = 1, 2
        DO 100 J = 1, NUMMU
          MU = MU_VALUES(J)
          IF (L .EQ. 2) MU = -MU
          WRITE (2,'(1X,F11.8)') MU
          WRITE (2,'(4(2X,E15.8))') EXTINCT, 0.0, 0.0, 0.0
          WRITE (2,'(4(2X,E15.8))') 0.0, EXTINCT, 0.0, 0.0
          WRITE (2,'(4(2X,E15.8))') 0.0, 0.0, EXTINCT, 0.0
          WRITE (2,'(4(2X,E15.8))') 0.0, 0.0, 0.0, EXTINCT
100   CONTINUE

      WRITE (2,*) 'C  EMISSION VECTOR'
      ABSORB = (1.0-ALBEDO)*EXTINCT
      DO 150 L = 1, 2
        DO 150 J = 1, NUMMU
          MU = MU_VALUES(J)
          IF (L .EQ. 2) MU = -MU
          WRITE (2,'(1X,F11.8,2X,4(1X,E15.8))') MU, ABSORB,0.0,0.0,0.0
150   CONTINUE

      RETURN
      END




      SUBROUTINE READ_SCAT_FILE (SCAT_FILE, NLEGEN, COEF,
     .                           EXTINCTION, SS_ALBEDO)
C        READ_SCAT_FILE reads in the scattering file and returns
C      the degree of the Legendre series (NLEGEN), the extinction
C      coefficient, the single scatter albedo, and the Legendre
C      coefficients.  Only coefficients for the six unique elements
C      (for randomly oriented particles with a plane of symmetry)
C      are returned.
      INTEGER NLEGEN
      REAL*8   COEF(6,*), EXTINCTION, SS_ALBEDO
      CHARACTER*(*)  SCAT_FILE
      INTEGER  L, K
      REAL*8   SCATCOEF
      CHARACTER*132  BUFFER

C           Skip over the comment lines at the beginning
      OPEN (UNIT=4, FILE=SCAT_FILE, STATUS='OLD')
100   CONTINUE
          READ (4,'(A)') BUFFER
      IF (BUFFER(1:1) .EQ. 'C') GOTO 100
      BACKSPACE 4

C           Input the extinction, scattering, and albedo
      READ (4,*) EXTINCTION
      READ (4,*) SCATCOEF
      READ (4,*) SS_ALBEDO
      READ (4,*) NLEGEN
C           The Legendre coefficients are input on four lines for
C             each L value.  Each line contains four values.
      DO 120 L = 1, NLEGEN+1
        READ (4,*) K, COEF(1,L), COEF(2,L), COEF(3,L), COEF(4,L),
     .                COEF(5,L), COEF(6,L)
120   CONTINUE
      SS_ALBEDO = SCATCOEF/EXTINCTION

      CLOSE (4)

      RETURN
      END







      SUBROUTINE SCATTERING (NUMMU, AZIORDER, NSTOKES, CONST,
     .                       MU_VALUES, NLEGEN, COEF)
      INTEGER  NSTOKES, NUMMU, AZIORDER, NLEGEN
      REAL*8   MU_VALUES(NUMMU), CONST, COEF(6,1)
      INTEGER  MAXLEG
      PARAMETER (MAXLEG=64)
      INTEGER  I1, I2, I, J1, J2, K, L1, L2, M
      INTEGER  NUMPTS, DOSUM(6)
      REAL*8   C, MU1, MU2,  DELPHI, COS_SCAT
      REAL*8   PHASE_MATRIX(4,4)
      REAL*8   SCAT_MATRIX(4,4,4*MAXLEG), BASIS_MATRIX(4,4,4*MAXLEG)
      REAL*8   ZERO, TWOPI
      PARAMETER (ZERO=0.0D0, TWOPI=2.0D0*3.1415926535897932384D0)


      NUMPTS = 2* 2**INT(LOG(FLOAT(NLEGEN+4))/LOG(2.0)+1.0)
      IF (AZIORDER .EQ. 0)  NUMPTS = 2*INT((NLEGEN+1)/2) + 4

C           Find how many Legendre series must be summed
      CALL NUMBER_SUMS (NSTOKES, NLEGEN, COEF, DOSUM)

      DO 50 I = 1, NLEGEN+1
        DO 50 K = 1, 6
          COEF(K,I) = CONST*COEF(K,I)
50    CONTINUE

      WRITE (2,*) 'C  SCATTERING MATRIX'

C       MU1 is the incoming direction, and MU2 is the outgoing direction.
      DO 150 L1 = 1, 2
        DO 150 J1 = 1, NUMMU
          MU1 = MU_VALUES(J1)
          IF (L1 .EQ. 2)  MU1 = -MU1
          DO 150 L2 = 1, 2
            DO 150 J2 = 1, NUMMU
              MU2 = MU_VALUES(J2)
              IF (L2 .EQ. 2)  MU2 = -MU2
C                   Only need to calculate phase matrix for half of
C                     the delphi's, the rest come from symmetry.
              DO 100 K = 1, NUMPTS/2 + 1
                DELPHI = (TWOPI*(K-1))/NUMPTS
                COS_SCAT = MU1*MU2 + DSQRT((1.-MU1**2)*(1.-MU2**2))*
     .                            DCOS(DELPHI)
                CALL SUM_LEGENDRE (NLEGEN, COEF, COS_SCAT,
     .                           DOSUM, PHASE_MATRIX)
                CALL ROTATE_PHASE_MATRIX (PHASE_MATRIX, MU1, MU2,
     .                 DELPHI, COS_SCAT, SCAT_MATRIX(1,1,K), NSTOKES)
                CALL MATRIX_SYMMETRY (NSTOKES, SCAT_MATRIX(1,1,K),
     .                              SCAT_MATRIX(1,1,NUMPTS-K+2) )
100           CONTINUE
            CALL FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES,
     .                   SCAT_MATRIX, BASIS_MATRIX)

111         FORMAT (1X,F11.8,2X,F11.8,3X,I3,1X,A)
113         FORMAT (1X,4(1X,E15.8))
            WRITE (2,111) MU1, MU2, 0, ' '
            DO 110 I2 = 1, NSTOKES
                WRITE (2,113) (BASIS_MATRIX(I2,I1,1), I1=1,NSTOKES)
110         CONTINUE
            DO 140 M = 1, AZIORDER
              WRITE (2,111) MU1, MU2, M, 'C'
              DO 120 I2 = 1, NSTOKES
                WRITE (2,113) (BASIS_MATRIX(I2,I1,M+1), I1=1,NSTOKES)
120           CONTINUE
              WRITE (2,111) MU1, MU2, M, 'S'
              DO 130 I2 = 1, NSTOKES
                WRITE (2,113)
     .               (BASIS_MATRIX(I2,I1,M+AZIORDER+1), I1=1,NSTOKES)
130           CONTINUE
140         CONTINUE

150   CONTINUE

      RETURN
      END









      SUBROUTINE NUMBER_SUMS (NSTOKES, NLEGEN, COEF, DOSUM)
      INTEGER  NSTOKES, NLEGEN, DOSUM(6)
      REAL*8   COEF(6,*)
      INTEGER  I, SCAT, CASE, SUMCASES(6,5)
      REAL*8   ZERO
      PARAMETER (ZERO=0.0D0)
      DATA     SUMCASES/ 1,0,0,0,0,0,  1,1,1,0,0,0, 1,1,1,1,0,0,
     .                   1,1,1,0,1,0,  1,1,1,1,1,1/

C           SCAT: 1 is Rayleigh, 2 is Mie, 3 is general
      SCAT = 1
      DO 100 I = 1, NLEGEN+1
         IF (COEF(4,I) .NE. ZERO)  SCAT = 2
100   CONTINUE
      DO 110 I = 1, NLEGEN+1
         IF (COEF(1,I) .NE. COEF(5,I) .OR.
     .       COEF(3,I) .NE. COEF(6,I))  SCAT = 3
110   CONTINUE

      IF (NSTOKES .EQ. 1) THEN
          CASE = 1
      ELSE IF (NSTOKES .LE. 3) THEN
          CASE = 2
          IF (SCAT .EQ. 3)  CASE = 4
      ELSE
          CASE = 2
          IF (SCAT .EQ. 2)  CASE = 3
          IF (SCAT .EQ. 3)  CASE = 5
      ENDIF
      DO 200 I = 1, 6
          DOSUM(I) = SUMCASES(I,CASE)
200   CONTINUE

      RETURN
      END



      SUBROUTINE SUM_LEGENDRE (NLEGEN, COEF, X, DOSUM, PHASE_MATRIX)
C       SUM_LEGENDRE sums the Legendre series for each element of the
C       phase matrix using X for the scattering angle.  There are
C       only six sets of Legendre coefficients in COEF since there are
C       six independant parameters for randomly oriented particles
C       with a plane of symmetry.  The constant array DOSUM controls which
C       of the six series are summed.  The ROW and COL arrays give the
C       location of each of the series in the phase matrix.
      INTEGER  NLEGEN, DOSUM(6)
      REAL*8   COEF(6,1), X, PHASE_MATRIX(4,4)
      INTEGER  I, J, K, L, M
      INTEGER  ROW(6), COL(6)
      REAL*8   SUM, PL, PL1, PL2
      DATA     ROW/1,1,3,3,2,4/, COL/1,2,3,4,2,4/

C           Sum the Legendre series
      DO 120 I = 1, 6
          SUM = 0.0
          IF (DOSUM(I) .EQ. 1) THEN
            PL1 = 1.0
            PL = 1.0
            DO 100 L = 0, NLEGEN
              M = L + 1
              IF (L .GT. 0)  PL = (2*L-1)*X*PL1/L - (L-1)*PL2/L
              SUM = SUM + COEF(I,M)*PL
              PL2 = PL1
              PL1 = PL
100         CONTINUE
          ENDIF
          PHASE_MATRIX(ROW(I),COL(I)) = SUM
120   CONTINUE
      PHASE_MATRIX(2,1) = PHASE_MATRIX(1,2)
      PHASE_MATRIX(4,3) = -PHASE_MATRIX(3,4)
      IF (DOSUM(5) .EQ. 0)  PHASE_MATRIX(2,2) = PHASE_MATRIX(1,1)
      IF (DOSUM(6) .EQ. 0)  PHASE_MATRIX(4,4) = PHASE_MATRIX(3,3)

      RETURN
      END






      SUBROUTINE ROTATE_PHASE_MATRIX (PHASE_MATRIX1, MU1, MU2,
     .                 DELPHI, COS_SCAT, PHASE_MATRIX2, NSTOKES)
C        ROTATE_PHASE_MATRIX applies the rotation of the polarization
C      basis from the incident plane into the scattering plane and
C      from the scattering plane to the outgoing plane.
C      MU1 is the incoming direction, and MU2 is the outgoing direction.
C      Currently, set up for a phase matrix for randomly oriented particles
C      with a plane of symmetry - only 6 unique parameters.
      INTEGER  NSTOKES
      REAL*8   PHASE_MATRIX1(4,4), PHASE_MATRIX2(4,4)
      REAL*8   MU1, MU2, DELPHI, COS_SCAT
      REAL*8   SIN_SCAT, SIN_THETA1, SIN_THETA2
      REAL*8   SINPHI, COSPHI, SIN1, SIN2, COS1, COS2
      REAL*8   SIN21, COS21, SIN22, COS22, A1, A2, A3, A4, B1, B2
      REAL*8   ZERO
      PARAMETER (ZERO=0.0D0)

      A1 = PHASE_MATRIX1(1,1)
      PHASE_MATRIX2(1,1) = A1
      IF (NSTOKES .EQ. 1) RETURN

      SIN_SCAT = DSQRT(1.-COS_SCAT**2)
      SIN_THETA1 = DSQRT(1.-MU1**2)
      SIN_THETA2 = DSQRT(1.-MU2**2)
      SINPHI = DSIN(DELPHI)
      COSPHI = DCOS(DELPHI)
      IF (SIN_SCAT .EQ. ZERO) THEN
          SIN1 = 0.0
          SIN2 = 0.0
          COS1 = 1.0
          COS2 = -1.0
      ELSE
          SIN1 = SIN_THETA2*SINPHI /SIN_SCAT
          SIN2 = SIN_THETA1*SINPHI /SIN_SCAT
          COS1 =  (SIN_THETA1*MU2 - SIN_THETA2*MU1*COSPHI)/SIN_SCAT
          COS2 =  (SIN_THETA2*MU1 - SIN_THETA1*MU2*COSPHI)/SIN_SCAT
      ENDIF
      SIN21 = 2.0*SIN1*COS1
      COS21 = 1.0 - 2.0*SIN1**2
      SIN22 = 2.0*SIN2*COS2
      COS22 = 1.0 - 2.0*SIN2**2

      IF (NSTOKES .GT. 1) THEN
        A2 = PHASE_MATRIX1(2,2)
        A3 = PHASE_MATRIX1(3,3)
        B1 = PHASE_MATRIX1(1,2)
        PHASE_MATRIX2(1,2) = B1*COS21
        PHASE_MATRIX2(2,1) = B1*COS22
        PHASE_MATRIX2(2,2) = A2*COS21*COS22 - A3*SIN21*SIN22
      ENDIF
      IF (NSTOKES .GT. 2) THEN
        PHASE_MATRIX2(1,3) = -B1*SIN21
        PHASE_MATRIX2(2,3) = -A2*SIN21*COS22 - A3*COS21*SIN22
        PHASE_MATRIX2(3,1) = B1*SIN22
        PHASE_MATRIX2(3,2) = A2*COS21*SIN22 + A3*SIN21*COS22
        PHASE_MATRIX2(3,3) = -A2*SIN21*SIN22 + A3*COS21*COS22
      ENDIF
      IF (NSTOKES .GT. 3) THEN
        A4 = PHASE_MATRIX1(4,4)
        B2 = PHASE_MATRIX1(3,4)
        PHASE_MATRIX2(1,4) = 0.0
        PHASE_MATRIX2(2,4) = -B2*SIN22
        PHASE_MATRIX2(3,4) = B2*COS22
        PHASE_MATRIX2(4,1) = 0.0
        PHASE_MATRIX2(4,2) = -B2*SIN21
        PHASE_MATRIX2(4,3) = -B2*COS21
        PHASE_MATRIX2(4,4) = A4
      ENDIF

      RETURN
      END



      SUBROUTINE MATRIX_SYMMETRY (NSTOKES, MATRIX1, MATRIX2)
C        MATRIX_SYMMETRY performs a symmetry operation on a
C      phase matrix.  The operation consists of negating the
C      off-diagonal 2 by 2 blocks.  This operation is equivalent
C      to negating (mu) and (mu') or negating (phi'-phi).
      INTEGER  NSTOKES
      REAL*8   MATRIX1(4,4), MATRIX2(4,4)

      MATRIX2(1,1) = MATRIX1(1,1)
      IF (NSTOKES .GT. 1) THEN
        MATRIX2(1,2) = MATRIX1(1,2)
        MATRIX2(2,1) = MATRIX1(2,1)
        MATRIX2(2,2) = MATRIX1(2,2)
      ENDIF
      IF (NSTOKES .GT. 2) THEN
        MATRIX2(1,3) = -MATRIX1(1,3)
        MATRIX2(2,3) = -MATRIX1(2,3)
        MATRIX2(3,1) = -MATRIX1(3,1)
        MATRIX2(3,2) = -MATRIX1(3,2)
        MATRIX2(3,3) = MATRIX1(3,3)
      ENDIF
      IF (NSTOKES .GT. 3) THEN
        MATRIX2(1,4) = -MATRIX1(1,4)
        MATRIX2(2,4) = -MATRIX1(2,4)
        MATRIX2(3,4) = MATRIX1(3,4)
        MATRIX2(4,1) = -MATRIX1(4,1)
        MATRIX2(4,2) = -MATRIX1(4,2)
        MATRIX2(4,3) = MATRIX1(4,3)
        MATRIX2(4,4) = MATRIX1(4,4)
      ENDIF

      RETURN
      END



      SUBROUTINE FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES,
     .                           REAL_MATRIX, BASIS_MATRIX)
      INTEGER  AZIORDER, NUMPTS, NSTOKES
      REAL*8   REAL_MATRIX(4,4,1), BASIS_MATRIX(4,4,1)
      INTEGER  MAXLEG
      PARAMETER (MAXLEG=64)
      REAL*8   BASIS_VECTOR(4*MAXLEG), REAL_VECTOR(4*MAXLEG)
      INTEGER  NUMAZI, I, J, K, M

      NUMAZI = 2*AZIORDER+1
      DO 150 I = 1, NSTOKES
        DO 150 J = 1, NSTOKES
          DO 100 K = 1, NUMPTS
              REAL_VECTOR(K) = REAL_MATRIX(I,J,K)
100       CONTINUE
          CALL FOURIER_BASIS (NUMAZI, AZIORDER, NUMPTS,
     .                        +1, BASIS_VECTOR, REAL_VECTOR)
          DO 120 M = 1, NUMAZI
              BASIS_MATRIX(I,J,M) = BASIS_VECTOR(M)
120       CONTINUE
150   CONTINUE

      RETURN
      END






      SUBROUTINE FOURIER_BASIS (NUMBASIS, ORDER, NUMPTS,
     .                          DIRECTION, BASIS_VECTOR, REAL_VECTOR)
C       FOURIER_BASIS converts a vector between the Fourier basis
C      and azimuth space.  The Fourier basis functions are:
C      { 1, cos(x), cos(2x), . ., cos(Mx), sin(x), sin(2x), . . sin(Mx) }
C      where M is the order of the Fourier basis (ORDER).  The number
C      of elements in the basis in NUMBASIS.  Generally, NUMBASIS=2*ORDER+1,
C      but for even functions NUMBASIS=ORDER+1
C      DIRECTION is negative for conversion to real space, and positive
C      for conversion into Fourier basis.
C      NUMPTS is the number of points in the real space (must be power of 2).
      INTEGER  NUMBASIS, NUMPTS, ORDER, DIRECTION
      REAL*8   BASIS_VECTOR(*), REAL_VECTOR(*)
      INTEGER  I, BASISLEN
      REAL*8   SUM

      BASISLEN = MIN(ORDER,NUMPTS/2-1)
      IF (DIRECTION .LT. 0) THEN
          IF (ORDER .EQ. 0) THEN
            DO 90 I = 1, NUMPTS
              REAL_VECTOR(I) = BASIS_VECTOR(1)
90          CONTINUE
          ELSE
            DO 100 I = 1, NUMPTS
              REAL_VECTOR(I) = 0.0
100         CONTINUE
            REAL_VECTOR(1) = BASIS_VECTOR(1)
            DO 110 I = 1, BASISLEN
              REAL_VECTOR(2*I+1) = BASIS_VECTOR(I+1)/2
110         CONTINUE
            IF (NUMBASIS .GT. ORDER+1) THEN
              DO 120 I = 1, BASISLEN
                  REAL_VECTOR(2*I+2) = BASIS_VECTOR(I+ORDER+1)/2
120           CONTINUE
            ENDIF
            CALL FFT1DR (REAL_VECTOR, NUMPTS, -1)
          ENDIF

      ELSE

          IF (ORDER .EQ. 0) THEN
            SUM = 0.0
            DO 200 I = 1, NUMPTS
              SUM = SUM + REAL_VECTOR(I)
200         CONTINUE
            BASIS_VECTOR(1) = SUM/NUMPTS
          ELSE
            CALL FFT1DR (REAL_VECTOR, NUMPTS, +1)
            BASIS_VECTOR(1) = REAL_VECTOR(1)/NUMPTS
            DO 210 I = 1, BASISLEN
                BASIS_VECTOR(I+1) = 2*REAL_VECTOR(2*I+1)/NUMPTS
210         CONTINUE
            DO 215 I = BASISLEN+1, ORDER
                BASIS_VECTOR(I+1) = 0.0D0
215         CONTINUE
            IF (NUMBASIS .GT. ORDER+1) THEN
                DO 220 I = 1, BASISLEN
                  BASIS_VECTOR(I+ORDER+1) = 2*REAL_VECTOR(2*I+2)/NUMPTS
220             CONTINUE
                DO 225 I = BASISLEN+1, ORDER
                    BASIS_VECTOR(I+ORDER+1) = 0.0D0
225             CONTINUE
            ENDIF
          ENDIF
      ENDIF

      RETURN
      END



      SUBROUTINE FFT1DR (DATA, N, ISIGN)
C        Real 1D FFT.  N must be a power of two.
C      If ISIGN=+1 then real to complex conjugate FFT is done.
C      If ISIGN=-1 then complex conjugate to real FFT is done.
C      The Nyquist frequency component is returned in the first imaginary
C      element.   No normalization (by N) is performed.
      INTEGER N, ISIGN
      REAL*8  DATA(*)
      INTEGER MAXN
      PARAMETER (MAXN=512)
      REAL*8  PHASE(4*MAXN)
      INTEGER MN
      REAL*8  NYQUIST(2)
      SAVE    MN, PHASE

      IF (MN .LT. N) THEN
          MN = N
          IF (MN .GT. MAXN)  STOP 'Phase array too small' 
          CALL MAKEPHASE (PHASE, MN)
      ENDIF

      IF (ISIGN .GT. 0) THEN
C           Forward transform:  real to complex-conjugate
        CALL FFTC (DATA, N/2, PHASE(1))
        CALL FIXREAL (DATA, NYQUIST, N/2, +1, PHASE(1))
        DATA(2) = NYQUIST(1)
      ELSE
C           Inverse transform:  complex-conjugate to real
        NYQUIST(1) = DATA(2)
        CALL FIXREAL (DATA, NYQUIST, N/2, -1, PHASE(2*MN+1))
        CALL FFTC (DATA, N/2, PHASE(2*MN+1))
      ENDIF

      RETURN
      END





      SUBROUTINE FFTC (DATA,N,PHASE)
      INTEGER N
      REAL*8  DATA(*)
      REAL*8  PHASE(*)
      INTEGER I, IREV, M, J, K, C, M0, M1
      INTEGER JMAX, POWER, IPH
      REAL*8  TMPR, TMPI, PHR, PHI

      IF (N .LE. 1) RETURN
      IREV=0
      DO 90 I = 0, N-1
          IF (I .GT. IREV) THEN
              M0 = 2*I+1
              M1 = 2*IREV+1
              TMPR = DATA(M0)
              TMPI = DATA(M0+1)
              DATA(M0) = DATA(M1)
              DATA(M0+1) = DATA(M1+1)
              DATA(M1) = TMPR
              DATA(M1+1) = TMPI
          ENDIF
          M = N
50        CONTINUE
              M = M/2
              IF (IREV.LT.M) GO TO 70
              IREV = IREV - M
          IF (M .GT. 1) GOTO 50
70        IREV = IREV + M
90    CONTINUE


99    CONTINUE
      JMAX = N
      POWER = 1
200   CONTINUE
          JMAX = JMAX/2
          M0 = 1
          M1 = POWER*2+1
          DO 250 J = 1, JMAX
              IPH = 2*POWER
              DO 220 K = 1,POWER
                  PHR = PHASE(IPH-1)
                  PHI = PHASE(IPH)
                  IPH = IPH + 2
                  TMPR = PHR*DATA(M1)-PHI*DATA(M1+1)
                  TMPI = PHI*DATA(M1)+PHR*DATA(M1+1)
                  DATA(M1) = DATA(M0) - TMPR
                  DATA(M1+1) = DATA(M0+1) - TMPI
                  DATA(M0) = DATA(M0) + TMPR
                  DATA(M0+1) = DATA(M0+1) + TMPI
                  M0 = M0 + 2
                  M1 = M1 + 2
220           CONTINUE
              M0 = M0 + POWER*2
              M1 = M1 + POWER*2
250       CONTINUE
          POWER = 2*POWER
      IF (JMAX .GT. 1) GOTO 200

      RETURN
      END




      SUBROUTINE MAKEPHASE (PHASE, NMAX)
      INTEGER NMAX
      REAL*8  PHASE(*)
      INTEGER I, J, N
      DOUBLE PRECISION PI, F
      PARAMETER (PI=3.1415926535897932D0)

      J = 1
      N = 1
100   CONTINUE
        F = PI/N
        DO 110 I = 0, N-1
          PHASE(J) = DCOS(F*DFLOAT(I))
          PHASE(J+1) = DSIN(F*DFLOAT(I))
          J = J + 2
110     CONTINUE
        N = 2*N
      IF (N .LT. NMAX) GOTO 100

      J = 2*NMAX+1
      N = 1
200   CONTINUE
        F = -PI/N
        DO 210 I = 0, N-1
          PHASE(J) = DCOS(F*DFLOAT(I))
          PHASE(J+1) = DSIN(F*DFLOAT(I))
          J = J + 2
210     CONTINUE
        N = 2*N
      IF (N .LT. NMAX) GOTO 200
      RETURN
      END




      SUBROUTINE FIXREAL (DATA, NYQUIST, N, ISIGN, PHASE)
      INTEGER N, ISIGN
      REAL*8  DATA(*), NYQUIST(2)
      REAL*8  PHASE(*)
      INTEGER I, IPH, M, MC
      REAL*8  TMP0R, TMP0I, TMP1R, TMP1I, PHR, PHI

      IPH = 2*N+1
      M = 3
      MC = 2*N-1
      IF (ISIGN .GT. 0) THEN
        NYQUIST(1) = DATA(1)-DATA(2)
        NYQUIST(2) = 0
        DATA(1) = DATA(1)+DATA(2)
        DATA(2) = 0
        DO 100 I = 2, N/2+1
          PHR = PHASE(IPH)
          PHI = PHASE(IPH+1)
          IPH = IPH + 2
          TMP0R = DATA(M)+DATA(MC)
          TMP0I = DATA(M+1)-DATA(MC+1)
          TMP1R = -PHI*(DATA(M)-DATA(MC)) - PHR*(DATA(M+1)+DATA(MC+1))
          TMP1I = PHR*(DATA(M)-DATA(MC)) - PHI*(DATA(M+1)+DATA(MC+1))
          DATA(M) = 0.5*(TMP0R-TMP1R)
          DATA(M+1) = 0.5*(TMP0I-TMP1I)
          DATA(MC) = 0.5*(TMP0R+TMP1R)
          DATA(MC+1) = -0.5*(TMP0I+TMP1I)
          M = M + 2
          MC = MC - 2
100     CONTINUE
      ELSE
        DATA(2) = DATA(1) - NYQUIST(1)
        DATA(1) = DATA(1) + NYQUIST(1)
        DO 200 I = 2, N/2+1
          PHR = PHASE(IPH)
          PHI = PHASE(IPH+1)
          IPH = IPH + 2
          TMP0R = DATA(M)+DATA(MC)
          TMP0I = DATA(M+1)-DATA(MC+1)
          TMP1R = PHI*(DATA(M)-DATA(MC)) + PHR*(DATA(M+1)+DATA(MC+1))
          TMP1I = -PHR*(DATA(M)-DATA(MC)) + PHI*(DATA(M+1)+DATA(MC+1))
          DATA(M) = TMP0R-TMP1R
          DATA(M+1) = TMP0I-TMP1I
          DATA(MC) = TMP0R+TMP1R
          DATA(MC+1) = -(TMP0I+TMP1I)
          M = M + 2
          MC = MC - 2
200     CONTINUE
      ENDIF

      RETURN
      END







      SUBROUTINE GAUSS_LEGENDRE_QUADRATURE
     .                          (NUM, ABSCISSAS, WEIGHTS)
C        Generates the abscissas and weights for an even 2*NUM point
C      Gauss-Legendre quadrature.  Only the NUM positive points are returned.
      INTEGER  NUM
      REAL*8   ABSCISSAS(1), WEIGHTS(1)
      INTEGER  N, I, J, L
      REAL*8   X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-14)

      N = 2*NUM
      DO 130 J = 1, NUM
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO 120 L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
120       CONTINUE
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        ABSCISSAS(NUM+1-J) = X
        WEIGHTS(NUM+1-J) = 2/((1-X*X)*DPL*DPL)
130   CONTINUE

      RETURN
      END




      SUBROUTINE LOBATTO_QUADRATURE
     .                          (NUM, ABSCISSAS, WEIGHTS)
C        Generates the abscissas and weights for an even 2*NUM point
C      Gauss-Legendre quadrature.  Only the NUM positive points are returned.
      INTEGER  NUM
      REAL*8   ABSCISSAS(*), WEIGHTS(*)
      INTEGER  N, N1, I, J, L
      REAL*8   X, XP, PL, PL1, PL2, DPL, D2PL, CI, TINY
      PARAMETER (TINY=3.0D-14)

      N = 2*NUM
      N1 = N-1
      CI = 0.50
      IF (MOD(N,2) .EQ. 1) CI = 1.00
      DO 130 J = 1, NUM-1
        X = SIN(3.141592654*(J-CI)/(N-.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO 120 L = 2, N1
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
120       CONTINUE
          DPL = N1*(X*PL-PL1)/(X*X-1)
          D2PL = (2.D0*X*DPL-N1*(N1+1)*PL) / (1D0-X*X)
          XP = X
          X = XP - DPL/D2PL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        ABSCISSAS(J) = X
        WEIGHTS(J) = 2.0D0/(N*N1*PL*PL)
130   CONTINUE
      ABSCISSAS(NUM) = 1.D0
      WEIGHTS(NUM) = 2.D0/(N*N1)

      RETURN
      END



