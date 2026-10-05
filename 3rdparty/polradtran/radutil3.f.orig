 
 
      SUBROUTINE LAMBERT_SURFACE (NSTOKES, NUMMU, MODE,
     .                            MU_VALUES, QUAD_WEIGHTS,
     .                            GROUND_ALBEDO,
     .                            REFLECT, TRANS, SOURCE)
C         LAMBERT_SURFACE makes the reflection matrix for a Lambertian
C      surface with albedo (GROUND_ALBEDO).  Also makes the diagonal
C      transmission and zero source function.
      INTEGER NSTOKES, NUMMU, MODE
      REAL*8  MU_VALUES(NUMMU), QUAD_WEIGHTS(NUMMU)
      REAL*8  GROUND_ALBEDO
      REAL*8  REFLECT(NSTOKES,NUMMU, NSTOKES,NUMMU, 2)
      REAL*8  TRANS(NSTOKES,NUMMU, NSTOKES,NUMMU, 2)
      REAL*8  SOURCE(NSTOKES,NUMMU, 2)
      INTEGER J1, J2, N
 
 
      N = NSTOKES*NUMMU
      CALL MZERO (2*N, N, REFLECT)
      CALL MZERO (2*N, 1, SOURCE)
      CALL MIDENTITY (N, TRANS(1,1,1,1,1))
      CALL MIDENTITY (N, TRANS(1,1,1,1,2))
C           The Lambertian ground reflects the flux equally in all direction
C             and completely unpolarizes the radiation
      IF (MODE .EQ. 0) THEN
        DO J1 = 1, NUMMU
          DO J2 = 1, NUMMU
            REFLECT(1,J1, 1,J2, 2) = 2.0*GROUND_ALBEDO*
     .                           MU_VALUES(J2)*QUAD_WEIGHTS(J2)
          ENDDO
        ENDDO
      ENDIF
 
      RETURN
      END
 




      SUBROUTINE LAMBERT_RADIANCE (NSTOKES, NUMMU, MODE,
     .                 SRC_CODE, GROUND_ALBEDO, GROUND_TEMP,
     .                 WAVELENGTH, DIRECT_SFC_FLUX, RADIANCE)
C        LAMBERT_RADIANCE calculates the ground radiance of a Lambertian
C      surface.  The radiance is due to thermal radiation and reflected
C      direct solar radiation if there is a solar source (SRC_CODE=1,3).
C      TOTAL_DEPTH is the total optical depth in the direction of the
C      direct source (sun), with the path length effect already included.
      INTEGER NSTOKES, NUMMU, MODE, SRC_CODE
      REAL*8  GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH
      REAL*8  DIRECT_SFC_FLUX
      REAL*8  RADIANCE(NSTOKES,NUMMU)
      INTEGER J, N
      REAL*8  TMP, THERMAL, PLANCK, PI
      PARAMETER (PI = 3.1415926535897932384D0)

      N = NSTOKES*NUMMU
      CALL MZERO (N, 1, RADIANCE)
      IF (MODE .EQ. 0) THEN
C           Thermal radiation going up
        IF (SRC_CODE .EQ. 2 .OR. SRC_CODE .EQ. 3) THEN
          CALL PLANCK_FUNCTION (GROUND_TEMP, 'R', WAVELENGTH, PLANCK)
          THERMAL = (1.0-GROUND_ALBEDO)*PLANCK
          DO J = 1, NUMMU
            RADIANCE(1,J) = THERMAL
          ENDDO
        ENDIF

C           Direct solar reflection (unpolarized)
        IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
          TMP = DIRECT_SFC_FLUX * GROUND_ALBEDO/PI
          DO J = 1, NUMMU
            RADIANCE(1,J) = RADIANCE(1,J) + TMP
          ENDDO
        ENDIF
      ENDIF

      RETURN
      END





      SUBROUTINE FRESNEL_SURFACE (NSTOKES, NUMMU,
     .                            MU_VALUES, INDEX, 
     .                            REFLECT, TRANS, SOURCE)
C         FRESNEL_REFLECT makes the reflection matrix for a
C      plane surface with index of refraction (INDEX) using
C      the Fresnel reflection formulae.  Also makes the diagonal
C      transmission and zero source function.
      INTEGER NSTOKES, NUMMU
      REAL*8  MU_VALUES(NUMMU)
      REAL*8  REFLECT(NSTOKES,NUMMU, NSTOKES,NUMMU, 2)
      REAL*8  TRANS(NSTOKES,NUMMU, NSTOKES,NUMMU, 2)
      REAL*8  SOURCE(NSTOKES,NUMMU, 2)
      COMPLEX*16  INDEX
      INTEGER J, N
      REAL*8      COSI, R1, R2, R3, R4
      COMPLEX*16  EPSILON, D, RH, RV

      N = NSTOKES*NUMMU
      CALL MZERO (2*N, N, REFLECT)
      CALL MZERO (2*N, 1, SOURCE)
      CALL MIDENTITY (N, TRANS(1,1,1,1,1))
      CALL MIDENTITY (N, TRANS(1,1,1,1,2))

      EPSILON = INDEX**2
      DO J = 1, NUMMU
          COSI = MU_VALUES(J)
          D = CDSQRT(EPSILON - 1.0D0 + COSI**2)
          RH = (COSI - D) / (COSI + D)
          RV =(EPSILON*COSI - D) / (EPSILON*COSI + D)
          R1 = (CDABS(RV)**2 + CDABS(RH)**2 )/2.0D0
          R2 = (CDABS(RV)**2 - CDABS(RH)**2 )/2.0D0
          R3 = DREAL(RV*CONJG(RH))
          R4 = DIMAG(RV*CONJG(RH))
          REFLECT(1,J, 1,J, 2) = R1
          IF (NSTOKES .GT. 1) THEN
            REFLECT(1,J, 2,J, 2) = R2
            REFLECT(2,J, 1,J, 2) = R2
            REFLECT(2,J, 2,J, 2) = R1
          ENDIF
          IF (NSTOKES .GT. 2) THEN
            REFLECT(3,J, 3,J, 2) = R3
          ENDIF
          IF (NSTOKES .GT. 3) THEN
            REFLECT(3,J, 4,J, 2) = -R4
            REFLECT(4,J, 3,J, 2) = R4
            REFLECT(4,J, 4,J, 2) = R3
          ENDIF
      ENDDO

      RETURN
      END





      SUBROUTINE FRESNEL_RADIANCE (NSTOKES, NUMMU, MODE,
     .                   MU_VALUES, INDEX, GROUND_TEMP,
     .                   WAVELENGTH, RADIANCE)
C        FRESNEL_RADIANCE calculates the ground radiance of a plane
C      surface using the Fresnel formulae.  The radiance is due only to
C      thermal radiation (this subroutine cannot do specular reflection).
      INTEGER NSTOKES, NUMMU, MODE
      REAL*8  GROUND_TEMP, WAVELENGTH, MU_VALUES(NUMMU)
      REAL*8  RADIANCE(NSTOKES,NUMMU)
      COMPLEX*16  INDEX
      INTEGER J, N
      REAL*8  ZERO,  PLANCK, COSI, R1, R2
      COMPLEX*16  EPSILON, D, RH, RV
      PARAMETER (ZERO=0.0D0)

C           Thermal radiation going up
      N = NSTOKES*NUMMU
      CALL MZERO (N, 1, RADIANCE)
      IF (MODE .EQ. 0) THEN
        CALL PLANCK_FUNCTION (GROUND_TEMP, 'R', WAVELENGTH, PLANCK)
        EPSILON = INDEX**2
        DO J = 1, NUMMU
          COSI = MU_VALUES(J)
          D = CDSQRT(EPSILON - 1.0D0 + COSI**2)
          RH = (COSI - D) / (COSI + D)
          RV =(EPSILON*COSI - D) / (EPSILON*COSI + D)
          R1 = (CDABS(RV)**2 + CDABS(RH)**2 )/2.0D0
          R2 = (CDABS(RV)**2 - CDABS(RH)**2 )/2.0D0
          RADIANCE(1,J) = (1.0 - R1) * PLANCK
          IF (NSTOKES .GT. 1)  RADIANCE(2,J) = -R2 * PLANCK
        ENDDO
      ENDIF

      RETURN
      END





      SUBROUTINE THERMAL_RADIANCE (NSTOKES, NUMMU, MODE, 
     .                             TEMPERATURE, ALBEDO,
     .                             WAVELENGTH, RADIANCE)
C        THERMAL_RADIANCE returns a polarized radiance vector for
C      thermal emission at WAVELENGTH (microns) for a body with
C      ALBEDO and with a TEMPERATURE (Kelvins).  The radiance is in
C      the units:  Watts / (meter^2 ster micron).  The emission is
C      isotropic and unpolarized: only the first term in azimuth Fourier
C      series is nonzero, and all MU terms are the same.
      INTEGER  NSTOKES, NUMMU, MODE
      REAL*8   TEMPERATURE, WAVELENGTH, ALBEDO
      REAL*8   RADIANCE(NSTOKES,NUMMU,2)
      INTEGER  J, N
      REAL*8   PLANCK, THERMAL

      N = NSTOKES*NUMMU
      CALL MZERO (2*N, 1, RADIANCE)
      IF (MODE .EQ. 0) THEN
        CALL PLANCK_FUNCTION (TEMPERATURE, 'R', WAVELENGTH, PLANCK)
        THERMAL = (1.0-ALBEDO)*PLANCK
        DO J = 1, NUMMU
            RADIANCE(1,J,1) = THERMAL
            RADIANCE(1,J,2) = THERMAL
        ENDDO
      ENDIF
      RETURN
      END




      SUBROUTINE PLANCK_FUNCTION (TEMP, UNITS, WAVELENGTH, PLANCK)
C        Calculates the Planck blackbody radiance in
C      [Watts /(meter^2 ster micron)] for a temperature in [Kelvins]
C      at a wavelength in [microns].  If using temperature units then
C      the Planck function is simple the temperature.
      REAL*8  TEMP, WAVELENGTH, PLANCK
      CHARACTER*1  UNITS

      IF (UNITS .EQ. 'T') THEN
          PLANCK = TEMP
      ELSE
          IF (TEMP .GT. 0.0) THEN
              PLANCK = 1.1911D8 / WAVELENGTH**5
     .            / (DEXP(1.4388D4/(WAVELENGTH*TEMP)) - 1)
          ELSE
              PLANCK = 0.0
          ENDIF
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
      DO J = 1, NUM
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
          ENDDO
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        ABSCISSAS(NUM+1-J) = X
        WEIGHTS(NUM+1-J) = 2/((1-X*X)*DPL*DPL)
      ENDDO

      RETURN
      END


      SUBROUTINE DOUBLE_GAUSS_QUADRATURE
     .                          (NUM, ABSCISSAS, WEIGHTS)
C        Generates the abscissas and weights for an even 2*NUM point
C      Gauss-Legendre quadrature.  Only the NUM positive points are returned.
      INTEGER  NUM
      REAL*8   ABSCISSAS(*), WEIGHTS(*)
      INTEGER  N, K, I, J, L
      REAL*8   X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-14)

      N = NUM
      K = (N+1)/2
      DO J = 1, K
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
          ENDDO
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        ABSCISSAS(J) = (1D0-X)/2D0
        ABSCISSAS(NUM+1-J) = (1D0+X)/2D0
        WEIGHTS(NUM+1-J) = 1.0D0/((1.0D0-X*X)*DPL*DPL)
        WEIGHTS(J)       = 1.0D0/((1.0D0-X*X)*DPL*DPL)
      ENDDO

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
      DO J = 1, NUM-1
        X = SIN(3.141592654*(J-CI)/(N-.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO L = 2, N1
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
          ENDDO
          DPL = N1*(X*PL-PL1)/(X*X-1)
          D2PL = (2.D0*X*DPL-N1*(N1+1)*PL) / (1D0-X*X)
          XP = X
          X = XP - DPL/D2PL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        ABSCISSAS(J) = X
        WEIGHTS(J) = 2.0D0/(N*N1*PL*PL)
      ENDDO
      ABSCISSAS(NUM) = 1.D0
      WEIGHTS(NUM) = 2.D0/(N*N1)

      RETURN
      END





