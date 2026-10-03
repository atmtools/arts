C      SUBROUTINE RADTRAN (NSTOKES, NUMMU, AZIORDER, MAX_DELTA_TAU,
C     .                    SRC_CODE, QUAD_TYPE, DELTAM,
C     .                    DIRECT_FLUX, DIRECT_MU,
C     .                    GROUND_TEMP, GROUND_TYPE,
C     .                    GROUND_ALBEDO, GROUND_INDEX,
C     .                    SKY_TEMP, WAVELENGTH,
C     .                    NUM_LAYERS, HEIGHT, TEMPERATURES,
C     .                    GAS_EXTINCT, SCAT_FILES,
C     .                    NOUTLEVELS, OUTLEVELS,
C     .                    MU_VALUES, UP_FLUX, DOWN_FLUX,
C     .                    UP_RAD, DOWN_RAD)
C
C        RADTRAN solves the plane-parallel polarized radiative transfer
C    equation for an inhomogenous atmosphere with randomly oriented
C    particles.  The sources of radiation are solar direct beam and thermal
C    emission.  The ground surface has Lambertian or Fresnel reflection and
C    emission.
C        Input is the relevant parameters for the atmospheric layers and
C    the boundary conditions for a number of different cases.  Output is
C    a Fourier azimuth series of the upwelling and downwelling radiances 
C    at selected levels at the discrete angles.
C        The heights and temperatures are specified at the layer interfaces,
C    i.e. N+1 values for N layers.  The gaseous extinction and scattering
C    files are specified for the N layers.  The layers are listed from the
C    top to the bottom of the atmosphere: HEIGHT(1) is top of the highest
C    layer and HEIGHT(N+1) is bottom of the lowest layer, SCAT_FILES(1) is
C    the scattering file of the top layer and SCAT_FILES(N) is the scattering
C    file of the bottom layer.
C
C
C    Parameter         Type             Description
C
C  NSTOKES           INTEGER       Number of Stokes parameters: 1 for I
C                                    (no polarization), 2 for I,Q,
C                                    3 for I,Q,U,  4 for I,Q,U,V.
C  NUMMU             INTEGER       Number of quadrature angles
C                                    (per hemisphere).
C  AZIORDER          INTEGER       Order of Fourier azimuth series:
C                                    0 is azimuthially symmetric case.
C  MAX_DELTA_TAU     REAL          Initial layer thickness for doubling;
C                                    governs accuracy, 10E-5 should be
C                                    adequate.  Don't go beyond half the
C                                    real precision, i.e. 10e-8 for REAL*8.
C  SRC_CODE          INTEGER       Radiation sources included:
C                                    none=0, solar=1, thermal=2, both=3
C  QUAD_TYPE         CHAR*1        Type of quadrature used: 
C                                    (G-gaussian, D-double gaussian, 
C                                     L-Lobatto, E-extra-angle).
C                                    If extra-angles then end of
C                                    MU_VALUES(<=NUMMU) contains the extra
C                                    angles and rest is zero.
C  DELTAM            CHAR*1        Delta-M scaling flag (Y or N)
C                                    If DELTAM='Y' then the scattering
C                                    properties are delta-M scaled when read in.
C  DIRECT_FLUX       REAL          Flux on horizontal plane from direct
C                                    (solar) source;  units W/(m*m)/um or K.
C  DIRECT_MU         REAL          Cosine of solar zenith angle
C
C  GROUND_TEMP       REAL          Ground surface temperature in Kelvin
C  GROUND_TYPE       CHAR*1        Type of ground surface:
C                                    L for Lambertian, F for Fresnel.
C                                    Only Lambertian allowed for solar source.
C  GROUND_ALBEDO     REAL          Albedo of Lambertian surface
C  GROUND_INDEX      COMPLEX       Index of refraction of Fresnel surface
C  SKY_TEMP          REAL          Temperature of blackbody radiation
C                                    incident on atmosphere from above
C  WAVELENGTH        REAL          Wavelength of radiation in microns.
C
C  NUM_LAYERS        INTEGER       Number of atmosphere layers input
C  HEIGHT            REAL array    Height of layer interfaces from top down
C                                    Units are inverse of units of extinction
C                                    and scattering, e.g. km.
C  TEMPERATURES      REAL array    Temperature (Kelvins) of layer interfaces
C  GAS_EXTINCT       REAL array    Gaseous (nonscattering) extinction of layers
C                                    For processes not in scattering file
C  SCAT_FILES        CHAR*64 array Names of scattering files for layers
C                                    String format 'RAIN.SCA', for no
C                                    scattering use ' '.  See example for
C                                    format of scattering file.
C
C  NOUTLEVELS        INTEGER       Number of output levels
C  OUTLEVELS         INTEGER       The levels numbers to output at,
C                                    from 1 at top to NUM_LAYERS+1 at bottom.
C
C  MU_VALUES         REAL array    Output quadrature angle values
C                                    (also input for QUAD_TYPE='E')
C  UP_FLUX           REAL array    Upward flux for each Stokes parameter
C                                    at each output level 
C                                    UP_FLUX(NSTOKES,NOUTLEVELS)
C  DOWN_FLUX         REAL array    Downward flux (NSTOKES,NOUTLEVELS)
C  UP_RAD            REAL array    Upward radiances
C                                    (NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)
C  DOWN_RAD          REAL array    Downward radiances
C                                    (NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)
C



      SUBROUTINE RADTRAN (NSTOKES, NUMMU, AZIORDER,  MAX_DELTA_TAU,
     .                    SRC_CODE, QUAD_TYPE, DELTAM,
     .                    DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH,
     .                    NUM_LAYERS, HEIGHT, TEMPERATURES,
     .                    GAS_EXTINCT, SCAT_FILES,
     .                    NOUTLEVELS, OUTLEVELS,
     .                    MU_VALUES, UP_FLUX, DOWN_FLUX,
     .                    UP_RAD, DOWN_RAD)
      INTEGER   NSTOKES, NUMMU, AZIORDER
      INTEGER   NUM_LAYERS, SRC_CODE
      INTEGER   NOUTLEVELS, OUTLEVELS(*)
      REAL*8    DIRECT_FLUX, DIRECT_MU
      REAL*8    GROUND_TEMP, GROUND_ALBEDO
      COMPLEX*16  GROUND_INDEX
      REAL*8    SKY_TEMP
      REAL*8    WAVELENGTH, MAX_DELTA_TAU
      REAL*8    HEIGHT(*), TEMPERATURES(*), GAS_EXTINCT(*)
      REAL*8    MU_VALUES(*)
      REAL*8    UP_FLUX(*), DOWN_FLUX(*)
      REAL*8    UP_RAD(*), DOWN_RAD(*)
      CHARACTER*1  QUAD_TYPE, DELTAM, GROUND_TYPE
      CHARACTER*64 SCAT_FILES(*)


      INTEGER   MAXV, MAXM, MAXLM, MAXLEG, MAXLAY, MAXSBUF, MAXDBUF
      PARAMETER (MAXV=64, MAXM=4096, MAXLM=101*MAXM)
      PARAMETER (MAXLEG=1024, MAXLAY=200)
      PARAMETER (MAXSBUF=16*MAXLAY*2*MAXM , MAXDBUF=16*MAXLAY*2*MAXV)

      REAL*8    PI, TWOPI, ZERO
      PARAMETER (PI = 3.1415926535897932384D0, TWOPI=2.0D0*PI)
      PARAMETER (ZERO=0.0D0)

      INTEGER   NUMLEGEN, NLEGLIM, MODE, LAYER, NUM_DOUBLES
      INTEGER   SCAT_NUMS(MAXLAY), SCAT_NUM
      INTEGER   I, J, K, N, NA, KRT, KS, L, LI
      LOGICAL   SYMMETRIC
      REAL*8    EXTINCTIONS(MAXLAY), ALBEDOS(MAXLAY)
      REAL*8    EXPFACTOR, LINFACTOR
      REAL*8    PLANCK0, PLANCK1, TMP
      REAL*8    ALBEDO, EXTINCTION, EXTINCT, SCATTER
      REAL*8    ZDIFF, DELTA_Z, F, NUM_SUB_LAYERS, TAU
      REAL*8    QUAD_WEIGHTS(MAXV)
      REAL*8    LEGENDRE_COEF(6,MAXLEG)
      REAL*8    SCATBUF(MAXSBUF), DIRECTBUF(MAXDBUF)
      REAL*8    SCATTER_MATRIX(4*MAXM)
      REAL*8    DIRECT_LEVEL_FLUX(MAXLAY+1)
      REAL*8    DIRECT_VECTOR(2*MAXV), EXP_SOURCE(2*MAXV)
      REAL*8    THERMAL_VECTOR(2*MAXV), LIN_SOURCE(2*MAXV)
      REAL*8    REFLECT1(2*MAXM),UPREFLECT(2*MAXM),DOWNREFLECT(2*MAXM)
      REAL*8    TRANS1(2*MAXM),  UPTRANS(2*MAXM),  DOWNTRANS(2*MAXM)
      REAL*8    SOURCE1(2*MAXV), UPSOURCE(2*MAXV), DOWNSOURCE(2*MAXV)
      REAL*8    REFLECT(2*MAXLM)
      REAL*8    TRANS(2*MAXLM)
      REAL*8    SOURCE(2*MAXV*(MAXLAY+1))
      REAL*8    GND_RADIANCE(MAXV), SKY_RADIANCE(2*MAXV)
      CHARACTER*64 SCAT_FILE


      IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
        IF (GROUND_TYPE .NE. 'L') THEN
          WRITE (*,*) 'Solar case requires Lambertian surface'
          STOP
        ENDIF
      ENDIF

      SYMMETRIC = .FALSE.
      IF (NSTOKES .LE. 2 )  SYMMETRIC = .TRUE.

      N = NSTOKES*NUMMU
      IF (N .GT. MAXV) THEN
          WRITE (*,'(1X,A,I3)')
     .     'Vector size exceeded.  Maximum size :', MAXV
          STOP
      ELSE IF (N*N .GT. MAXM) THEN
          WRITE (*,'(1X,A,I3)')
     .     'Matrix size exceeded.  Maximum size :', MAXM
          STOP
      ENDIF
      IF ((SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) .AND.
     .    MAXDBUF .LT. (AZIORDER+1)*2*N*NUM_LAYERS) THEN
          WRITE (*,'(1X,A,I3)') 'Direct source buffer size exceeded.'
          STOP
      ENDIF
      IF (NUM_LAYERS .GT. MAXLAY) THEN
          WRITE (*,'(1X,A,I3)')
     .     'Too many layers.  Maximum number :', MAXLAY
          STOP
      ENDIF
      IF ((NUM_LAYERS+1)*N*N .GT. MAXLM) THEN
          WRITE (*,'(1X,A,A,I8)') 'Matrix layer size exceeded.',
     .     '  Maximum number :', MAXLM
          STOP
      ENDIF


C           Make the desired quadrature abscissas and weights
      IF (QUAD_TYPE(1:1) .EQ. 'D') THEN
          CALL DOUBLE_GAUSS_QUADRATURE
     .                       (NUMMU, MU_VALUES, QUAD_WEIGHTS)
          NLEGLIM = 2*NUMMU - 3
      ELSE IF (QUAD_TYPE(1:1) .EQ. 'L') THEN
          CALL LOBATTO_QUADRATURE
     .                       (NUMMU, MU_VALUES, QUAD_WEIGHTS)
          NLEGLIM = 4*NUMMU - 5
      ELSE IF (QUAD_TYPE(1:1) .EQ. 'E') THEN
          J = NUMMU
          DO I = NUMMU, 1, -1
              IF (MU_VALUES(I) .NE. 0.0) THEN
                  QUAD_WEIGHTS(I) = 0.0
                  J = I - 1
              ENDIF
          ENDDO
          CALL GAUSS_LEGENDRE_QUADRATURE
     .                       (J, MU_VALUES, QUAD_WEIGHTS)
          NLEGLIM = 4*J - 3
      ELSE
          CALL GAUSS_LEGENDRE_QUADRATURE
     .                       (NUMMU, MU_VALUES, QUAD_WEIGHTS)
          NLEGLIM = 4*NUMMU - 3
      ENDIF
      NLEGLIM = MAX(NLEGLIM,1)




C       Make all of the scattering matrices ahead of time
C         and store them in memory
      SCAT_NUM = 0
      SCAT_FILE = '&&&'
C           Loop through the layers
      DO LAYER = 1, NUM_LAYERS
C                   Special case for a non-scattering layer
          IF (SCAT_FILES(LAYER) .EQ. ' ') THEN
              SCAT_FILE = SCAT_FILES(LAYER)
              EXTINCT = 0.0
              SCATTER = 0.0
          ELSE
C                   If there is a new scattering file then do it
            IF (SCAT_FILES(LAYER) .NE. SCAT_FILE)  THEN
              SCAT_FILE = SCAT_FILES(LAYER)
              SCAT_NUM = SCAT_NUM + 1
              IF (SCAT_NUM*(AZIORDER+1)*2*(NUMMU*NSTOKES)**2 
     .              .GT. MAXSBUF) THEN
                WRITE (*,'(1X,A,I3)') 
     .            'Scattering matrix buffer size exceeded.'
                STOP
              ENDIF
C                   Read scattering file in
              CALL READ_SCAT_FILE (SCAT_FILE, DELTAM, NUMMU, NUMLEGEN,
     .                             LEGENDRE_COEF, EXTINCT, SCATTER)
              IF (NUMLEGEN .GT. MAXLEG) THEN
                  WRITE (*,*)  'Too many Legendre terms.'
                  STOP
              ENDIF
C                   Truncate the Legendre series to enforce normalization
              IF (NUMLEGEN .GT. NLEGLIM) THEN
                  WRITE (*,*) 'Truncating Legendre series for file:',
     .                         SCAT_FILE, NUMLEGEN, NLEGLIM
                  NUMLEGEN = NLEGLIM
              ENDIF
C                   Make the scattering matrix
              CALL SCATTERING (NUMMU, AZIORDER, NSTOKES,
     .                     MU_VALUES, QUAD_WEIGHTS,
     .                     NUMLEGEN, LEGENDRE_COEF,
     .                     SCAT_NUM, SCATBUF)
C                   Make the direct (solar) pseudo source
              IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
                  CALL DIRECT_SCATTERING (NUMMU, AZIORDER, NSTOKES,
     .                    MU_VALUES, NUMLEGEN, LEGENDRE_COEF,
     .                    DIRECT_MU, SCAT_NUM, DIRECTBUF)
              ENDIF
            ENDIF
          ENDIF

          SCAT_NUMS(LAYER) = SCAT_NUM
          EXTINCTIONS(LAYER) = EXTINCT + MAX(GAS_EXTINCT(LAYER),0.0D0)
          IF (EXTINCTIONS(LAYER) .GT. 0.0) THEN
            ALBEDOS(LAYER) = SCATTER/EXTINCTIONS(LAYER)
          ELSE
            ALBEDOS(LAYER) = 0.0
          ENDIF
      ENDDO


C           Compute the direct beam flux at each level
      IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
        TAU = 0.0
        DIRECT_LEVEL_FLUX(1) = DIRECT_FLUX
        DO LAYER = 1, NUM_LAYERS
          TAU = TAU + EXTINCTIONS(LAYER)/DIRECT_MU*
     .                ABS(HEIGHT(LAYER)-HEIGHT(LAYER+1))
          DIRECT_LEVEL_FLUX(LAYER+1) = DIRECT_FLUX*DEXP(-TAU)
        ENDDO
      ENDIF




C       Loop through each azimuth mode
      DO MODE = 0, AZIORDER
        SCAT_NUM = 0
C     ------------------------------------------------------
C           Loop through the layers

        DO LAYER = 1, NUM_LAYERS
C                   Calculate the layer thickness
          ZDIFF = ABS(HEIGHT(LAYER) - HEIGHT(LAYER+1))
          EXTINCTION = EXTINCTIONS(LAYER)
          ALBEDO = ALBEDOS(LAYER)

          IF (SCAT_NUMS(LAYER) .NE. SCAT_NUM) THEN
            SCAT_NUM = SCAT_NUMS(LAYER)
C                   Get the scattering matrix from the buffer
            CALL GET_SCATTERING (NSTOKES, NUMMU, MODE, AZIORDER,
     .                           SCAT_NUM, SCATBUF, SCATTER_MATRIX)
C                   Check the normalization of the scattering matrix
            IF (MODE .EQ. 0) THEN
                CALL CHECK_NORM (NSTOKES, NUMMU, QUAD_WEIGHTS,
     .                           SCATTER_MATRIX)
            ENDIF
C                   Get the direct (solar) vector from the buffer
            IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
                CALL GET_DIRECT (NSTOKES, NUMMU, MODE, AZIORDER,
     .                           SCAT_NUM, DIRECTBUF, DIRECT_VECTOR)
            ENDIF
          ENDIF


C                   Compute the thermal emission at top and bottom of layer
          IF (SRC_CODE .EQ. 2 .OR. SRC_CODE .EQ. 3) THEN
C                   Calculate the thermal source for end of layer
              CALL THERMAL_RADIANCE (NSTOKES, NUMMU, MODE,
     .                       TEMPERATURES(LAYER+1), ALBEDO,
     .                       WAVELENGTH, THERMAL_VECTOR)
              PLANCK1 = THERMAL_VECTOR(1)
C                   Calculate the thermal source for beginning of layer
              CALL THERMAL_RADIANCE (NSTOKES, NUMMU, MODE,
     .                       TEMPERATURES(LAYER), ALBEDO,
     .                       WAVELENGTH, THERMAL_VECTOR)
              PLANCK0 = THERMAL_VECTOR(1)
          ELSE
              PLANCK0 = 0.0
              PLANCK1 = 0.0
          ENDIF

 
          KRT = 1 + 2*N*N*(LAYER-1)
          KS = 1 + 2*N*(LAYER-1)
          IF (ALBEDO .EQ. ZERO) THEN
C                   If the layer is purely absorbing then quickly
C                     make the reflection and transmission matrices
C                     and source vector instead of doubling.
              CALL NONSCATTER_LAYER (NSTOKES, NUMMU, MODE,
     .                           ZDIFF*EXTINCTION, MU_VALUES,
     .                           PLANCK0, PLANCK1,
     .                           REFLECT(KRT), TRANS(KRT), SOURCE(KS))
          ELSE

C                   Find initial thickness of sublayer and
C                     the number of times to double
              F = MAX(EXTINCTION*ZDIFF,1.0D-7)
              F = LOG(F/MAX_DELTA_TAU) /LOG(2.)
              NUM_DOUBLES = 0
              IF (F .GT. 0.0)  NUM_DOUBLES = INT(F) + 1
              NUM_SUB_LAYERS = 2.0**NUM_DOUBLES
              DELTA_Z = ZDIFF / NUM_SUB_LAYERS

C                   For a solar source make the pseudo source vector
C                     and initialize it
              IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
                  TMP = DIRECT_LEVEL_FLUX(LAYER) *ALBEDO 
     .                      /(4.0D0*PI*DIRECT_MU)
                  CALL MSCALARMULT (2*N, 1, TMP, DIRECT_VECTOR,SOURCE1)
                  CALL INITIAL_SOURCE (NSTOKES, NUMMU, N,
     .                      DELTA_Z, MU_VALUES, EXTINCTION,
     .                      SOURCE1,  EXP_SOURCE)
                  EXPFACTOR = DEXP(-EXTINCTION*DELTA_Z/DIRECT_MU)
              ENDIF

C                   Initialize the thermal source vector
              IF (SRC_CODE .EQ. 2 .OR. SRC_CODE .EQ. 3) THEN
                  CALL INITIAL_SOURCE (NSTOKES, NUMMU, N,
     .                        DELTA_Z, MU_VALUES, EXTINCTION,
     .                        THERMAL_VECTOR,  LIN_SOURCE)
                  IF (PLANCK0 .EQ. 0.0) THEN
                      LINFACTOR = 0.0
                  ELSE
                      LINFACTOR = (PLANCK1/PLANCK0-1.0)/NUM_SUB_LAYERS
                  ENDIF
              ENDIF

C                Generate the local reflection and transmission matrices
              CALL INITIALIZE (NSTOKES, NUMMU, N, DELTA_Z, MU_VALUES, 
     .                         EXTINCTION, ALBEDO, SCATTER_MATRIX, 
     .                         REFLECT1, TRANS1)

C                   Double up to the thickness of the layer
              CALL DOUBLING_INTEGRATION (N, NUM_DOUBLES, 
     .               SRC_CODE, SYMMETRIC, REFLECT1, TRANS1,
     .               EXP_SOURCE, EXPFACTOR, LIN_SOURCE, LINFACTOR,
     .               REFLECT(KRT), TRANS(KRT), SOURCE(KS))

          ENDIF

        ENDDO
C            End of layer loop


C           Get the surface reflection and transmission matrices
C             and the surface radiance
        KRT = 1 + 2*N*N*(NUM_LAYERS)
        KS = 1 + 2*N*(NUM_LAYERS)
        IF (GROUND_TYPE .EQ. 'F') THEN
C               For a Fresnel surface
          CALL FRESNEL_SURFACE (NSTOKES, NUMMU, 
     .                        MU_VALUES, GROUND_INDEX, 
     .                        REFLECT(KRT), TRANS(KRT), SOURCE(KS))
C                The radiance from the ground is thermal
          CALL FRESNEL_RADIANCE (NSTOKES, NUMMU, MODE,
     .                  MU_VALUES, GROUND_INDEX, GROUND_TEMP,
     .                  WAVELENGTH, GND_RADIANCE)
        ELSE
C               For a Lambertian surface
          CALL LAMBERT_SURFACE (NSTOKES, NUMMU, MODE,
     .                       MU_VALUES, QUAD_WEIGHTS, GROUND_ALBEDO,
     .                       REFLECT(KRT), TRANS(KRT), SOURCE(KS))
C                The radiance from the ground is thermal and reflected direct
          CALL LAMBERT_RADIANCE (NSTOKES, NUMMU, MODE,
     .                   SRC_CODE, GROUND_ALBEDO, GROUND_TEMP,
     .                   WAVELENGTH, DIRECT_LEVEL_FLUX(NUM_LAYERS+1),
     .                   GND_RADIANCE)
        ENDIF

C           Assume the radiation coming from above is blackbody radiation
        CALL THERMAL_RADIANCE (NSTOKES, NUMMU, MODE, SKY_TEMP, ZERO,  
     .                         WAVELENGTH,  SKY_RADIANCE)


 
C               For each desired output level (1 thru NL+2) add layers
C               above and below level and compute internal radiance
        DO I = 1, NOUTLEVELS
            LAYER = MIN( MAX( OUTLEVELS(I), 1), NUM_LAYERS+2)
	    CALL MZERO (2*N, N, UPREFLECT)
	    CALL MZERO (2*N, N, DOWNREFLECT)
            CALL MIDENTITY (N, UPTRANS(1))
            CALL MIDENTITY (N, UPTRANS(1+N*N))
            CALL MIDENTITY (N, DOWNTRANS(1))
            CALL MIDENTITY (N, DOWNTRANS(1+N*N))
	    CALL MZERO (2*N, 1, UPSOURCE)
	    CALL MZERO (2*N, 1, DOWNSOURCE)
            DO L = 1, LAYER-1
                KRT = 1 + 2*N*N*(L-1)
                KS = 1 + 2*N*(L-1)
                IF (L .EQ. 1) THEN
                  CALL MCOPY (2*N,N, REFLECT(KRT), UPREFLECT)
                  CALL MCOPY (2*N,N, TRANS(KRT), UPTRANS)
                  CALL MCOPY (2*N,1, SOURCE(KS), UPSOURCE)
                ELSE
                  CALL MCOPY (2*N,N, UPREFLECT, REFLECT1)
                  CALL MCOPY (2*N,N, UPTRANS, TRANS1)
                  CALL MCOPY (2*N,1, UPSOURCE, SOURCE1)
                  CALL COMBINE_LAYERS (N, REFLECT1, TRANS1, SOURCE1,
     .                        REFLECT(KRT), TRANS(KRT), SOURCE(KS),
     .                        UPREFLECT, UPTRANS, UPSOURCE)
                ENDIF
            ENDDO
            DO L = LAYER, NUM_LAYERS+1
                KRT = 1 + 2*N*N*(L-1)
                KS = 1 + 2*N*(L-1)
                IF (L .EQ. LAYER) THEN
                  CALL MCOPY (2*N,N, REFLECT(KRT), DOWNREFLECT)
                  CALL MCOPY (2*N,N, TRANS(KRT), DOWNTRANS)
                  CALL MCOPY (2*N,1, SOURCE(KS), DOWNSOURCE)
                ELSE
                  CALL MCOPY (2*N,N, DOWNREFLECT, REFLECT1)
                  CALL MCOPY (2*N,N, DOWNTRANS, TRANS1)
                  CALL MCOPY (2*N,1, DOWNSOURCE, SOURCE1)
                  CALL COMBINE_LAYERS (N, REFLECT1, TRANS1, SOURCE1,
     .                        REFLECT(KRT), TRANS(KRT), SOURCE(KS),
     .                        DOWNREFLECT, DOWNTRANS, DOWNSOURCE)
                ENDIF
            ENDDO
            NA = N*(AZIORDER+1)
	    CALL INTERNAL_RADIANCE (N, UPREFLECT, UPTRANS, UPSOURCE,
     .                           DOWNREFLECT, DOWNTRANS, DOWNSOURCE,
     .                           SKY_RADIANCE, GND_RADIANCE,
     .                           UP_RAD(1+MODE*N+(I-1)*NA),
     .                           DOWN_RAD(1+MODE*N+(I-1)*NA))
        ENDDO
 
 
      ENDDO 
C         End of azimuth mode loop
 


C           Integrate mu times the radiance to find the fluxes
      DO L = 1, NOUTLEVELS
        NA = (AZIORDER+1) *NSTOKES*NUMMU
        DO I = 1, NSTOKES
          K = I+NSTOKES*(L-1)
          UP_FLUX(K) = 0.0
          DOWN_FLUX(K) = 0.0
          DO J = 1, NUMMU
            UP_FLUX(K) = UP_FLUX(K)
     .               + TWOPI*QUAD_WEIGHTS(J) * MU_VALUES(J)
     .               * UP_RAD(I+NSTOKES*(J-1)+NA*(L-1))
            DOWN_FLUX(K) = DOWN_FLUX(K)
     .               + TWOPI*QUAD_WEIGHTS(J) * MU_VALUES(J)
     .               * DOWN_RAD(I+NSTOKES*(J-1)+NA*(L-1))
          ENDDO
        ENDDO
C           Add in direct beam fluxes
        IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
          K = NSTOKES*(L-1)
          LI = OUTLEVELS(L)
          DOWN_FLUX(K+1) = DOWN_FLUX(K+1) + DIRECT_LEVEL_FLUX(LI)
        ENDIF
      ENDDO

      RETURN
      END


