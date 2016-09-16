C        RADTRANO solves the plane-parallel polarized radiative transfer
C    equation for an inhomogenous atmosphere with particles oriented 
C    in zenith angle but randomly oriented in azimuth (azimuthal symmetry).
C    Thermal emission is the only source of radiation.  The radiance is
C    computed in units of Watts/meter^2 micron ster.  The ground surface 
C    has Lambertian or Fresnel reflection and emission.
C        Input is the relevant parameters for the atmospheric layers and
C    the boundary conditions.  Output is the upward and downward fluxes
C    and radiances at the discrete quadrature angles for the specified levels.
C        The heights and temperatures are specified at the layer interfaces,
C    i.e. N+1 values for N layers.  The gaseous extinction and scattering
C    files are specified for the N layers.  The layers are listed from the
C    top to the bottom of the atmosphere: HEIGHT(1) is top of the highest
C    layer and HEIGHT(N+1) is bottom of the lowest layer, SCAT_FILES(1) is
C    the scattering file of the top layer and SCAT_FILES(N) is the scattering
C    file of the bottom layer.  The format of the oriented scattering files
C    is given below.
C
C
C    Parameter         Type             Description
C
C  NSTOKES           INTEGER       Number of Stokes parameters: 1 for I
C                                    (no polarization), 2 for I,Q.
C  NUMMU             INTEGER       Number of quadrature angles
C                                    (per hemisphere).
C  MAX_DELTA_TAU     REAL          Initial layer thickness for doubling;
C                                    governs accuracy, 10E-5 should be
C                                    adequate.  Don't go beyond half the
C                                    real precision, i.e. 10e-8 for REAL*8.
C  QUAD_TYPE         CHAR*1        Type of quadrature used: 
C                                    (G-gaussian, D-double gaussian, 
C                                    L-Lobatto, E-extra-angle).
C                                    If extra-angles then end of
C                                    MU_VALUES(<=NUMMU) contains the extra
C                                    angles and rest is zero.
C  GROUND_TEMP       REAL          Ground surface temperature in Kelvin
C  GROUND_TYPE       CHAR*1        Type of ground surface:
C                                    L for Lambertian, F for Fresnel.
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
C  SCAT_FILES        CHAR*64 array Names of oriented scattering files for layers
C                                    String format 'PLATE.DDA', for no
C                                    scattering use ' '.  See below for
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
C                                    (NSTOKES,NUMMU,NOUTLEVELS)
C  DOWN_RAD          REAL array    Downward radiances
C                                    (NSTOKES,NUMMU,NOUTLEVELS)
C
C
C             Format of Scattering Files
C
C  Nmu  Maz  'quadtype'   } Nmu is number of quadrature angles in each
C                           hemisphere.  Maz is azimuth series order.
C                           'quadtype' is the quadrature type (e.g. 'GAUSSIAN').
C  C  SCATTERING MATRIX
C                        _   (2*Nmu)^2 (2*Maz+1) sets of scattering matrices
C  mu'  mu   m  c/s       |
C  P11  P12  P13  P14     |  m is azimuth mode: 0, 1, 2, ...
C  P21  P22  P23  P24     |    c/s is C for cosine term, S for sine term
C  P31  P32  P33  P34     |    azimuth ordering is: 0, 1 C, 1 S, 2 C, 2 S, ...
C  P41  P42  P43  P44    _|  mu is outgoing cosine of zenith angle
C                            mu' is incoming cosine of zenith angle
C                            Order of indexing is:  m (fastest), mu, mu'
C  C  EXTINCTION MATRIX
C                        _   2Nmu sets of extinction matrices
C  mu'                    |
C  K11  K12  K13  K14     |  mu' is incoming cosine of zenith angle
C  K21  K22  K23  K24     |
C  K31  K32  K33  K34     |
C  K41  K42  K43  K44    _|
C
C  C  EMISSION VECTOR
C
C  mu'  S1  S2  S3  S4    }  2Nmu sets of emission vectors
C
C  Further description:
C      All angles (mu's) are stored by hemisphere: first from 0 to 1, then 
C  from 0 to -1.  There are no blank lines in the file.  There may be any
C  number of comment lines (beginning with "C") before the first data line,
C  but only the three comment lines shown (scattering, extinction, emission)
C  beyond that.  There are spaces between all elements in the format, so that
C  Fortran free format reading will work.
C      The (I,Q,U,V) Stokes basis is used for the polarization.  The frame
C  of reference for the polarization is meridional plane, i.e. the
C  vertical polarization vector is in the plane defined by the ray
C  direction and the vertical (z) axis.  Because of the azimuthal symmetry
C  of the medium, the scattering matrix depends only on the incident and
C  outgoing zenith angles and the difference between the incident and
C  outgoing azimuth angles (delphi).  The delphi dependence is represented
C  as a Fourier series, because the radiative transfer of azimuthal modes
C  is decoupled.  For this code only the m=0 (lowest, symmetric) mode is used.
C  The units of the scattering matrix, extinction matrix, and emission vector
C  are units of inverse length;  K11 = S1 + Integ{P11}, where Integ{} is 
C  integration over azimuth and zenith angles.




      SUBROUTINE RADTRANO (NSTOKES, NUMMU, NUUMMU, MAX_DELTA_TAU,
     .               QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,
     .               GROUND_ALBEDO, GROUND_INDEX,
     .               SKY_TEMP, WAVELENGTH,
     .               NUM_LAYERS, HEIGHT, TEMPERATURES,
c     .               GAS_EXTINCT, SCAT_FILES,
     .               GAS_EXTINCT,
     .               NSL, SCATLAYERS,
     .               EXTINCT_MATRIX, EMIS_VECTOR, SCATTER_MATRIX,
c     .               NOUTLEVELS, OUTLEVELS,
c     .               MU_VALUES, UP_FLUX, DOWN_FLUX,
     .               MU_VALUES,
     .               UP_RAD, DOWN_RAD)
      INTEGER   NSTOKES, NUMMU, NUM_LAYERS, NSL
c      INTEGER   NOUTLEVELS, OUTLEVELS(NOUTLEVELS)
      REAL*8    GROUND_TEMP, GROUND_ALBEDO
      COMPLEX*16  GROUND_INDEX
      REAL*8    SKY_TEMP
      REAL*8    WAVELENGTH, MAX_DELTA_TAU
      REAL*8    HEIGHT(NUM_LAYERS+1), TEMPERATURES(NUM_LAYERS+1)
      REAL*8    GAS_EXTINCT(NUM_LAYERS)
      REAL*8    MU_VALUES(NUMMU)
c      REAL*8    UP_FLUX(NSTOKES*NUM_LAYERS+1)
c      REAL*8    DOWN_FLUX(NSTOKES*NUM_LAYERS+1)
      REAL*8    UP_RAD(NSTOKES*NUMMU*NUM_LAYERS+1)
      REAL*8    DOWN_RAD(NSTOKES*NUMMU*NUM_LAYERS+1)
      CHARACTER*1  QUAD_TYPE, GROUND_TYPE
c      CHARACTER*64 SCAT_FILES(*)
      REAL*8    SCATLAYERS(NUM_LAYERS)
      REAL*8    EXTINCT_MATRIX(NSTOKES,NSTOKES,NUMMU,2,NSL)
      REAL*8    EMIS_VECTOR(NSTOKES,NUMMU,2,NSL)
      REAL*8    SCATTER_MATRIX(NSTOKES,NUMMU,NSTOKES,NUMMU,4,NSL)

      INTEGER   MAXV, MAXM, MAXLAY, MAXLM
      PARAMETER (MAXV=64, MAXM=4096, MAXLAY=200, MAXLM=201*4096)

      REAL*8    PI, TWOPI, ZERO
      PARAMETER (PI = 3.1415926535897932384D0, TWOPI=2.0D0*PI)
      PARAMETER (ZERO=0.0D0)

      INTEGER   LAYER, NUM_DOUBLES
      INTEGER   I, J, K, L, N, KRT, KS, TSL
      LOGICAL   SYMMETRIC
      REAL*8    LINFACTOR
      REAL*8    PLANCK0, PLANCK1
      REAL*8    ZDIFF, DELTA_Z, F, NUM_SUB_LAYERS, EXTINCT
      REAL*8    QUAD_WEIGHTS(MAXV)
c      REAL*8    SCATTER_MATRIX(4*MAXM)
      REAL*8    LIN_SOURCE(2*MAXV)
c      REAL*8    EXTINCT_MATRIX(16*2*MAXV), EMIS_VECTOR(4*2*MAXV)
      REAL*8    REFLECT1(2*MAXM),UPREFLECT(2*MAXM),DOWNREFLECT(2*MAXM)
      REAL*8    TRANS1(2*MAXM),  UPTRANS(2*MAXM),  DOWNTRANS(2*MAXM)
      REAL*8    SOURCE1(2*MAXV), UPSOURCE(2*MAXV), DOWNSOURCE(2*MAXV)
      REAL*8    REFLECT(2*MAXLM)
      REAL*8    TRANS(2*MAXLM)
      REAL*8    SOURCE(2*MAXV*(MAXLAY+1))
      REAL*8    GND_RADIANCE(MAXV), SKY_RADIANCE(2*MAXV)
c      CHARACTER*64 SCAT_FILE


c     this is dangerous to do, as we use nstokes also as array-size
c     determining parameter. resetting it will cause problems later on
c     shaping the arrays and correct value extraction. so, don't do
c     this. just check whether (and fail if) condition is not fulfilled.
c      NSTOKES = MIN(NSTOKES,2)
      IF (NSTOKES .GT. 2) THEN
          WRITE (*,'(1X,A,I3)')
     .     'Number of Stokes parameters exceeded.  Maximum size :', 2
          STOP
      ENDIF

      SYMMETRIC = .TRUE.
      N = NSTOKES*NUMMU
      IF (N .GT. MAXV) THEN
          WRITE (*,'(1X,2(A,I4))')
     .     'Vector size exceeded.  Maximum size :', MAXV,
     .     '.  Yours is ', N
          STOP
      ENDIF
      IF (N*N .GT. MAXM) THEN
          WRITE (*,'(1X,A,I5,2(A,I4),A,I5)')
     .     'Matrix size exceeded.  Maximum size :', MAXM,
     .     '.  Yours is ', N, '*', N, ' = ', N*N
          STOP
      ENDIF
      IF (NUM_LAYERS .GT. MAXLAY) THEN
          WRITE (*,'(1X,A,3(A,I4))') 'Number of layers exceeded.',
     .     '  Maximum number :', MAXLAY,
     .     '.  Yours is ', NUM_LAYERS
          STOP
      ENDIF
      IF ((NUM_LAYERS+1)*N*N .GT. MAXLM) THEN
          WRITE (*,'(1X,A,A,I5,3(A,I3),A,I5)')
     .     'Matrix layer size exceeded.',
     .     '  Maximum number (num_layers+1)*(nstokes*nummu)^2:', MAXLM,
     .     '.  Yours is (', NUM_LAYERS, '+1)*(', NSTOKES, '*', NUMMU,
     .     ')^2 = ', (NUM_LAYERS+1)*N*N
          STOP
      ENDIF


C           Make the desired quadrature abscissas and weights
c      WRITE(*,'(3A)') '>',QUAD_TYPE,'<'
      IF ( NUUMMU .GT. 0 ) THEN
        DO I = NUMMU, NUUMMU, -1
          QUAD_WEIGHTS(I) = 0.0
        ENDDO
      ENDIF

      J = NUMMU-NUUMMU
      IF ( QUAD_TYPE(1:1) .EQ. 'D' ) THEN
c        WRITE(*,'(A)') 'double gauss'
        CALL DOUBLE_GAUSS_QUADRATURE(J, MU_VALUES, QUAD_WEIGHTS)
      ELSE IF ( QUAD_TYPE(1:1) .EQ. 'L' ) THEN
c        WRITE(*,'(A)') 'lobatto'
        CALL LOBATTO_QUADRATURE(J, MU_VALUES, QUAD_WEIGHTS)
      ELSE
c        WRITE(*,'(A)') 'simple gauss'
        CALL GAUSS_LEGENDRE_QUADRATURE(J, MU_VALUES, QUAD_WEIGHTS)
      ENDIF




C     ------------------------------------------------------
C           Loop through the layers
C              Do doubling to make the reflection and transmission matrices
C              and soure vectors for each layer, which are stored.

      IF (NSL .GT. 0) THEN
          CALL CHECK_NORM (NSTOKES, NUMMU, NSL,
     .                     QUAD_WEIGHTS,
     .                     SCATTER_MATRIX,
     .                     EXTINCT_MATRIX, EMIS_VECTOR)
      END IF
c      WRITE(*,'(A)') 'norm check done'

      DO LAYER = 1, NUM_LAYERS
c          WRITE(*,'(A,I4)') 'processing layer', LAYER
C                   Calculate the layer thickness
          ZDIFF = ABS(HEIGHT(LAYER) - HEIGHT(LAYER+1))
          GAS_EXTINCT(LAYER) = MAX(GAS_EXTINCT(LAYER),0.0D0)

C                   Do the stuff for thermal source in layer
C                   Calculate the thermal source for end of layer
          CALL PLANCK_FUNCTION (TEMPERATURES(LAYER+1), 'R',
     .                          WAVELENGTH, PLANCK1)
C                   Calculate the thermal source for beginning of layer
          CALL PLANCK_FUNCTION (TEMPERATURES(LAYER), 'R',
     .                          WAVELENGTH, PLANCK0)

          KRT = 1 + 2*N*N*(LAYER-1)
          KS = 1 + 2*N*(LAYER-1)

          TSL = NINT(SCATLAYERS(LAYER))
          IF (TSL .LT. 1) THEN
c              WRITE(*,'(A)') 'non-scatt layer'
C                   If the layer is purely absorbing then quickly
C                     make the reflection and transmission matrices
C                     and source vector instead of doubling.
              CALL NONSCATTER_LAYER (NSTOKES, NUMMU, 
     .                          ZDIFF*GAS_EXTINCT(LAYER), MU_VALUES,
     .                          PLANCK0, PLANCK1,
     .                          REFLECT(KRT), TRANS(KRT), SOURCE(KS))
          ELSE

C                   Find initial thickness of sublayer and
C                     the number of times to double
c              WRITE(*,'(A)') 'scatt layer'
              EXTINCT = EXTINCT_MATRIX(1,1,1,1,TSL)+GAS_EXTINCT(LAYER)
              F =DLOG(MAX(EXTINCT*ZDIFF,1.0D-7)/MAX_DELTA_TAU)/LOG(2.)
              NUM_DOUBLES = 0
              IF (F .GT. 0.0)  NUM_DOUBLES = INT(F) + 1
              NUM_SUB_LAYERS = 2.0**NUM_DOUBLES
              DELTA_Z = ZDIFF / NUM_SUB_LAYERS

C                   Initialize the source vector
              CALL INITIAL_SOURCE (NSTOKES, NUMMU,
     .                    DELTA_Z, MU_VALUES,
     .                    PLANCK0, EMIS_VECTOR(1,1,1,TSL),
     .                    GAS_EXTINCT(LAYER),
     .                    LIN_SOURCE)
              IF (PLANCK0 .EQ. 0.0) THEN
                  LINFACTOR = 0.0
              ELSE
                  LINFACTOR = (PLANCK1/PLANCK0-1.0D0) /NUM_SUB_LAYERS
              ENDIF

C                Generate the local reflection and transmission matrices
              CALL INITIALIZE (NSTOKES, NUMMU,
     .                         DELTA_Z, MU_VALUES, QUAD_WEIGHTS,
     .                         GAS_EXTINCT(LAYER),
     .                         EXTINCT_MATRIX(1,1,1,1,TSL),
     .                         SCATTER_MATRIX(1,1,1,1,1,TSL),
     .                         REFLECT1, TRANS1)

C                   Double up to the thickness of the layer
              CALL DOUBLING_INTEGRATION (N, NUM_DOUBLES, SYMMETRIC,
     .                     REFLECT1, TRANS1, LIN_SOURCE, LINFACTOR,
     .                     REFLECT(KRT), TRANS(KRT), SOURCE(KS))
          ENDIF

      ENDDO
C            End of layer loop

c      WRITE(*,'(A)') 'layer looping finished'

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
        CALL FRESNEL_RADIANCE (NSTOKES, NUMMU,
     .                  MU_VALUES, GROUND_INDEX, GROUND_TEMP,
     .                  WAVELENGTH, GND_RADIANCE)
      ELSE
C               For a Lambertian surface
        CALL LAMBERT_SURFACE (NSTOKES, NUMMU, 0,
     .                       MU_VALUES, QUAD_WEIGHTS, GROUND_ALBEDO,
     .                       REFLECT(KRT), TRANS(KRT), SOURCE(KS))
C                The radiance from the ground is thermal and reflected direct
        CALL LAMBERT_RADIANCE (NSTOKES, NUMMU, 
     .         GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH, GND_RADIANCE)
      ENDIF

C           Assume the radiation coming from above is blackbody radiation
      CALL THERMAL_RADIANCE (NSTOKES, NUMMU, SKY_TEMP, ZERO,  
     .                       WAVELENGTH,  SKY_RADIANCE)


C         For each desired output level (1 thru NL+2) add layers 
C           above and below level and compute internal radiance.
C           OUTLEVELS gives the desired output levels.
      DO I = 1, NUM_LAYERS+1
        LAYER = MIN( MAX( I, 1), NUM_LAYERS+2)
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
        CALL INTERNAL_RADIANCE (N, UPREFLECT, UPTRANS, UPSOURCE,
     .                          DOWNREFLECT, DOWNTRANS, DOWNSOURCE, 
     .                          SKY_RADIANCE, GND_RADIANCE,
     .                          UP_RAD(1+(I-1)*N), DOWN_RAD(1+(I-1)*N))
      ENDDO



C           Integrate the mu times the radiance to find the fluxes
c jm: we don't care about the fluxes so far. so, just skip this.
c      DO L = 1, NUM_LAYERS+1
c        DO I = 1, NSTOKES
c          K = I+NSTOKES*(L-1)
c          UP_FLUX(K) = 0.0
c          DOWN_FLUX(K) = 0.0
c          DO J = 1, NUMMU
c            UP_FLUX(K) = UP_FLUX(K)
c     .               + TWOPI*QUAD_WEIGHTS(J) * MU_VALUES(J)
c     .               * UP_RAD(I+NSTOKES*(J-1)+N*(L-1))
c            DOWN_FLUX(K) = DOWN_FLUX(K)
c     .               + TWOPI*QUAD_WEIGHTS(J) * MU_VALUES(J)
c     .               * DOWN_RAD(I+NSTOKES*(J-1)+N*(L-1))
c          ENDDO
c        ENDDO
c      ENDDO

      RETURN
      END

