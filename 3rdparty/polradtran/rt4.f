      PROGRAM RT4
C       RT4 solves the plane parallel case of polarized monochromatic
C     radiative transfer for azimuthally symmetric scattering media.
C       This model is briefly described in Evans, K. F., and G. L. Stephens, 
C     1995: Microwave radiative transfer through clouds composed of 
C     realistically shaped ice crystals. Part II: Remote Sensing of Ice Clouds
C     J. Atmos. Sci., v. 52, 2058-2072.
C
C       This model now allow the output of radiances at any level in the
C     input layer file.  The output format has also been changed.
C
C        Frank Evans,  University of Colorado, Boulder,  May, 1995
C
C
C                          Method
C         This model is a modification of the model for randomly oriented
C     particles (RT3).  Because of the additional complexities of the 
C     doubling method associated with the polarization of the direct beam,
C     the solar source parts have been removed, so thermal emission is
C     the only source of radiation.  In a plane-parallel geometry with
C     isotropic thermal emission and particles having random azimuthal 
C     orientation, the radiation field is azimuthally symmetric.  This 
C     symmetry also implies that the U and V Stokes parameters are zero.
C
C     There may be many arbitrary layers which are uniform and infinite in 
C     horizontal extent and may be any thickness.  The properties of the 
C     layers are read from a file.  The single scattering properties are 
C     read in from other files containing the Stokes scattering (Meuller) 
C     matrix, extinction matrix, and emission vector for appropriate sets 
C     of discrete quadrature angles.  Linear thermal emission within each 
C     layer is calculated.  Thermal emission and reflection from a 
C     Lambertian or Fresnel ground surface is incorporated.
C         The doubling and adding technique is used to solve the plane-
C     parallel radiative transfer equation.  Each input layer is divided
C     into a number of homogeneous sublayers with each sublayer being
C     thin enough for the finite difference initialization to be accurate.
C     Infinitesimal generator initialization is used to relate the scattering
C     matrix to the reflection and transmission matrices. The sublayers
C     are integrated with the doubling algorithm.  For each desired output
C     level the transmission, reflection, and source of the layers above 
C     and below the level are combined with the adding algorithm.  The
C     internal radiances are computed from the properties of the layers
C     above and below and the incident radiance from the boundaries.
C
C                          Operation
C         First a subroutine is called to get user input.  A bunch of
C     parameters are input with (hopefully) self-explanatory prompts.
C     See radtran4.f for explanation of input parameters.  Note that
C     letter inputs (except for filenames) must be in uppercase and
C     only the first letter is used.
C     The parameters for the layers are read in from the layer file.
C     All files, input and output, are Fortran formatted (text) files.
C     Each line of that file contains the height, temperature, gaseous
C     extinction, and scattering file name. The height and temperature
C     are specified for the interfaces between the layers, while the
C     gaseous extinction, and the scattering file are specified for each
C     layer.  The layers should start at the top and go down.  The
C     scattering file name is a Fortran string in single quotes.
C     The format of the scattering file is documented in radtran4.f and 
C     can be found from the subroutine GET_SCAT_FILE or from an example.  
C     The oriented scattering file must have the same quadrature angle setup
C     as used in the code.
C
C         The operation is considerably simplified over that of RT3.  There
C     is only the lowest azimuthal Fourier mode to consider.  The single
C     scattering properties are read directly in a form that can be used
C     so there is no conversion from the phase function of scattering angle
C     to the quadrature coordinates.  For each layer the single scattering
C     properties are read in and the normalization is check (to make sure
C     that scattering plus absorption equals extinction). For each layer the 
C     infinitesimal generator initialization is used to calculate the 
C     local reflection and transmission matrices and source vectors for a 
C     very thin sublayer.  A doubling subroutine is then called to calculate 
C     the matrices and vectors for the whole layer.  If the layer doesn't 
C     scatter then a subroutine calculates the reflection and transmission 
C     matrices and source vector rather than using initialization and 
C     doubling.  The reflection and transmission matrices and source
C     vectors for each layer are stored in memory.
C
C         After the doubling has been done to find the properties of all 
C     the layers, there is a loop over the output levels.  For each
C     output level an adding subroutine is called to combine the layers
C     above and below the output level.  Then an internal radiance subroutine
C     is called to compute the radiance at the output level from the
C     reflection and transmission matrices and source vectors for the
C     medium above and for below and the incident radiance.  There is 
C     assumed to be thermal radiance from above and thermal radiance 
C     from the lower surface.  The reflection from the lower surface
C     is simply treated as another layer in the medium (with unity 
C     transmission and no source).  The radiances for each quadrature 
C     direction for the m=0 azimuth mode are output, as are the integrated 
C     fluxes, for each output level.  Note, that these output levels
C     have to be at input layer boundaries.  
C
C         There are four types of numerical quadrature schemes available.
C     Gaussian, double Gaussian, and Lobatto are standard.  The 'extra-angle' 
C     is the same as gaussian quadrature but with extra angles added in. The 
C     weights for the extra angles are zero, so the radiances calculated at 
C     the gaussian quadrature angles are uneffected.  The radiances at the 
C     extra angles are effectively interpolated.  
C
C         The output is a text file that contains the parameter values at
C     the beginning followed by the radiance and flux values.  The 
C     polarized radiance may be output as I and Q Stokes parameters or 
C     as V and H polarizations.  The radiance may be also converted to 
C     brightness temperature, though it is always first computed using 
C     the Planck function, in Watts/(meter^2 ster micron).  The brightness
C     temperature may be effective blackbody (UNITS=T) or Rayleigh-Jeans (R).
C     The output Stokes parameter are listed together for each angle 
C     and height.  A mu of +2 or -2 indicates the hemispheric flux value.  
C     Positive mu values are downwelling, and negative are upwelling angles.
C
C                        Program Structure
C         The radiative transfer program is contained in six files.
C     rt4.f has the main program, the user input routine, the layer reading
C     routine, and the output routine.  radtran4.f has the RADTRANO subroutine
C     which performs all of the radiative transfer calculation. It may be
C     called separately, and its parameter passing is documented in
C     radtran4.f.  radscat4.f reads in the oriented scattering file and 
C     checks the normalization.  radintg4.f contains routines that perform 
C     the initialization, doubling, and adding.  radutil4.f has the Lambertian 
C     and Fresnel surface, Planck function, and quadrature routines.  
C     radmat.f contains general purpose matrix routines.
C     The location of the individual subroutines is given below.
C
C         The Fortran used for the program is basically Fortran 77 with
C     a few major differences:  long descriptive variable and subroutine
C     names and the use of ENDDO.  In addition the floating point variables 
C     are declared REAL*8.  The program will compile and work under 
C     VMS and most Unix Fortrans. 
C
C                        Data Storage
C         The basis for the radiance vectors has NSTOKES*NUMMU elements.
C     The first dimension is the polarization vector, made up of the Stokes
C     parameters either just I or I and Q.  The second dimension is the 
C     quadrature angles for the mu (cosine zenith angle) parameter with 
C     NUMMU elements.
C         The scattering matrix variable is actually four matrices:
C     P++, P+-, P-+, and P-- in that order.  Note: positive angles are
C     in the direction of increasing optical depth (down).
C         All real numbers are REAL*8.  The vectors and matrices are
C     one-dimensional arrays in the main program, and adjustable array
C     sizing is used in the subroutines.  The main storage requirement
C     is the scattering, reflection, and transmission matrices for all 
C     the layers (6*MAXLAY*MAXM real numbers).  This is wasteful of 
C     memory, but is done for speed at computing radiance at many levels.
C
C                        Limits and Caveats
C         The various array size limits are in Fortran parameter statements.
C     MAXV is the maximum size of the basis vector (NSTOKES*NUMMU), while
C     MAXM is the square of MAXV.  MAXLAY is the maximum number of layers.
C     MAXLM is the size of the largest arrays which hold the reflection 
C     and transmission matrices for each layer.  The following table gives 
C     the location of array size parameters:
C         rt4              MAXV, MAXLAY
C         radtran4.f       MAXV, MAXM, MAXLM, MAXLAY
C         radintg4.f       MAXV, MAXM
C
C     The fractional accuracy of the output radiances is about the size
C     of the MAX_DELTA_TAU parameter.  It is highly recommended that
C     double precision (15 decimals) rather than single precision
C     (7 decimals) be used.
C
C     It is assumed that the scattering and extinction matrices that are 
C     read in are symmetric between up and down (P+-=P-+).  If this is
C     not true then the SYMMETRIC logical in RADTRANO should be set to FALSE.
C
C     Routine locations:
C          File         Routines
C        rt4.f          READ_LAYERS, USER_INPUT, OUTPUT_FILE
C        radtran4.f     RADTRANO
C        radutil4.f     LAMBERT_SURFACE, LAMBERT_RADIANCE,
C                       FRESNEL_SURFACE, FRESNEL_RADIANCE,
C                       THERMAL_RADIANCE, PLANCK_FUNCTION,
C                       GAUSS_LEGENDRE_QUADRATURE,
C                       DOUBLE_GAUSS_QUADRATURE, LOBATTO_QUADRATURE
C        radscat4.f     GET_SCAT_FILE, CHECK_NORM
C        radintg4.f     INITIALIZE, INITIAL_SOURCE,
C                       NONSCATTER_LAYER, INTERNAL_RADIANCE,
C                       DOUBLING_INTEGRATION, COMBINE_LAYERS
C        radmat.f       MCOPY, MADD, MSUB, MSCALARMULT, MZERO, MDIAG,
C                       MIDENTITY, MTRANSPOSE, MMULT, MINVERT
C
C

      INTEGER   MAXV, MAXLAY
      PARAMETER (MAXV=64)
      PARAMETER (MAXLAY=200)

      INTEGER   NSTOKES, NUMMU
      INTEGER   NUM_LAYERS
      INTEGER   NOUTLEVELS, OUTLEVELS(MAXLAY)
      REAL*8    GROUND_TEMP, GROUND_ALBEDO
      COMPLEX*16 GROUND_INDEX
      REAL*8    SKY_TEMP, WAVELENGTH, MAX_DELTA_TAU
      REAL*8    MU_VALUES(MAXV)
      REAL*8    HEIGHT(MAXLAY), TEMPERATURES(MAXLAY)
      REAL*8    GAS_EXTINCT(MAXLAY)
      REAL*8    UP_RAD(MAXV*(MAXLAY+1)), DOWN_RAD(MAXV*(MAXLAY+1))
      REAL*8    UP_FLUX(4*(MAXLAY+1)), DOWN_FLUX(4*(MAXLAY+1))
      CHARACTER QUAD_TYPE*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
      CHARACTER*64 LAYER_FILE, OUT_FILE
      CHARACTER*64 SCAT_FILES(MAXLAY)



      CALL USER_INPUT (NSTOKES, NUMMU, MU_VALUES,
     .                 LAYER_FILE, OUT_FILE,
     .                 QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,
     .                 GROUND_ALBEDO, GROUND_INDEX,
     .                 SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                 NOUTLEVELS, OUTLEVELS)


      CALL READ_LAYERS (LAYER_FILE, MAXLAY, NUM_LAYERS,
     .                  HEIGHT, TEMPERATURES,
     .                  GAS_EXTINCT, SCAT_FILES)

      MAX_DELTA_TAU = 1.0E-6
      CALL RADTRANO (NSTOKES, NUMMU, MAX_DELTA_TAU,
     .               QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,
     .               GROUND_ALBEDO, GROUND_INDEX,
     .               SKY_TEMP, WAVELENGTH,
     .               NUM_LAYERS, HEIGHT, TEMPERATURES,
     .               GAS_EXTINCT, SCAT_FILES,
     .               NOUTLEVELS, OUTLEVELS,
     .               MU_VALUES, UP_FLUX, DOWN_FLUX,
     .               UP_RAD, DOWN_RAD)


      CALL OUTPUT_FILE (NSTOKES, NUMMU,
     .                  LAYER_FILE, OUT_FILE,
     .                  QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,
     .                  GROUND_ALBEDO, GROUND_INDEX,
     .                  SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                  NUM_LAYERS, HEIGHT,
     .                  NOUTLEVELS, OUTLEVELS, 
     .                  MU_VALUES, UP_FLUX, DOWN_FLUX,
     .                  UP_RAD, DOWN_RAD)

      END







      SUBROUTINE READ_LAYERS (LAYER_FILE, MAXLAY, NUM_LAYERS,
     .                        HEIGHT, TEMPERATURES,
     .                        GAS_EXTINCT, SCAT_FILES)
      INTEGER  MAXLAY, NUM_LAYERS
      REAL*8   HEIGHT(*), TEMPERATURES(*)
      REAL*8   GAS_EXTINCT(*)
      CHARACTER*(*)  LAYER_FILE, SCAT_FILES(*)
      INTEGER   I

C           Read in height, temperature, gaseous extinction, and
C                 scattering file for the layers
      OPEN (UNIT=1, FILE=LAYER_FILE, STATUS='OLD')
      I = 1
100   CONTINUE
          READ (1,*,ERR=990,END=110) HEIGHT(I), TEMPERATURES(I),
     .                GAS_EXTINCT(I), SCAT_FILES(I)
          I = I + 1
          IF (I .EQ. MAXLAY) THEN
              WRITE (*,*) 'Too many layers'
              STOP
          ENDIF
      GOTO 100
110   CONTINUE
      CLOSE(1)
      NUM_LAYERS = I - 2
      RETURN

990   CONTINUE
      WRITE (*,*) 'Error reading layers data file'
      RETURN
      END





      SUBROUTINE USER_INPUT (NSTOKES, NUMMU, MU_VALUES,
     .                    LAYER_FILE, OUT_FILE,
     .                    QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NOUTLEVELS, OUTLEVELS)
      INTEGER  NSTOKES, NUMMU
      INTEGER  NOUTLEVELS, OUTLEVELS(*)
      REAL*8   MU_VALUES(*), GROUND_TEMP, GROUND_ALBEDO
      REAL*8   SKY_TEMP, WAVELENGTH
      COMPLEX*16  GROUND_INDEX
      CHARACTER*(*) QUAD_TYPE, UNITS, OUTPOL, GROUND_TYPE
      CHARACTER*(*) LAYER_FILE, OUT_FILE
      INTEGER  I

      WRITE (*,'(1X,A)')  'Number of Stokes parameters (1 - 2) : '
      READ (*,*)  NSTOKES
      WRITE (*,'(1X,A)')  'Number of quadrature directions : '
      READ (*,*)  NUMMU
      WRITE (*,'(1X,A)') 'Type of quadrature : '
      WRITE (*,'(1X,A)')
     .  '(Gaussian, Double-Gauss, Lobatto, Extra-angles) : '
      READ (*,'(A)')  QUAD_TYPE
      IF (QUAD_TYPE(1:1) .EQ. 'E') THEN
          WRITE (*,*) 'Enter extra quadrature mu values (end with 0):'
          I = NUMMU
55        CONTINUE
              WRITE (*,'(1X,A)') 'Mu value : '
              READ (*,*) MU_VALUES(I)
              I = I - 1
          IF (MU_VALUES(I+1) .NE. 0.0) GOTO 55
      ENDIF

      WRITE (*,'(1X,A)') 'Layers data file name : '
      READ (*,'(A)') LAYER_FILE

      WRITE (*,'(1X,A)')  'Ground temperature : '
      READ (*,*)  GROUND_TEMP
      WRITE (*,'(1X,A)')  'Ground type (Lambertian or Fresnel) : '
      READ (*,'(A)')  GROUND_TYPE
      IF (GROUND_TYPE(1:1) .EQ. 'F') THEN
          WRITE (*,'(1X,A)')
     .              'Complex index of refraction of ground : '
          READ (*,*)  GROUND_INDEX
      ELSE
          WRITE (*,'(1X,A)') 'Ground albedo : '
          READ (*,*)  GROUND_ALBEDO
      ENDIF
      WRITE (*,'(1X,A)')  'Sky temperature : '
      READ (*,*)  SKY_TEMP

      WRITE (*,'(1X,A)')  'Wavelength (microns) : '
      READ (*,*)  WAVELENGTH
      WRITE (*,'(1X,A)') 'Output radiance units :'
      WRITE (*,'(1X,A,A)') '(W-W/m^2 um sr, ',
     .     'T-EBB brightness temperature, R-Rayleigh-Jeans Tb) : '
      READ (*,'(A)')  UNITS
      WRITE (*,'(1X,A)') 'Output polarization (IQ or VH) : '
      READ (*,'(A)') OUTPOL

      WRITE (*,'(1X,A)')  'Number of output levels : '
      READ (*,*)  NOUTLEVELS
      WRITE (*,'(1X,A)')  'Output level numbers : '
      READ (*,*)  (OUTLEVELS(I), I = 1,NOUTLEVELS)
 
      WRITE (*,'(1X,A)') 'Output data file name : '
      READ (*,'(A)') OUT_FILE

      RETURN
      END







      SUBROUTINE OUTPUT_FILE (NSTOKES, NUMMU, 
     .                    LAYER_FILE, OUT_FILE,
     .                    QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NUM_LAYERS, HEIGHT,
     .                    NOUTLEVELS, OUTLEVELS, 
     .                    MU_VALUES, UP_FLUX, DOWN_FLUX,
     .                    UP_RAD, DOWN_RAD)
      INTEGER  NSTOKES, NUMMU, NUM_LAYERS
      INTEGER  NOUTLEVELS, OUTLEVELS(*)
      REAL*8   GROUND_TEMP, GROUND_ALBEDO
      REAL*8   SKY_TEMP, WAVELENGTH
      REAL*8   MU_VALUES(NUMMU)
      REAL*8   HEIGHT(NUM_LAYERS+1)
      REAL*8   UP_FLUX(NSTOKES,NOUTLEVELS)
      REAL*8   DOWN_FLUX(NSTOKES,NOUTLEVELS)
      REAL*8   UP_RAD(NSTOKES,NUMMU,NOUTLEVELS)
      REAL*8   DOWN_RAD(NSTOKES,NUMMU,NOUTLEVELS)
      COMPLEX*16  GROUND_INDEX
      CHARACTER*(*) LAYER_FILE, OUT_FILE
      CHARACTER  QUAD_TYPE*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
      CHARACTER*32 QUAD_NAME, UNITS_NAME, GROUND_NAME
      CHARACTER*64 FORM1
      INTEGER  I, J, L, LI


      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NUMMU*NOUTLEVELS, 
     .                     WAVELENGTH, 0, UP_RAD)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NUMMU*NOUTLEVELS, 
     .                     WAVELENGTH, 0, DOWN_RAD)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS, 
     .                     WAVELENGTH, 1, UP_FLUX)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS, 
     .                     WAVELENGTH, 1, DOWN_FLUX)

      QUAD_NAME = 'GAUSSIAN'
      IF (QUAD_TYPE .EQ. 'D')  QUAD_NAME = 'DOUBLEGAUSS'
      IF (QUAD_TYPE .EQ. 'L')  QUAD_NAME = 'LOBATTO'
      IF (QUAD_TYPE .EQ. 'E')  QUAD_NAME = 'EXTRA-ANGLES'
      UNITS_NAME = 'WATTS/(M^2 MICRON STER)'
      IF (UNITS .EQ. 'T') UNITS_NAME = 'KELVINS - EBB'
      IF (UNITS .EQ. 'R') UNITS_NAME = 'KELVINS - RJ'
      GROUND_NAME = 'LAMBERTIAN'
      IF (GROUND_TYPE .EQ. 'F')  GROUND_NAME = 'FRESNEL'

      OPEN (UNIT=3, FILE=OUT_FILE, STATUS='UNKNOWN')

C           Output the parameters
      WRITE (3,'(A,I3,A,I3,A,I3,A,I1)')
     .                'C  NUMMU=', NUMMU,  '  NUMAZI=',1,
     .                '  AZIORDER=',0, '  NSTOKES=',NSTOKES
      WRITE (3,'(A,A32)')
     .                'C  LAYER_FILE=',    LAYER_FILE
      WRITE (3,'(A,I1,A,A16)')
     .                'C  SRC_CODE=',      2,
     .                '   QUAD_TYPE=',     QUAD_NAME
      WRITE (3,'(A,F8.2,A,A16)')
     .                'C  GROUND_TEMP=',   GROUND_TEMP,
     .                '   GROUND_TYPE=',   GROUND_NAME
      IF (GROUND_TYPE(1:1) .EQ. 'F') THEN
          WRITE (3,'(A,2F9.4,A,F8.2)')
     .                'C  GROUND_INDEX=',  GROUND_INDEX,
     .                '   SKY_TEMP=',      SKY_TEMP
      ELSE
          WRITE (3,'(A,F8.5,A,F8.2)')
     .                'C  GROUND_ALBEDO=', GROUND_ALBEDO,
     .                '   SKY_TEMP=',      SKY_TEMP
      ENDIF
      WRITE (3,'(A,E12.6)') 'C  WAVELENGTH=',    WAVELENGTH
      WRITE (3,'(A,A25,A,A2)') 'C  UNITS='     ,    UNITS_NAME,
     .                '   OUTPUT_POLARIZATION=', OUTPOL  


      IF (UNITS(1:1) .EQ. 'T') THEN
          FORM1 = '(F8.3,1X,F8.5,2(1X,F7.2),:)'
      ELSE
          FORM1 = '(F8.3,1X,F8.5,2(1X,E13.6),:)'
      ENDIF
 
      IF (OUTPOL .EQ. 'VH') THEN
        WRITE (3,'(A,A)') 'C    Z       MU    FLUX/RADIANCE (V,H)'
      ELSE
        WRITE (3,'(A,A)') 'C    Z       MU    FLUX/RADIANCE (I,Q)'
      ENDIF
 
      DO L = 1, NOUTLEVELS
        LI = OUTLEVELS(L)
C               Output fluxes at this level
        WRITE (3,FORM1) HEIGHT(LI), -2.0,
     .          (SNGL(UP_FLUX(I,L)),I=1,NSTOKES)
        WRITE (3,FORM1) HEIGHT(LI), +2.0,
     .          (SNGL(DOWN_FLUX(I,L)),I=1,NSTOKES)
 
C           For each zenith at this level output the Stokes parameters.
C             Output upwelling radiance: -1 < mu < 0
        DO J = NUMMU, 1, -1
          WRITE (3,FORM1) HEIGHT(LI), -MU_VALUES(J),
     .                    (SNGL(UP_RAD(I,J,L)),I=1,NSTOKES)
        ENDDO
C             Output downwelling radiance: 0 < mu < 1
        DO J = 1, NUMMU
          WRITE (3,FORM1) HEIGHT(LI), MU_VALUES(J),
     .                    (SNGL(DOWN_RAD(I,J,L)),I=1,NSTOKES)
        ENDDO
      ENDDO

      CLOSE (3)

      RETURN
      END





      SUBROUTINE CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUT,
     .                           WAVELEN, FLUXCODE, OUTPUT)
C       Converts the output radiance or flux arrays to VH polarization
C     and effective blackbody temperature if desired.  OUTPOL='VH'
C     converts the polarization basis of the first two Stokes parameters
C     to vertical/horizontal polarization.  If UNITS='T' the radiance is
C     converted to effective blackbody brightness temperature, and if
C     UNITS='R' the radiance is converted to Rayleigh-Jeans brightness
C     temperature.  If the output is flux then FLUXCODE=1, and the flux 
C     is divided by pi before converting to brightness temperature.
      INTEGER NSTOKES, NOUT, FLUXCODE
      REAL*8  WAVELEN, OUTPUT(NSTOKES,NOUT)
      CHARACTER UNITS*1, OUTPOL*2
      INTEGER I, J
      REAL*8  IV, IH, RAD, TEMP

      DO J = 1, NOUT      
C           Convert to Vertical and Horizontal polarization if desired
        IF (OUTPOL .EQ. 'VH') THEN
          IV = 0.5*(OUTPUT(1,J) + OUTPUT(2,J))
          IH = 0.5*(OUTPUT(1,J) - OUTPUT(2,J))
          OUTPUT(1,J) = IV
          OUTPUT(2,J) = IH
        ENDIF
C           Convert to brightness temperature
        IF (UNITS .EQ. 'T' .OR. UNITS .EQ. 'R') THEN
          DO I = 1, NSTOKES
            RAD = OUTPUT(I,J)
            IF (OUTPOL .EQ. 'VH' .AND. I .LE. 2)  RAD = 2.0*RAD
            IF (FLUXCODE .EQ. 1)  RAD = RAD/ACOS(-1.0)
            IF (UNITS .EQ. 'R') THEN
              TEMP = RAD * WAVELEN**4 * 1.4388D4/1.1911D8
            ELSE
              IF (RAD .GT. 0.0) THEN
                TEMP = 1.4388D4 /
     .            (WAVELEN*DLOG(1.0+ 1.1911D8/(RAD*WAVELEN**5)))
              ELSE IF (RAD .EQ. 0.0) THEN
                TEMP = 0.0D0
              ELSE
                TEMP = -1.4388D4 /
     .            (WAVELEN*DLOG(1.0+ 1.1911D8/(-RAD*WAVELEN**5)))
              ENDIF
            ENDIF
            OUTPUT(I,J) = TEMP
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END



