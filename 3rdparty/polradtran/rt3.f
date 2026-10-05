      PROGRAM RT3
C       RT3 solves the  plane parallel case of polarized monochromatic
C     radiative transfer for isotropic media.  This model is described
C     in K. F. Evans and G. L. Stephens, 1991: A New Polarized Atmospheric 
C     Radiative Transfer Model, J. Quant. Spectrosc. Radiat. Transfer, 
C     v. 46, no. 5, pp. 413-423, 1991.
C
C       This model now allow the output of radiances at any level in the
C     input layer file.  The output format has also been changed.  May 1995
C       Delta-M scaling added. July 1996
C
C        Frank Evans,  University of Colorado, Boulder,  July, 1996
C
C                          Method
C         The radiation field may have full angular dependence (polar and
C     azimuthal angles).  The scattering particles are assumed to be
C     randomly oriented (isotropic media) and to have a plane of symmetry.
C     The layers are uniform and infinite in horizontal extent and may be 
C     any thickness.  The properties of the layers are read from a file.
C     Solar and/or thermal sources of radiation are treated.
C         For randomly oriented particles with a plane of symmetry there are
C     six unique elements (out of 16) in the scattering phase matrix.
C     The phase matrix information is read in from a file containing the
C     coefficients for Legendre polynomial expansions for the six elements.
C     Linear thermal emission within each layer is calculated.
C     Thermal emission and reflection from a Lambertian or Fresnel
C     ground surface is incorporated.
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
C     See radtran3.f for explanation of input parameters.  Note that 
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
C     The format of the scattering file can be found from the subroutine
C     READ_SCAT_FILE or from an example.
C
C         The calculation is performed sequentially for each azimuth mode.
C     Before the loop over the azimuth modes, the Fourier modes of the
C     phase matrix for each scattering file is made. The phase matrices
C     are stored temporily either in memory.  Then there is a loop for 
C     the azimuth modes, and inside that there is a loop for each layer. 
C     For each layer the phase matrix and pseudo-source vector is retrieved, 
C     and the normalization of the phase matrix is checked.  The thermal
C     source vector is made if needed.  The parameters for the doubling 
C     of the sources are made.  The source is linear in optical depth 
C     for the thermal case and exponential in optical depth for the 
C     solar case.  For each layer the infinitesimal generator initialization
C     is used to calculate the local reflection and transmission matrices
C     and source vectors for a very thin sublayer.  A doubling subroutine
C     is then called to calculate the matrices and vectors for the whole
C     layer.  If the layer doesn't scatter then a subroutine calculates the
C     reflection and transmission matrices and source vector rather than
C     using initialization and doubling.  The reflection and transmission 
C     matrices and source vectors for each layer are stored in memory.
C
C         After the doubling has been done to find the properties of all 
C     the layers, there is a loop over the output levels.  For each
C     output level an adding subroutine is called to combine the layers
C     above and below the output level.  Then an internal radiance subroutine
C     is called to compute the radiance at the output level from the
C     reflection and transmission matrices and source vectors for the
C     medium above and for below and the incident radiance.  There is 
C     assumed to be thermal radiance from above and thermal and/or 
C     reflected direct solar radiance from the lower surface.  The 
C     reflection from the lower surface is simply treated as another 
C     layer in the medium (with unity transmission and no source).  The
C     RADTRAN subroutine computes the Fourier azimuthal modes for the
C     discrete quadrature zenith angles at each output level, which must
C     be at the input layer boundaries.
C
C         There are four types of numerical quadrature schemes available.
C     Gaussian, double Gaussian, and Lobatto are standard.  The 'extra-angle' 
C     is the same as gaussian quadrature but with extra angles added in. The 
C     weights for the extra angles are zero, so the radiances calculated at 
C     the gaussian quadrature angles are uneffected.  The radiances at the 
C     extra angles are effectively interpolated.  The user defined quadrature
C     was removed because of instabilities.
C
C         The output is a text file that contains the parameter values at
C     the beginning followed by the radiance and flux values.  The 
C     polarized radiance may be output as I, Q, U, V Stokes parameters or 
C     as V, H, U, V.  The radiance may be also converted to 
C     brightness temperature, though it is always first computed using 
C     the Planck function, in Watts/(meter^2 ster micron).  The brightness
C     temperature may be effective blackbody (UNITS=T) or Rayleigh-Jeans (R).
C     The output Stokes parameter are listed together for each angle 
C     and height.  A mu of +2 or -2 indicates the hemispheric flux value.  
C     Positive mu values are downwelling, and negative are upwelling angles.
C     For the output the Fourier azimuth series is summed and the radiances
C     output at a specified number of evenly spaced azimuthal angles between
C     0 and 180 degrees inclusive (NUMAZIMUTHS=3 means 0, 90, 180 degrees).
C
C                        Program Structure
C         The radiative transfer program is contained in six files.
C     rt3.f has the main program, the user input routine, the layer reading
C     routine, and the output routine. radtran3.f has the RADTRAN subroutine
C     which performs all of the radiative transfer calculation. It may be
C     called separately, and its parameter passing is documented
C     in radtran3.f.  radscat3.f has routines that convert the Legendre series
C     scattering file data into the phase matrices.  radintg3.f contains
C     routines that perform the initialization, doubling, and adding.
C     radutil3.f has the Lambertian and Fresnel surface, Planck function, and
C     quadrature routines.  radmat.f contains general purpose matrix routines.
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
C     parameters.  As described in the paper the Stokes vector is the 
C     cosine azimuth modes of I and Q, and the sine modes of U and V 
C     (Ic,Qc,Us,Vs).  The second dimension is the quadrature angles for the 
C     mu (cosine theta) parameter with NUMMU elements.
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
C     MAXM is the square of MAXV.  MAXA is the maximum number of azimuth
C     modes.  MAXLEG is the maximum number of terms in the Legendre series.
C     MAXLAY is the maximum number of layers.  MAXLM is the size of the
C     largest arrays which hold the reflection and transmission matrices
C     for each layer.  The following table gives the location of array 
C     size parameters:
C         rt3.f             MAXV, MAXA, MAXLAY
C         radtran3.f        MAXV, MAXM, MAXLM, MAXLEG, MAXLAY
C         radscat3.f        MAXLEG
C         radintg3.f        MAXV, MAXM
C     The RADTRAN subroutine also has array sizes for the scattering
C     matrices temporary buffer (MAXSBUF) and the direct source temporary
C     buffer (MAXDBUF).  
C
C     The fractional accuracy of the output radiances is about the size
C     of the MAX_DELTA_TAU parameter.  It is highly recommended that
C     double precision (15 decimals) rather than single precision
C     (7 decimals) be used.
C
C     It is important that the phase matrix be normalized so energy is
C     conserved, by having enough quadrature angles for the number of
C     Legendre terms in the scattering file.  This limit on the number of
C     Legendre terms is enforced by truncating the series at the appropriate
C     degree (4*NUMMU-3 for Gaussian quadrature) and printing out a warning
C     message giving the scattering file name.  The Legendre series truncation
C     should prevent the normalization check from failing, but the
C     normalization check is not always sufficient for energy to be conserved.
C     If truncation occurs then the true shape of the phase function is
C     distorted, which may not be significant if the truncation is modest.
C     To avoid truncation the number of quadrature angles may be increased,
C     but note that the CPU time goes as the cube of number of angles.
C     Delta-M scaling may be used to avoid truncation for highly peaked 
C     phase functions.  Note, that for collimated (solar) problems the 
C     delta-M solution will be accurate for flux, but may still have 
C     radiances which oscillate unphysically (though less than with 
C     truncating the original phase function). 
C     
C
C     Specular surface reflection is not treated so only the Lambertian
C     surface may be used with a solar source.
C
C
C     Routine locations:
C          File         Routines
C        rt3.f          READ_LAYERS, USER_INPUT, OUTPUT_FILE
C        radtran3.f     RADTRAN
C        radutil3.f     LAMBERT_SURFACE, LAMBERT_RADIANCE,
C                       FRESNEL_SURFACE, FRESNEL_RADIANCE,
C                       THERMAL_RADIANCE, PLANCK_FUNCTION,
C                       GAUSS_LEGENDRE_QUADRATURE,
C                       DOUBLE_GAUSS_QUADRATURE, LOBATTO_QUADRATURE
C        radscat3.f     READ_SCAT_FILE, SCATTERING, GET_SCATTER,
C                       CHECK_NORM, DIRECT_SCATTERING, GET_DIRECT
C                       (SUM_LEGENDRE, NUMBER_SUMS,
C                        ROTATE_PHASE_MATRIX, SUM_MATRIX)
C        radintg3.f     INITIALIZE, INITIAL_SOURCE,
C                       NONSCATTER_LAYER, INTERNAL_RADIANCE,
C                       DOUBLING_INTEGRATION, COMBINE_LAYERS
C        radmat.f       MCOPY, MADD, MSUB, MSCALARMULT, MZERO, MDIAG,
C                       MIDENTITY, MTRANSPOSE, MMULT, MINVERT
C
C

      INTEGER   MAXV, MAXA, MAXLAY
      PARAMETER (MAXV=64, MAXA=32)
      PARAMETER (MAXLAY=200)

      INTEGER   NSTOKES, NUMMU, AZIORDER
      INTEGER   NUM_LAYERS, SRC_CODE
      INTEGER   NOUTLEVELS, OUTLEVELS(MAXLAY), NUMAZIMUTHS
      REAL*8    GROUND_TEMP, GROUND_ALBEDO
      COMPLEX*16  GROUND_INDEX
      REAL*8    SKY_TEMP, WAVELENGTH, MAX_DELTA_TAU
      REAL*8    DIRECT_FLUX, DIRECT_MU
      REAL*8    MU_VALUES(MAXV)
      REAL*8    HEIGHT(MAXLAY), TEMPERATURES(MAXLAY)
      REAL*8    GAS_EXTINCT(MAXLAY)
      REAL*8    UP_FLUX(4*MAXLAY), DOWN_FLUX(4*MAXLAY)
      REAL*8    UP_RAD(MAXLAY*MAXA*MAXV), DOWN_RAD(MAXLAY*MAXA*MAXV)
      CHARACTER QUAD_TYPE*1, DELTAM*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
      CHARACTER*64 LAYER_FILE, OUT_FILE
      CHARACTER*64 SCAT_FILES(MAXLAY)



      CALL USER_INPUT (NSTOKES, NUMMU, AZIORDER, MU_VALUES,
     .                    SRC_CODE, LAYER_FILE, OUT_FILE,
     .                    QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS)


      IF (1+AZIORDER .GT. MAXA) THEN
          WRITE (*,*) 'Maximum number of azimuth modes exceeded.'
          STOP
      ENDIF


      CALL READ_LAYERS (LAYER_FILE, MAXLAY, NUM_LAYERS,
     .                  HEIGHT, TEMPERATURES,
     .                  GAS_EXTINCT, SCAT_FILES)


      MAX_DELTA_TAU = 1.0E-6
      CALL RADTRAN (NSTOKES, NUMMU, AZIORDER, MAX_DELTA_TAU,
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



      CALL OUTPUT_FILE (NSTOKES, NUMMU, AZIORDER,
     .                    SRC_CODE, LAYER_FILE, OUT_FILE,
     .                    QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NUM_LAYERS, HEIGHT,
     .                    NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,
     .                    MU_VALUES, UP_FLUX, DOWN_FLUX,
     .                    UP_RAD, DOWN_RAD)

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






      SUBROUTINE USER_INPUT (NSTOKES, NUMMU, AZIORDER, MU_VALUES,
     .                    SRC_CODE, LAYER_FILE, OUT_FILE,
     .                    QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS)
      INTEGER  NSTOKES, NUMMU, AZIORDER, SRC_CODE
      INTEGER  NOUTLEVELS, OUTLEVELS(*), NUMAZIMUTHS
      REAL*8   MU_VALUES(*), GROUND_TEMP, GROUND_ALBEDO
      REAL*8   SKY_TEMP, WAVELENGTH
      REAL*8   DIRECT_FLUX, DIRECT_MU
      COMPLEX*16  GROUND_INDEX
      CHARACTER QUAD_TYPE*1, DELTAM*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
      CHARACTER*(*) LAYER_FILE, OUT_FILE
      REAL*8   THETA
      INTEGER  I

      WRITE (*,'(1X,A)')  'Number of Stokes parameters (1 - 4) : '
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
50        CONTINUE
              WRITE (*,'(1X,A)') 'Mu value : '
              READ (*,*) MU_VALUES(I)
              I = I - 1
          IF (MU_VALUES(I+1) .NE. 0.0) GOTO 50
      ENDIF

      WRITE (*,'(1X,A)')  'Order of azimuth expansion (0,1,...) : '
      READ (*,*)  AZIORDER

      WRITE (*,'(1X,A)') 'Layers data file name : '
      READ (*,'(A)') LAYER_FILE

      WRITE (*,'(1X,A)') 'Delta-M scaling (Y or N) : '
      READ (*,'(A)')  DELTAM

      WRITE (*,'(1X,A)')
     .   'Source code (none=0, solar=1, thermal=2, both=3) : '
      READ (*,*)  SRC_CODE
      SRC_CODE = MIN0( MAX0( SRC_CODE, 0), 3)

      DIRECT_MU = 1.0
      IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
          WRITE (*,'(1X,A)')  'Direct flux (W/(m*m)/um or K) : '
          READ (*,*)  DIRECT_FLUX
          WRITE (*,'(1X,A)')
     .         'Direct flux direction (zenith angle) (deg) : '
          READ (*,*)  THETA
          DIRECT_MU = DABS(DCOS(0.017453292D0*(THETA)))
      ENDIF

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
      WRITE (*,'(1X,A)')  'Number of output azimuths : '
      READ (*,*)  NUMAZIMUTHS

      WRITE (*,'(1X,A)') 'Output data file name : '
      READ (*,'(A)') OUT_FILE

      RETURN
      END





      SUBROUTINE OUTPUT_FILE (NSTOKES, NUMMU, AZIORDER,
     .                    SRC_CODE, LAYER_FILE, OUT_FILE,
     .                    QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NUM_LAYERS, HEIGHT,
     .                    NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,
     .                    MU_VALUES, UP_FLUX, DOWN_FLUX,
     .                    UP_RAD, DOWN_RAD)
      INTEGER  NSTOKES, NUMMU, NUMAZI, AZIORDER, SRC_CODE, NUM_LAYERS
      INTEGER  NOUTLEVELS, OUTLEVELS(*), NUMAZIMUTHS
      REAL*8   GROUND_TEMP, GROUND_ALBEDO
      REAL*8   SKY_TEMP, WAVELENGTH
      REAL*8   DIRECT_FLUX, DIRECT_MU
      REAL*8   HEIGHT(NUM_LAYERS+1)
      REAL*8   MU_VALUES(NUMMU)
      REAL*8   UP_FLUX(NSTOKES,NOUTLEVELS)
      REAL*8   DOWN_FLUX(NSTOKES,NOUTLEVELS)
      REAL*8   UP_RAD(NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)
      REAL*8   DOWN_RAD(NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)
      COMPLEX*16  GROUND_INDEX
      CHARACTER*(*) LAYER_FILE, OUT_FILE
      CHARACTER QUAD_TYPE*1, DELTAM*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
      CHARACTER*32 QUAD_NAME, UNITS_NAME, GROUND_NAME
      CHARACTER*64 FORM1
      INTEGER  I, J, K, L, LI, M, N
      REAL*4   OUT(4), PHI, PHID, PI
      PARAMETER (PI=3.1415926535897932384D0)


      N = NUMMU*(AZIORDER+1)*NOUTLEVELS
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, N, 
     .                     WAVELENGTH, 0, UP_RAD)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, N, 
     .                     WAVELENGTH, 0, DOWN_RAD)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS, 
     .                     WAVELENGTH, 1, UP_FLUX)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS, 
     .                     WAVELENGTH, 1, DOWN_FLUX)

      NUMAZI = 2*AZIORDER+1
      IF (NSTOKES .LE. 2) NUMAZI = AZIORDER+1
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
     .                'C  NUMMU=', NUMMU,  '  NUMAZI=',NUMAZI,
     .                '  AZIORDER=',AZIORDER, '  NSTOKES=',NSTOKES
      WRITE (3,'(A,A32,A,A1)')
     .                'C  LAYER_FILE=',    LAYER_FILE,
     .                '   DELTA-M=',DELTAM
      WRITE (3,'(A,I1,A,A16)')
     .                'C  SRC_CODE=',      SRC_CODE,
     .                '   QUAD_TYPE=',     QUAD_NAME
      IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
          WRITE (3,'(A,E11.5,A,F8.6)')
     .                'C  DIRECT_FLUX=',   DIRECT_FLUX,
     .                '   DIRECT_MU=',     DIRECT_MU
      ENDIF
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
          FORM1 = '(F8.3,1X,F5.1,1X,F8.5,4(1X,F7.2),:)'
      ELSE
          FORM1 = '(F8.3,1X,F5.1,1X,F8.5,4(1X,E13.6),:)'
      ENDIF
 
      IF (OUTPOL .EQ. 'VH') THEN
        WRITE (3,'(A)') 
     .    'C    Z      PHI     MU    FLUX/RADIANCE (V,H,U,V)'
      ELSE
        WRITE (3,'(A)') 
     .    'C    Z      PHI     MU    FLUX/RADIANCE (I,Q,U,V)'
      ENDIF
 
      DO L = 1, NOUTLEVELS
        LI = OUTLEVELS(L)
C               Output fluxes at this level
        WRITE (3,FORM1) HEIGHT(LI), 0., -2.0,
     .        (SNGL(UP_FLUX(I,L)),I=1,NSTOKES)
        WRITE (3,FORM1) HEIGHT(LI), 0., +2.0,
     .        (SNGL(DOWN_FLUX(I,L)),I=1,NSTOKES)
 
C               For each azimuth and zenith at this level sum the Fourier
C               azimuth series appropriate for the particular Stokes parameter
C               and output the radiance.
        DO K = 1, NUMAZIMUTHS
          IF (NUMAZIMUTHS .EQ. 1) THEN
            PHID = 0.0
          ELSE
            PHID = 180.0*FLOAT(K-1)/(NUMAZIMUTHS-1)
          ENDIF
          PHI = PI*PHID/180.0
C               Output upwelling radiance: -1 < mu < 0
          DO J = NUMMU, 1, -1
            DO I = 1, NSTOKES
              OUT(I) = 0.0
              DO M = 0, AZIORDER
                IF (I .LE. 2) THEN
                  OUT(I) = OUT(I) + COS(M*PHI)*UP_RAD(I,J,M+1,L)
                ELSE
                  OUT(I) = OUT(I) + SIN(M*PHI)*UP_RAD(I,J,M+1,L)
                ENDIF
              ENDDO
            ENDDO
            WRITE (3,FORM1) HEIGHT(LI), PHID, -MU_VALUES(J),
     .                      (OUT(I),I=1,NSTOKES)
          ENDDO
C               Output downwelling radiance: 0 < mu < 1
          DO J = 1, NUMMU
            DO I = 1, NSTOKES
              OUT(I) = 0.0
              DO M = 0, AZIORDER
                IF (I .LE. 2) THEN
                  OUT(I) = OUT(I) + COS(M*PHI)*DOWN_RAD(I,J,M+1,L)
                ELSE
                  OUT(I) = OUT(I) + SIN(M*PHI)*DOWN_RAD(I,J,M+1,L)
                ENDIF
              ENDDO
            ENDDO
            WRITE (3,FORM1) HEIGHT(LI), PHID,  MU_VALUES(J),
     .                      (OUT(I),I=1,NSTOKES)
          ENDDO
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



