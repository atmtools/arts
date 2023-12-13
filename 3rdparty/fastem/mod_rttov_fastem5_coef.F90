!
MODULE mod_rttov_fastem5_coef
! Description:
!> @file
!!   Contains data for the FASTEM-4,5,6 MW sea surface emissivity models
!
!> @brief
!!   Contains data for the FASTEM-4,5,6 MW sea surface emissivity models
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
  USE parkind1, Only : fp => jprb
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  PUBLIC :: FresnelVariables_type
  PUBLIC :: PermittivityVariables_type
  
  PUBLIC :: fp
  
  ! Environment setup
  ! -----------------
  ! Module use
!  INTEGER, PUBLIC, PARAMETER  :: fp = SELECTED_REAL_KIND(15)
  REAL(fp), PUBLIC, PARAMETER :: ZERO     = 0.0_fp
  REAL(fp), PUBLIC, PARAMETER :: POINT_5  = 0.5_fp
  REAL(fp), PUBLIC, PARAMETER :: ONE      = 1.0_fp
  REAL(fp), PUBLIC, PARAMETER :: TWO      = 2.0_fp
  REAL(fp), PUBLIC, PARAMETER :: THREE    = 3.0_fp
  REAL(fp), PUBLIC, PARAMETER :: PI = 3.141592653589793238462643383279_fp
  REAL(fp), PUBLIC, PARAMETER :: DEGREES_TO_RADIANS = PI/180.0_fp
  REAL(fp), PUBLIC, PARAMETER :: transmittance_limit_lower = 0.00001_fp
  REAL(fp), PUBLIC, PARAMETER :: transmittance_limit_upper = 0.9999_fp

  REAL(fp), PUBLIC, PARAMETER :: e0_4 = 0.0088419_fp         ! Value used in FASTEM-4
  REAL(fp), PUBLIC, PARAMETER :: e0_5 = 0.00885418781762_fp  ! from Paul van Delst (used in FASTEM-5)
  ! minimum and maximum frequency
  REAL( fp ), PUBLIC, PARAMETER ::  min_f = 1.4_fp
  REAL( fp ), PUBLIC, PARAMETER ::  max_f = 200.0_fp
  ! minimum and maximum wind speed
  REAL( fp ), PUBLIC, PARAMETER ::  min_wind = 0.3_fp
  REAL( fp ), PUBLIC, PARAMETER ::  max_wind = 35.0_fp

  ! The fitting coefficients for the JCSDA permittivity model
  ! see the ref. (An Improved Fast Microwave Sea Surface Emissivity Model, FASTEM4,
  ! Technical Report, UK Met. Office)
  REAL(fp), PUBLIC, PARAMETER :: A_COEF(0:38) = (/ 3.8_fp, 0.0248033_fp, 87.9181727_fp, &
    -0.4031592248_fp,    0.0009493088010_fp,  -0.1930858348E-05_fp, -0.002697_fp,       &
    -7.3E-06_fp,        -8.9E-06_fp,           5.723_fp,             0.022379_fp,       &
    -0.00071237_fp,     -6.28908E-03_fp,     1.76032E-04_fp,        -9.22144E-05_fp,    &
     0.1124465_fp,      -0.0039815727_fp,    0.00008113381_fp,      -0.00000071824242_fp,&
    -2.39357E-03_fp,     3.1353E-05_fp,       -2.52477E-07_fp,       0.003049979018_fp, &
    -3.010041629E-05_fp, 0.4811910733E-05_fp, -0.4259775841E-07_fp,  0.149_fp,          &
    -8.8E-04_fp,        -1.05E-04_fp,          2.033E-02_fp,         1.266E-04_fp,  &
     2.464E-06_fp,      -1.849E-05_fp,         2.551E-07_fp,        -2.551E-08_fp,  &
     0.182521_fp,       -1.46192E-03_fp,       2.09324E-05_fp,      -1.28205E-07_fp/)

  ! fitting coefficients for the large-scale correction - FASTEM-5
  REAL(fp),  PUBLIC, PARAMETER :: Lcoef5(36) = (/ &
  -5.994667E-02_fp, 9.341346E-04_fp,-9.566110E-07_fp, 8.360313E-02_fp,-1.085991E-03_fp, &
   6.735338E-07_fp,-2.617296E-02_fp, 2.864495E-04_fp,-1.429979E-07_fp,-5.265879E-04_fp, &
   6.880275E-05_fp,-2.916657E-07_fp,-1.671574E-05_fp, 1.086405E-06_fp,-3.632227E-09_fp, &
   1.161940E-04_fp,-6.349418E-05_fp, 2.466556E-07_fp,-2.431811E-02_fp,-1.031810E-03_fp, &
   4.519513E-06_fp, 2.868236E-02_fp, 1.186478E-03_fp,-5.257096E-06_fp,-7.933390E-03_fp, &
  -2.422303E-04_fp, 1.089605E-06_fp,-1.083452E-03_fp,-1.788509E-05_fp, 5.464239E-09_fp, &
  -3.855673E-05_fp, 9.360072E-07_fp,-2.639362E-09_fp, 1.101309E-03_fp, 3.599147E-05_fp, &
  -1.043146E-07_fp /)
  
  ! fitting coefficients for the large-scale correction - FASTEM-4
  REAL(fp),  PUBLIC, PARAMETER :: Lcoef4(36) = (/ &
  -9.197134E-02_fp, 8.310678E-04_fp,-6.065411E-07_fp, 1.350073E-01_fp,-1.032096E-03_fp, &
   4.259935E-07_fp,-4.373322E-02_fp, 2.545863E-04_fp, 9.835554E-08_fp,-1.199751E-03_fp, &
   1.360423E-05_fp,-2.088404E-08_fp,-2.201640E-05_fp, 1.951581E-07_fp,-2.599185E-10_fp, &
   4.477322E-04_fp,-2.986217E-05_fp, 9.406466E-08_fp,-7.103127E-02_fp,-4.713113E-05_fp, &
   1.754742E-06_fp, 9.720859E-02_fp, 1.374668E-04_fp,-2.591771E-06_fp,-2.687455E-02_fp, &
  -3.677779E-05_fp, 7.548377E-07_fp,-3.049506E-03_fp,-5.412826E-05_fp, 2.285387E-07_fp, &
  -2.201640E-05_fp, 1.951581E-07_fp,-2.599185E-10_fp, 2.297488E-03_fp, 3.787032E-05_fp, &
  -1.553581E-07_fp /)


  ! fitting coefficients for the small-scale correction
  REAL(fp),  PUBLIC, PARAMETER :: Scoef(8) = (/ &
    -5.0208480E-06_fp,   2.3297951E-08_fp,   4.6625726E-08_fp,  -1.9765665E-09_fp, &
    -7.0469823E-04_fp,   7.5061193E-04_fp,   9.8103876E-04_fp,   1.5489504E-04_fp /)

! August 3, 2011
! 0.5
  ! fitting coefficients for the downward radiation using transmittance - FASTEM-5
  REAL(fp),  PUBLIC, PARAMETER :: t_c5(45) = (/ &
    0.199277E+00_fp, 0.166155E+00_fp, 0.153272E-01_fp, 0.399234E+01_fp,-0.130968E+01_fp, &
   -0.874716E+00_fp,-0.169403E+01_fp,-0.260998E-01_fp, 0.540443E+00_fp,-0.282483E+00_fp, &
   -0.219994E+00_fp,-0.203438E-01_fp, 0.351731E+00_fp, 0.208641E+01_fp,-0.693299E+00_fp, &
    0.867861E-01_fp, 0.619020E-01_fp, 0.595251E-02_fp,-0.475191E+01_fp,-0.430134E-01_fp, &
    0.248524E+01_fp, 0.388242E-01_fp, 0.194901E+00_fp,-0.425093E-01_fp, 0.607698E+01_fp, &
   -0.313861E+01_fp,-0.103383E+01_fp,-0.377867E+01_fp, 0.180284E+01_fp, 0.699556E+00_fp, &
   -0.506455E-01_fp,-0.262822E+00_fp, 0.703056E-01_fp, 0.362055E+01_fp,-0.120318E+01_fp, &
   -0.124971E+01_fp, 0.154014E-01_fp, 0.759848E-01_fp,-0.268604E-01_fp,-0.802073E+01_fp, &
    0.324658E+01_fp, 0.304165E+01_fp, 0.100000E+01_fp, 0.200000E-01_fp, 0.300000E+00_fp /)
  
  ! fitting coefficients for the downward radiation using transmittance - FASTEM-4
  REAL(fp),  PUBLIC, PARAMETER :: t_c4(45) = (/ &
  -0.675700E-01_fp, 0.214600E+00_fp,-0.363000E-02_fp, 0.636730E+01_fp, 0.900610E+00_fp, &
  -0.524880E+00_fp,-0.370920E+01_fp,-0.143310E+01_fp, 0.397450E+00_fp, 0.823100E-01_fp, &
  -0.255980E+00_fp, 0.552000E-02_fp, 0.208000E+01_fp, 0.244920E+01_fp,-0.456420E+00_fp, &
  -0.224900E-01_fp, 0.616900E-01_fp,-0.344000E-02_fp,-0.507570E+01_fp,-0.360670E+01_fp, &
   0.118750E+01_fp, 0.124950E+00_fp, 0.121270E+00_fp, 0.714000E-02_fp, 0.736620E+01_fp, &
  -0.114060E+00_fp,-0.272910E+00_fp,-0.504350E+01_fp,-0.336450E+00_fp, 0.161260E+00_fp, &
  -0.154290E+00_fp,-0.141070E+00_fp,-0.809000E-02_fp, 0.395290E+01_fp, 0.958580E+00_fp, &
  -0.159080E+00_fp, 0.368500E-01_fp, 0.307100E-01_fp, 0.810000E-03_fp,-0.619960E+01_fp, &
  -0.172580E+01_fp, 0.641360E+00_fp, 0.100000E+01_fp, 0.200000E-01_fp, 0.300000E+00_fp /)

  ! fitting coefficients for the azimuth dependence, S1 and S4 trained by
  REAL(fp), PUBLIC, PARAMETER :: b_coef(120) = (/ &
   3.307255E-04_fp,-2.901276E-06_fp,-1.475497E-04_fp, 1.288152E-06_fp, 1.004010E-04_fp, &
  -2.671158E-07_fp, 4.363154E-06_fp,-9.817795E-09_fp,-4.777876E-05_fp, 3.051852E-08_fp, &
   1.369383E-03_fp,-2.215847E-05_fp,-8.099833E-04_fp, 1.767702E-05_fp,-5.977649E-06_fp, &
  -1.784656E-07_fp,-9.355531E-07_fp, 5.495131E-08_fp,-3.479300E-05_fp,-3.751652E-07_fp, &
   2.673536E-04_fp,-1.378890E-06_fp,-8.660113E-05_fp, 2.871488E-07_fp, 1.361118E-05_fp, &
  -1.622586E-08_fp,-1.232439E-07_fp,-3.067416E-09_fp,-1.835366E-06_fp, 8.098728E-09_fp, &
   1.255415E-04_fp,-5.145201E-07_fp,-8.832514E-06_fp,-5.105879E-09_fp, 2.734041E-05_fp, &
  -3.398604E-07_fp, 3.417435E-06_fp,-7.043251E-09_fp, 1.497222E-05_fp,-6.832110E-09_fp, &
  -2.315959E-03_fp,-1.023585E-06_fp, 5.154471E-05_fp, 9.534546E-06_fp,-6.306568E-05_fp, &
  -4.378498E-07_fp,-2.132017E-06_fp, 1.612415E-08_fp,-1.929693E-06_fp,-6.217311E-09_fp, &
  -1.656672E-04_fp, 6.385099E-07_fp, 2.290074E-06_fp, 1.103787E-07_fp,-5.548757E-06_fp, &
   5.275966E-08_fp,-4.653774E-07_fp, 1.427566E-09_fp,-3.197232E-06_fp,-4.048557E-09_fp, &
  -1.909801E-04_fp,-3.387963E-07_fp, 4.641319E-05_fp, 4.502372E-07_fp,-5.055813E-05_fp, &
   2.104201E-07_fp,-4.121861E-06_fp,-1.633057E-08_fp,-2.469888E-05_fp, 4.492103E-08_fp, &
  -4.582853E-03_fp,-5.373940E-06_fp, 9.713047E-04_fp, 1.783009E-05_fp,-4.539091E-04_fp, &
   7.652954E-07_fp,-6.708905E-06_fp, 2.148401E-08_fp, 8.054350E-05_fp, 3.069258E-07_fp, &
  -6.405746E-05_fp,-9.694284E-08_fp, 1.914498E-05_fp, 1.336975E-07_fp,-4.561696E-06_fp, &
   3.769169E-08_fp,-6.105244E-07_fp, 2.433761E-10_fp,-3.961735E-06_fp, 1.995636E-08_fp, &
   1.350148E-06_fp, 3.678149E-07_fp, 1.261701E-05_fp,-2.011440E-07_fp,-2.361347E-05_fp, &
   2.943147E-08_fp,-1.304551E-07_fp,-1.119368E-09_fp, 8.469458E-06_fp,-2.292171E-09_fp, &
   1.419156E-03_fp,-3.838338E-06_fp, 8.222562E-05_fp,-1.106098E-06_fp,-5.482327E-05_fp, &
   3.083137E-07_fp, 4.418828E-06_fp,-1.302562E-08_fp, 3.768883E-05_fp,-5.012753E-08_fp, &
  -9.396649E-06_fp, 2.764698E-07_fp, 1.745336E-05_fp,-1.427031E-07_fp,-3.879930E-06_fp, &
  -1.117458E-08_fp, 5.688281E-08_fp, 1.513582E-09_fp, 6.778764E-06_fp,-7.691286E-09_fp /)

    ! Frequency-dependent azimuth correction
    REAL(fp), PUBLIC, PARAMETER :: x(9) = (/ 0.0_fp, 1.4_fp, 6.8_fp, 10.7_fp, 19.35_fp, &
                                   37._fp, 89._fp, 150._fp, 200._fp/)
    REAL(fp), PUBLIC, PARAMETER :: y(9) = (/ 0.0_fp, 0.1_fp, 0.6_fp, 0.9_fp, 1._fp, &
                                   1.0_fp, 0.4_fp, 0.2_fp, 0.0_fp/)

  ! Coefficients for M.Kazumori azimuth model function (FASTEM-6)
  REAL(fp), PUBLIC, PARAMETER :: coef_mk_azi(6,6,2) = RESHAPE( (/ &
  4.401E-02, -1.636E+01,  1.478E+00, -4.800E-02,  3.202E-06, -6.002E-05,&    ! 06V OK
  4.379E-02, -1.633E+01,  1.453E+00, -4.176E-02,  5.561E-06, -4.644E-05,&    ! 10V OK
  5.009E-02, -1.638E+01,  1.520E+00, -3.994E-02,  1.330E-05,  1.113E-05,&    ! 19V OK
  5.165E-02, -1.638E+01,  1.543E+00, -4.066E-02,  1.494E-05,  1.010E-05,&    ! 23V interpolated
  5.553E-02, -1.638E+01,  1.602E+00, -4.246E-02,  1.903E-05,  7.524E-06,&    ! 37V OK
 -9.131E-05,  1.251E+00,  6.769E-01, -2.913E-02,  1.092E+00, -1.806E-04,&    ! 89V OK revised
 -1.234E-07, -8.179E-03, -1.040E+01,  4.477E-01,  0.000E+00,  3.390E-05,&    ! 06H OK
 -1.938E-05, -8.007E-03, -1.039E+01,  4.610E-01,  0.000E+00,  4.419E-05,&    ! 10H OK
  1.362E-04, -1.013E-03, -9.235E+00,  3.844E-01,  0.000E+00,  2.891E-04,&    ! 19H OK
  1.519E-04, -7.865E-04, -9.234E+00,  3.884E-01,  0.000E+00,  6.856E-04,&    ! 23H Interpolated
  1.910E-04, -2.224E-04, -9.232E+00,  3.982E-01,  0.000E+00,  1.673E-03,&    ! 37H OK
  3.554E-04,  5.226E-04,  9.816E-01, -7.783E-03,  0.000E+00,  2.437E+01 /), &! 89H OK revised
  (/6,6,2/))



  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
 ! =============================================================
  ! Routine for forward model foam reflectivity
  ! Function dependence is on zenith angle only
  ! so no TL or AD routine.
  ! See Eqns(18.44a) (18.44b) in
  !   Ulaby, F.T. et al. (1986) Microwave Remote Sensing, Active
  !     and Passive, vol.3, From Theory to Applications, pp1457.
  ! =============================================================
  REAL(fp), PUBLIC, PARAMETER :: FR_COEFF(5) = &
    (/ 0.07_fp, -1.748e-3_fp, -7.336e-5_fp, 1.044e-7_fp, -0.93_fp /)

  !> Structure holding calculation results for the Fresnel coefficient
  !! calculation
  TYPE :: FresnelVariables_type
!    PRIVATE
    ! The intermediate terms
    COMPLEX(fp) :: z1, z2
    ! The real and imaginary components
    REAL(fp)    :: rzRv,izRv  ! Vertical
    REAL(fp)    :: rzRh,izRh  ! Horizontal
  END TYPE FresnelVariables_type


  !> Structure holding calculation results for the Ellison et al (2003)
  !! permittivity model
  TYPE :: PermittivityVariables_type
!    PRIVATE
    REAL(fp) :: t, t_sq, t_cu                   ! Temperature in degC
    REAL(fp) :: f1, f2, del1, del2, tau1_k, tau2_k, es_k, e1_k
    REAL(fp) :: ces, ctau1, ctau2, ce1, delta, beta, sigma, S
  END TYPE PermittivityVariables_type

END MODULE mod_rttov_fastem5_coef
!
