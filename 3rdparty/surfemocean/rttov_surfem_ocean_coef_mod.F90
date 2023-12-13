! Description:
!> @file
!!   Contains data for the SURFEM-Ocean MW sea surface emissivity model
!
!> @brief
!!   Contains data for the SURFEM-Ocean MW sea surface emissivity model
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
!    Copyright 2022, EUMETSAT, All Rights Reserved.

MODULE rttov_surfem_ocean_coef_mod

  USE parkind1, ONLY : jprb, jpim
  USE rttov_const, ONLY : t0

  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  PUBLIC :: surfem_iso_net, surfem_aniso_net

  ! Dielectric parameterisation input variable limits
  REAL(jprb), PUBLIC, PARAMETER :: surfem_dielec_min_sst   = t0 - 40._jprb ! K, -40 deg C
  REAL(jprb), PUBLIC, PARAMETER :: surfem_dielec_max_sst   = t0 + 34._jprb ! K,  34 deg C
  REAL(jprb), PUBLIC, PARAMETER :: surfem_dielec_min_sss   = 0._jprb
  REAL(jprb), PUBLIC, PARAMETER :: surfem_dielec_max_sss   = 40._jprb

  ! Neural network training limits
  ! NB zenith angle limits are 0-89, but RTTOV will not allow anything outside this range
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_min_sst       = t0 - 2._jprb  ! K, -2 deg C
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_max_sst       = t0 + 32._jprb ! K, 32 deg C
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_min_sss       = 0._jprb
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_max_sss       = 40._jprb
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_min_windspeed = 0._jprb
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_max_windspeed = 50._jprb
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_min_freq      = 0.5_jprb
  REAL(jprb), PUBLIC, PARAMETER :: surfem_nn_max_freq      = 700._jprb

  ! (1/(2*pi*epsi_0) term in conductivity term in GHz.m/S
  REAL(jprb), PUBLIC, PARAMETER :: f0 = 17.97510_jprb

  !> es:  Static dielectric constant for pure water by Stogryn et al. 1995
  REAL(jprb), PUBLIC, PARAMETER :: es_coef(3) = (/ &
   3.70886E4_jprb   , 8.2168E1_jprb    , 4.21854E2_jprb  /) 

  !> ai coefficients from Table III in Meissner and Wentz 2004 
  REAL(jprb), PUBLIC, PARAMETER :: ai_coef(11) = (/ &
   5.7230_jprb      , 2.2379E-02_jprb  ,-7.1237E-04_jprb  , 5.0478_jprb      , -7.0315E-02_jprb , &
   6.0059E-04_jprb  , 3.6143_jprb      , 2.8841E-02_jprb  , 1.3652E-01_jprb  ,  1.4825E-03_jprb , &
   2.4166E-04_jprb /) 

  !> sigma35:   Salt water. Conductivity [s/m] by Stogryn et al. 1995
  REAL(jprb), PUBLIC, PARAMETER :: sigma35_coef(5) = (/ &
   2.903602_jprb    , 8.60700E-02_jprb , 4.738817E-04_jprb, 2.9910E-06_jprb  ,  4.3047E-09_jprb /) 
  !> R15:
  REAL(jprb), PUBLIC, PARAMETER :: R15_coef(5) = (/ &
   37.5109_jprb     , 5.45216_jprb     , 1.4409E-02_jprb   , 1004.75_jprb     ,  182.283_jprb    /) 
  !> alpha0:
  REAL(jprb), PUBLIC, PARAMETER :: alpha0_coef(5) = (/ &
   6.9431_jprb      , 3.2841_jprb      , 9.9486E-02_jprb   , 84.850_jprb      ,  69.024_jprb     /)
  !> alpha1:
  REAL(jprb), PUBLIC, PARAMETER :: alpha1_coef(3) = (/ &
   49.843_jprb     , 0.2276_jprb       , 0.198E-02_jprb  /)

  !> bi coefficients from Table VI in Meissner and Wentz 2004
  REAL(jprb), PUBLIC, PARAMETER :: bi_coef(13) = (/ &
 -3.56417E-03_jprb ,  4.74868E-06_jprb , 1.15574E-05_jprb , 2.39357E-03_jprb , -3.13530E-05_jprb, &
  2.52477E-07_jprb , -6.28908E-03_jprb , 1.76032E-04_jprb ,-9.22144E-05_jprb , -1.99723E-02_jprb, &
  1.81176E-04_jprb , -2.04265E-03_jprb , 1.57883E-04_jprb /)

  !> updated bi (i = 0--2) coefficients from Table VI in Meissner and Wentz 2012
  REAL(jprb), PUBLIC, PARAMETER :: bi_new1_coef(3) = (/ &
 -0.33330E-02_jprb , 4.74868E-06_jprb  , 0.0_jprb /)

  !> updated bi (i = 3--5) coefficients from Table VII in Meissner and Wentz 2012 
  !> (now uses 5 coefficients instead of 3), includes change of sign on coef number 4 
  !> reported in Meissner et al. 2016
  REAL(jprb), PUBLIC, PARAMETER :: bi_new2_coef(5) = (/ &
  0.23232E-02_jprb ,-0.79208E-04_jprb  , 0.36764E-05_jprb ,-0.35594E-06_jprb , 0.89795E-08_jprb /)

  !>Salt water parameters
  REAL(jprb), PUBLIC, PARAMETER :: c1_coef(2) = (/ &
  9.1873715E-04_jprb , 1.5012396E-04_jprb /)


!===========rough sea coefficients====================

    INTEGER(jpim), PUBLIC, PARAMETER :: surfem_ni = 7       !< Number of input variables (freq, theta, windspeed, 
                                                            !!   tskin, salinity, e_Nh, e_Nv)
    INTEGER(jpim), PUBLIC, PARAMETER :: surfem_nh = 180     !< Number of nodes in the hidden layer 
    INTEGER(jpim), PUBLIC, PARAMETER :: surfem_no_iso = 2   !< Number of net outputs for isotropic network
    INTEGER(jpim), PUBLIC, PARAMETER :: surfem_no_aniso = 8 !< Number of net outputs for anisotropic network

!> Structure holding isotropic neural network data
  TYPE surfem_iso_net
    REAL(jprb)  :: b1(surfem_nh)
    REAL(jprb)  :: W1(surfem_ni,surfem_nh)
    REAL(jprb)  :: b2(surfem_no_iso)
    REAL(jprb)  :: W2(surfem_nh,surfem_no_iso)
    REAL(jprb)  :: ti_gain(surfem_ni)
    REAL(jprb)  :: ti_xoffset(surfem_ni)
    REAL(jprb)  :: ti_ymin(1)
    REAL(jprb)  :: to_gain(surfem_no_iso)
    REAL(jprb)  :: to_xoffset(surfem_no_iso)
    REAL(jprb)  :: to_ymin(1)
  END TYPE surfem_iso_net

!> Structure holding anisotropic neural network data 
!> only difference is surfem_no_aniso is different to surfem_no_iso
  TYPE surfem_aniso_net
    REAL(jprb)  :: b1(surfem_nh)
    REAL(jprb)  :: W1(surfem_ni,surfem_nh)
    REAL(jprb)  :: b2(surfem_no_aniso)
    REAL(jprb)  :: W2(surfem_nh,surfem_no_aniso)
    REAL(jprb)  :: ti_gain(surfem_ni)
    REAL(jprb)  :: ti_xoffset(surfem_ni)
    REAL(jprb)  :: ti_ymin(1)
    REAL(jprb)  :: to_gain(surfem_no_aniso)
    REAL(jprb)  :: to_xoffset(surfem_no_aniso)
    REAL(jprb)  :: to_ymin(1)
  END TYPE surfem_aniso_net

!> Neural network 
TYPE(surfem_iso_net), PUBLIC, TARGET :: iso_net

! bias vector of first layer, vector [nh * 1]
! 180 elements in b1
DATA iso_net%b1 /    -5.5481738E+00_jprb,  3.9259262E+00_jprb,  1.2903522E-01_jprb,  5.7747632E+00_jprb, -3.4955900E+00_jprb, &
 1.9731688E+00_jprb,  2.6186546E+00_jprb,  9.6486797E-01_jprb,  6.4886240E+00_jprb,  1.7422257E+01_jprb, -6.2273779E-01_jprb, &
-5.6164188E+00_jprb,  8.8147782E-02_jprb,  1.9090410E+00_jprb, -3.8945939E+00_jprb, -2.5731869E+00_jprb, -3.4574890E+00_jprb, &
 2.9375877E+00_jprb, -2.7923329E+00_jprb, -5.9966955E-02_jprb, -2.5470903E-01_jprb, -3.1523754E+00_jprb,  1.6734683E+00_jprb, &
-1.6850978E+00_jprb,  4.9600477E-01_jprb,  7.1460770E-01_jprb, -1.6503954E+00_jprb,  1.7002254E+00_jprb,  6.6826480E+00_jprb, &
-3.3066180E-01_jprb, -2.3343158E+00_jprb, -2.8836894E+00_jprb, -2.8866139E+00_jprb, -7.5490631E-01_jprb,  1.2374446E+00_jprb, &
 1.9656897E+00_jprb, -4.2194173E+00_jprb,  4.6965353E-01_jprb, -2.5322403E-01_jprb, -1.3409111E+00_jprb,  1.4733955E-01_jprb, &
-6.4585518E+00_jprb,  7.7035274E-01_jprb, -3.4051257E+00_jprb,  8.4632954E-01_jprb, -2.4466905E+00_jprb,  8.7709190E-02_jprb, &
 8.2858345E-01_jprb,  2.6014941E+00_jprb,  5.7178763E-01_jprb,  9.5499461E-01_jprb, -3.5066699E+00_jprb,  5.7141639E+00_jprb, &
-1.2282062E+00_jprb, -7.0709510E-02_jprb, -1.0320635E+00_jprb,  3.6946125E+00_jprb, -5.1425969E+00_jprb,  5.8039958E+00_jprb, &
 9.6453709E-01_jprb, -8.7771914E-02_jprb,  4.4305016E-01_jprb, -3.9167391E+00_jprb,  5.3693636E-01_jprb,  3.9559449E-01_jprb, &
-3.0413259E+00_jprb,  1.5686924E+00_jprb,  1.4590518E+00_jprb, -2.5007534E-01_jprb, -5.1717349E-01_jprb, -9.4267242E-03_jprb, &
 3.1883431E+00_jprb, -5.1271489E+00_jprb,  1.7351447E+00_jprb, -4.3630587E+00_jprb,  6.2340466E+00_jprb, -1.2541247E+00_jprb, &
 1.4336153E+00_jprb,  3.6713399E-01_jprb,  4.5302409E+00_jprb,  6.6148257E+00_jprb,  1.9939988E+00_jprb,  1.9188821E+00_jprb, &
 9.7495215E-02_jprb,  3.0831234E+00_jprb, -6.7468324E-01_jprb,  3.4251731E-01_jprb, -5.9487855E+00_jprb,  1.0563463E+00_jprb, &
 8.9569879E+00_jprb,  1.8077609E+00_jprb,  6.3080465E-01_jprb, -4.1943947E-01_jprb,  3.2394749E+00_jprb, -2.9533144E-01_jprb, &
-1.0049421E+00_jprb, -5.7914343E+00_jprb, -2.1820610E+00_jprb,  7.2663183E-01_jprb, -2.5277851E+00_jprb, -3.6921785E+00_jprb, &
 1.7098985E+00_jprb,  8.6796439E-01_jprb,  3.4847648E+00_jprb,  1.1885506E-01_jprb, -2.4258304E+00_jprb,  7.8384750E-01_jprb, &
-1.5209127E+00_jprb, -2.0123938E+00_jprb,  8.7896217E+00_jprb,  4.0862815E+00_jprb,  1.2694204E+01_jprb,  5.6236062E+00_jprb, &
-1.4390538E-01_jprb,  6.1379520E-01_jprb, -8.3365086E+00_jprb,  1.7487500E+00_jprb, -6.4238095E-01_jprb, -9.0614889E-01_jprb, &
-1.0806230E+00_jprb, -1.5049265E-02_jprb, -5.2789586E-02_jprb,  2.0301096E+00_jprb,  3.7255898E+00_jprb,  3.1291891E+00_jprb, &
-7.6105773E-02_jprb, -6.6583291E-01_jprb,  1.0344484E+01_jprb,  1.2671066E-01_jprb,  1.6368579E+00_jprb,  2.1274504E+00_jprb, &
-4.3137736E+00_jprb,  4.1388680E+00_jprb,  2.6691611E+00_jprb,  2.3639628E+00_jprb,  2.4324766E-01_jprb, -4.8662650E+00_jprb, &
-1.7035996E+01_jprb, -8.1819474E-01_jprb,  1.0406554E+00_jprb,  1.0212438E+00_jprb, -5.4983094E+00_jprb, -1.2130910E+00_jprb, &
-9.9552048E+00_jprb,  2.5364106E+00_jprb, -5.4081412E+00_jprb,  1.4700785E+01_jprb,  2.2844158E+00_jprb, -2.2517876E+00_jprb, &
 6.3110574E+00_jprb,  1.8515557E+00_jprb,  4.3082534E+00_jprb,  9.8313725E+00_jprb, -1.1260268E+00_jprb, -3.2519601E+00_jprb, &
-3.0361690E+00_jprb, -2.7540413E-02_jprb, -9.5756237E+00_jprb, -6.6322780E+00_jprb, -9.8585120E-01_jprb,  1.1244647E+01_jprb, &
 8.4663933E+00_jprb, -7.9238526E+00_jprb, -2.3590854E+00_jprb, -3.3792361E+00_jprb, -3.8287921E-01_jprb, -3.9592118E+00_jprb, &
 4.3456090E+00_jprb, -9.1350570E-01_jprb,  1.5412949E+01_jprb, -1.9441747E+00_jprb, -2.4888740E+00_jprb,  1.7883404E+01_jprb, &
 8.5767523E+00_jprb, -1.5459566E+00_jprb,  2.0453268E+00_jprb, -3.4888357E+00_jprb, -1.5259789E+00_jprb, -5.5527006E+00_jprb, &
-1.2275883E+01_jprb /     
   

! weight matrix of first layer, matrix [nh (columns) * ni (rows)]
! 180 * 7 elements in W1 
DATA iso_net%W1 /  &     
 -9.2979879E-01_jprb,  4.4571587E-01_jprb, -1.7062251E+00_jprb, -8.3148321E-03_jprb,  8.8258917E-03_jprb,  6.2399664E-01_jprb, &
 -8.2878242E-01_jprb,  4.1461704E+00_jprb,  2.4964055E-01_jprb, -4.4311678E-01_jprb, -3.9028611E-02_jprb,  1.9626402E-02_jprb, &
 -1.4928297E-01_jprb, -1.7351189E-01_jprb,  6.7916759E-02_jprb, -9.9813995E-02_jprb, -4.5443032E-01_jprb,  2.7365763E-02_jprb, &
  3.9750499E-03_jprb,  2.3763076E-01_jprb,  8.2098597E-01_jprb, -1.4722619E-02_jprb, -1.3479225E+00_jprb,  4.8354880E+00_jprb, &
  8.7738758E-03_jprb, -4.2596351E-03_jprb,  2.2636867E-01_jprb, -4.2103878E-01_jprb,  8.4057173E-02_jprb, -4.7387358E-01_jprb, &
  2.7071709E-01_jprb, -2.9345094E-02_jprb,  1.3307743E-02_jprb, -3.5555809E-01_jprb, -4.4159311E+00_jprb,  1.4495187E+00_jprb, &
 -7.5977867E-01_jprb,  1.0110700E-01_jprb,  2.3144523E-02_jprb,  8.8411296E-03_jprb, -3.3403279E-01_jprb,  2.4711803E-01_jprb, &
 -5.1312804E-02_jprb, -1.1151146E+00_jprb,  2.4704847E+00_jprb, -6.3171809E-03_jprb, -1.5071587E-04_jprb,  4.2148163E-02_jprb, &
 -3.5330817E-01_jprb,  5.7952974E-01_jprb, -5.6962362E-01_jprb,  1.5617670E-01_jprb, -1.7928660E-02_jprb, -2.6088910E-03_jprb, &
 -1.0404523E-01_jprb,  3.0656032E-03_jprb,  5.0206411E-02_jprb,  2.8451973E-01_jprb,  2.0794085E+00_jprb, -1.1312285E-02_jprb, &
 -6.2077850E-03_jprb, -8.3137033E-02_jprb,  5.5121670E+00_jprb,  8.3190051E-02_jprb, -6.2685587E+00_jprb,  1.2066237E+01_jprb, &
 -2.9563476E-03_jprb, -1.2024903E-03_jprb,  4.3895227E-01_jprb, -1.0780342E+00_jprb,  1.3384641E-03_jprb, -2.8723239E-01_jprb, &
  6.8066845E-01_jprb, -7.9047835E-03_jprb,  4.7113356E-04_jprb, -9.2635657E-02_jprb, -1.1558706E+00_jprb, -5.2302845E+00_jprb, &
  5.3930237E-01_jprb,  1.9852186E-01_jprb, -1.4673135E-02_jprb, -1.9877275E-02_jprb,  7.6745861E-02_jprb, -1.7787875E-01_jprb, &
 -3.5355859E-02_jprb, -2.4744932E-01_jprb,  8.4479142E-01_jprb, -1.1038183E-02_jprb,  1.0534527E-03_jprb, -7.6277296E-02_jprb, &
 -6.5378632E-01_jprb,  1.9358369E+00_jprb,  5.0596129E-02_jprb, -5.7432370E-01_jprb, -5.3975200E-02_jprb, -3.8812749E-02_jprb, &
 -9.9739823E-01_jprb, -8.0326549E-01_jprb, -2.4383054E+00_jprb, -1.5206303E-01_jprb, -4.3148399E-01_jprb, -9.8161160E-02_jprb, &
 -2.7414470E-02_jprb, -1.3494500E+00_jprb, -2.1084344E+00_jprb, -1.1116440E+00_jprb,  3.0285878E-01_jprb,  3.8813813E-01_jprb, &
 -9.6929317E-02_jprb, -1.3996950E-02_jprb, -1.6737284E-01_jprb,  3.6000308E-01_jprb, -2.2898453E-01_jprb,  3.9213364E-01_jprb, &
 -1.9554059E+00_jprb,  2.1531371E-02_jprb, -1.8965496E-03_jprb, -8.9125974E-02_jprb, -7.1204184E-01_jprb,  3.2346468E+00_jprb, &
 -2.3472668E-01_jprb, -6.7730977E-01_jprb,  2.3401243E-02_jprb,  2.2066186E-02_jprb, -4.5424809E-02_jprb, -5.7101693E-04_jprb, &
  8.0897244E-02_jprb, -8.6784862E-02_jprb,  3.5464727E-01_jprb,  2.7974663E-03_jprb,  4.5356615E-03_jprb, -1.4308270E-02_jprb, &
 -3.7655189E+00_jprb, -2.5835578E-01_jprb, -1.5092543E-01_jprb, -8.9982824E-02_jprb,  6.3071881E-02_jprb, -4.3308953E-02_jprb, &
  5.0784082E-01_jprb, -2.5622071E-01_jprb,  1.9499358E-03_jprb,  1.7483992E-01_jprb, -2.5435567E+00_jprb, -1.3391246E-02_jprb, &
 -1.9420712E-03_jprb, -1.7160377E-01_jprb,  2.6361945E-01_jprb, -3.3691936E-02_jprb,  8.7087075E-01_jprb, -2.7508919E+00_jprb, &
  1.5650491E-02_jprb,  2.3903542E-03_jprb, -4.1324623E-01_jprb,  8.5276068E-01_jprb,  1.0516539E+00_jprb, -3.3535385E-01_jprb, &
  1.0954073E-01_jprb,  8.0777310E-02_jprb,  2.5212141E-02_jprb,  8.4196448E-01_jprb,  1.5010275E+00_jprb,  4.8585223E-01_jprb, &
  7.5301514E-01_jprb,  2.9529383E-01_jprb,  4.2477442E-02_jprb,  2.7668819E-02_jprb, -1.9328132E-01_jprb, -2.5883452E-01_jprb, &
  3.5496171E-02_jprb, -2.8820117E-01_jprb,  7.5778944E-01_jprb, -1.9802679E-03_jprb, -4.6233189E-04_jprb,  1.0419740E-01_jprb, &
 -1.2325865E+00_jprb,  1.6666296E-01_jprb,  2.1438626E-01_jprb, -1.9662614E-01_jprb,  1.5964113E-01_jprb,  5.5616525E-02_jprb, &
 -4.6228809E-02_jprb,  1.0292634E+00_jprb,  4.6557264E-03_jprb,  4.4933589E-01_jprb, -1.4034206E+00_jprb, -4.9456415E-04_jprb, &
  1.0745853E-03_jprb,  3.9598937E-01_jprb,  5.1946898E-01_jprb,  2.1976360E+00_jprb,  3.9240255E-01_jprb, -7.8704339E-01_jprb, &
  4.9924572E-02_jprb,  2.0881443E-02_jprb,  3.1885411E-02_jprb,  3.3566105E-01_jprb,  2.2376790E+00_jprb,  1.0192159E+00_jprb, &
  5.7181550E-02_jprb,  1.6352632E-02_jprb,  7.0687888E-02_jprb,  3.9149583E+00_jprb,  9.7118252E+00_jprb,  4.8293830E-01_jprb, &
 -9.3600977E-01_jprb, -6.9094258E-01_jprb, -1.0787445E-01_jprb,  1.8937550E-02_jprb, -2.1663415E-01_jprb, -1.4295595E+00_jprb, &
 -3.1099939E-01_jprb, -1.8654057E-02_jprb, -1.3377822E-01_jprb, -2.4743606E-02_jprb,  1.4474119E-03_jprb, -4.2031534E-01_jprb, &
 -2.3298025E+00_jprb, -2.9826216E+00_jprb,  3.3107815E-02_jprb,  6.5132501E-01_jprb, -5.1206533E-02_jprb, -2.1932144E-02_jprb, &
  3.4169608E-02_jprb, -1.8734167E-01_jprb, -2.1426634E-02_jprb, -3.7918600E-01_jprb,  4.2002889E-01_jprb, -2.4305402E-03_jprb, &
  2.9280851E-03_jprb, -4.7285640E-03_jprb, -3.1390404E+00_jprb, -4.3705507E-01_jprb,  5.1661650E-01_jprb, -6.5704158E-02_jprb, &
  9.6857810E-02_jprb,  9.4517490E-02_jprb, -4.5606367E-01_jprb,  2.3533582E-01_jprb, -2.2897180E-02_jprb, -5.6534599E-01_jprb, &
  2.6546113E-01_jprb,  2.7476675E-02_jprb,  1.4060951E-03_jprb,  1.1378572E+00_jprb, -1.4845594E+00_jprb,  1.9457110E+00_jprb, &
  4.0806453E-02_jprb, -5.7859371E-01_jprb, -5.8536964E-02_jprb, -3.8473931E-02_jprb, -1.0524793E+00_jprb, -8.2195268E-01_jprb, &
 -3.0121716E+00_jprb, -9.5233990E-02_jprb, -9.9346172E-01_jprb,  1.8606577E-02_jprb, -1.4377912E-02_jprb, -4.3904801E-02_jprb, &
 -2.4947284E-01_jprb, -7.0001092E-02_jprb,  4.1602810E-01_jprb,  3.7603031E-01_jprb, -1.4220068E-02_jprb, -3.0816121E-03_jprb, &
 -5.4705729E-01_jprb,  1.8739568E-01_jprb, -1.4947763E-01_jprb, -7.5608312E-01_jprb, -4.4073418E-02_jprb, -8.1712325E-03_jprb, &
 -2.9563576E-03_jprb, -3.2785903E-01_jprb, -7.1430632E-03_jprb, -8.2299071E-01_jprb, -2.5717835E-03_jprb,  3.2371193E-02_jprb, &
 -1.0848985E-01_jprb, -2.9422600E-02_jprb, -1.0551377E+00_jprb, -1.6372524E+00_jprb,  3.0196511E-01_jprb,  6.9227767E-01_jprb, &
  2.0394639E-01_jprb, -3.1746292E-02_jprb, -1.1676353E-02_jprb, -2.9033413E+00_jprb, -2.5570774E-01_jprb, -5.7096230E+00_jprb, &
  4.2305987E-01_jprb, -1.5475364E-01_jprb,  1.6798078E-02_jprb, -1.4993857E-02_jprb,  6.4242757E-02_jprb, -1.5786027E-01_jprb, &
  6.4541847E-01_jprb, -4.2286948E-01_jprb,  4.8489840E-02_jprb,  1.5806858E-02_jprb, -1.4042590E-02_jprb,  3.2794425E-02_jprb, &
  4.1871133E-02_jprb, -1.5185166E+00_jprb, -2.2120500E-01_jprb,  1.0565461E+00_jprb,  4.5187863E-02_jprb, -9.6709381E-02_jprb, &
 -2.6591743E+00_jprb, -5.3377191E+00_jprb,  2.9818048E-01_jprb,  6.5599532E-01_jprb, -7.3518215E-03_jprb, -8.4895008E-02_jprb, &
 -2.5191490E-03_jprb, -8.0881528E-01_jprb,  5.9587761E-01_jprb, -6.4435723E-02_jprb,  6.2299747E-01_jprb, -1.9091091E+00_jprb, &
  9.2290035E-04_jprb,  1.8096769E-03_jprb,  1.2153726E-01_jprb, -9.5780758E-02_jprb, -2.0464116E-01_jprb, -2.0007676E-01_jprb, &
  5.0562330E-02_jprb,  4.2543679E-02_jprb, -6.5675256E-03_jprb, -1.5208592E-01_jprb, -8.3026800E-01_jprb, -4.2637154E-02_jprb, &
 -1.7761359E-02_jprb, -2.9713148E-01_jprb, -2.0645428E-02_jprb, -5.4515034E-03_jprb, -9.9346491E-01_jprb,  7.9045873E-01_jprb, &
  1.0280711E+00_jprb, -1.0229325E+00_jprb, -3.9103399E-01_jprb,  9.2261227E-02_jprb,  7.4560034E-02_jprb,  3.8579725E+00_jprb, &
  5.5211981E+00_jprb,  2.3637511E-02_jprb,  1.4871169E-01_jprb,  3.4777176E-01_jprb,  4.7626590E-03_jprb, -3.5521106E-03_jprb, &
 -1.3568540E-01_jprb,  1.1687029E+00_jprb,  5.4443009E-01_jprb, -2.4545695E-01_jprb,  1.1245371E-01_jprb, -4.5982947E-02_jprb, &
 -8.8125449E-02_jprb,  2.4591063E-01_jprb,  1.3175241E-01_jprb,  2.0773224E-02_jprb, -1.8521673E-01_jprb, -1.1527595E+00_jprb, &
  5.6313631E-03_jprb,  4.2924693E-03_jprb,  1.2302489E-01_jprb, -3.3533704E+00_jprb,  4.7746255E+00_jprb, -3.0113040E-02_jprb, &
  1.0724819E+00_jprb,  8.7356767E-03_jprb,  2.5642492E-02_jprb,  6.6632662E-02_jprb,  1.7953296E-01_jprb, -1.1092410E+00_jprb, &
 -2.7215496E-02_jprb,  1.8294548E-01_jprb,  1.4429474E-03_jprb,  6.4980061E-03_jprb, -5.3948577E-01_jprb,  1.7014427E-01_jprb, &
 -1.0374110E-01_jprb, -4.2921459E-02_jprb, -1.1617897E-01_jprb,  1.0943411E-02_jprb,  2.4512365E-04_jprb,  4.2070047E-01_jprb, &
 -9.7117565E-01_jprb,  1.0182716E-02_jprb,  5.5804124E-01_jprb,  4.5140839E-01_jprb,  3.2470199E-02_jprb,  9.3234071E-03_jprb, &
  9.7838574E-01_jprb, -3.9948599E-01_jprb, -1.4606294E-01_jprb, -2.1172088E+00_jprb, -1.5195805E-02_jprb,  3.9726073E-03_jprb, &
  4.7688248E-04_jprb,  8.2993992E-02_jprb,  1.9979563E+00_jprb,  5.8788527E-02_jprb,  1.4300600E+00_jprb, -4.2042164E+00_jprb, &
 -9.2430923E-04_jprb,  3.0364616E-03_jprb,  3.3918820E-02_jprb,  1.5358186E-01_jprb,  5.4441250E+00_jprb,  1.8000032E-01_jprb, &
  5.5105031E-01_jprb,  8.5757623E-03_jprb,  1.8400092E-02_jprb,  3.2333442E-02_jprb,  2.6270002E-01_jprb, -5.2908915E-02_jprb, &
 -3.6401674E-02_jprb,  6.1747335E-01_jprb,  3.9771111E-04_jprb, -1.2030089E-03_jprb, -1.1917454E-01_jprb, -5.0649838E-01_jprb, &
  5.1652739E-02_jprb, -2.5015322E-01_jprb,  6.6440742E-01_jprb,  1.5296889E-03_jprb,  2.5891434E-03_jprb,  8.1299763E-02_jprb, &
 -1.7084275E+00_jprb,  7.9470727E-01_jprb,  6.0773983E-01_jprb, -4.9363612E-01_jprb, -1.1723786E-01_jprb,  4.9607568E-02_jprb, &
 -1.2793774E+00_jprb, -1.3306503E+00_jprb,  1.0390464E-01_jprb,  2.6731256E-01_jprb, -2.3160785E+00_jprb, -2.0838671E-02_jprb, &
  5.5952614E-03_jprb,  2.9460316E-01_jprb, -1.1575763E-01_jprb, -2.6303193E-01_jprb, -1.1050511E+00_jprb, -6.2795343E-02_jprb, &
  1.1644727E-03_jprb, -5.9895860E-03_jprb,  9.2784891E-01_jprb, -1.3624272E+00_jprb,  2.9974855E-01_jprb,  9.3754176E-01_jprb, &
 -1.1596077E+00_jprb,  5.6377298E-03_jprb,  2.2665416E-03_jprb, -3.5853236E-01_jprb,  5.8272201E-01_jprb, -2.4490651E-02_jprb, &
  8.9316297E-01_jprb, -2.8773242E+00_jprb,  3.7788737E-03_jprb,  1.6113553E-03_jprb, -2.8475986E-01_jprb,  5.6880343E-01_jprb, &
  1.0340627E+00_jprb, -2.0797560E-01_jprb,  5.6435599E-02_jprb,  9.7616174E-02_jprb,  3.0595238E-02_jprb,  8.8429854E-01_jprb, &
  1.5972117E+00_jprb, -4.2689664E-03_jprb, -7.9457931E-01_jprb,  1.5695643E+00_jprb, -4.8436964E-04_jprb, -6.6750715E-04_jprb, &
 -1.9851849E-01_jprb, -5.1283713E-01_jprb, -2.9819096E-02_jprb,  5.2947095E-01_jprb,  5.3530538E-01_jprb, -2.4821000E-03_jprb, &
  9.8354304E-04_jprb, -8.4915489E-02_jprb, -7.5824008E-01_jprb, -1.9320861E-02_jprb,  7.7387820E-01_jprb, -1.2523638E+00_jprb, &
  5.8457695E-03_jprb, -2.1181860E-04_jprb,  1.8870666E-02_jprb,  6.7158743E-01_jprb,  8.3464990E-02_jprb,  4.7357035E-01_jprb, &
 -2.2917511E-02_jprb,  1.4636141E-02_jprb, -8.9632346E-04_jprb,  8.6534858E-01_jprb, -7.6886546E-01_jprb,  2.9149816E+00_jprb, &
  1.7365248E-01_jprb,  9.6564375E-01_jprb,  3.2821645E-01_jprb,  8.5883631E-02_jprb,  9.3572728E-02_jprb,  4.6459762E-01_jprb, &
  7.4685093E-03_jprb,  9.8836857E-01_jprb, -3.4214659E+00_jprb,  5.1091374E-03_jprb,  5.7598939E-03_jprb,  8.4117868E-01_jprb, &
  4.2585451E-02_jprb,  1.4738944E+00_jprb,  6.4405319E-02_jprb, -1.9416046E-01_jprb,  1.5502632E-01_jprb,  4.2862797E-02_jprb, &
  5.4174292E-01_jprb,  6.6159988E-01_jprb, -4.7496584E+00_jprb, -3.3934944E-01_jprb,  2.9499805E-01_jprb, -7.3731487E-03_jprb, &
  3.4259282E-02_jprb,  1.6457787E-01_jprb,  3.9720553E-04_jprb,  2.1871525E-01_jprb, -1.8518453E+00_jprb,  5.6368414E+00_jprb, &
 -2.5667487E-02_jprb, -2.8239602E-05_jprb,  1.2701574E+00_jprb, -6.7686608E-01_jprb, -1.8317044E-01_jprb, -8.9662036E-02_jprb, &
 -6.5550811E-02_jprb,  4.1773595E-02_jprb, -2.1952845E-03_jprb,  1.2168807E+00_jprb, -3.7250253E-01_jprb, -2.6904104E-02_jprb, &
  2.3497769E-01_jprb,  1.2329240E+00_jprb,  3.3855034E-02_jprb, -9.3408806E-03_jprb, -1.0031460E+00_jprb, -2.5263504E-01_jprb, &
 -2.4194940E-01_jprb,  1.3149530E-01_jprb,  3.9417321E-01_jprb,  7.8823933E-02_jprb, -2.3309398E-03_jprb,  1.6480540E-01_jprb, &
  2.0104613E-01_jprb, -7.1768815E-02_jprb, -3.7008693E+00_jprb,  3.9778344E-01_jprb,  4.6561350E-03_jprb,  2.6308628E-03_jprb, &
  6.5948234E-01_jprb,  6.1144307E-01_jprb,  2.2215043E+00_jprb,  1.0054307E+00_jprb,  3.4702999E-02_jprb,  1.6652618E-02_jprb, &
  7.2188187E-02_jprb,  3.8852323E+00_jprb,  9.6230482E+00_jprb,  2.2514366E+00_jprb,  2.4200557E-01_jprb, -1.9853490E-01_jprb, &
  1.5928052E-01_jprb,  2.3456124E-02_jprb,  3.8848921E-01_jprb,  3.6236636E-01_jprb, -1.7117008E-01_jprb, -2.1801923E+00_jprb, &
  3.4307230E-02_jprb,  7.8118381E-03_jprb, -2.3502325E-03_jprb,  2.3473237E-01_jprb,  8.3641878E-01_jprb, -4.5338604E-01_jprb, &
 -6.0928010E-01_jprb, -1.3964309E-01_jprb, -2.2895488E-01_jprb, -1.3436763E-01_jprb,  5.4702356E-01_jprb, -5.4748473E-01_jprb, &
  2.0483505E+00_jprb, -5.7822858E-01_jprb, -7.2169964E-03_jprb,  4.6500021E-02_jprb,  9.7677687E-03_jprb,  1.0243914E+00_jprb, &
  1.4932022E+00_jprb, -4.2693010E-01_jprb, -8.7946060E-01_jprb, -1.9584795E-03_jprb,  3.5110213E-02_jprb, -6.0163396E-03_jprb, &
  9.0905557E-01_jprb, -8.0336854E-01_jprb, -3.0183107E-01_jprb, -6.4924766E-01_jprb, -1.5097909E-01_jprb, -3.1918935E-01_jprb, &
 -1.6682236E-01_jprb,  3.3265715E-01_jprb, -7.6006474E-01_jprb, -4.2331599E+00_jprb,  2.1714556E+00_jprb,  9.9230030E-02_jprb, &
  5.2796620E-02_jprb,  2.0401657E-03_jprb, -2.0898005E-02_jprb,  6.3001552E-01_jprb,  1.5159840E-01_jprb, -8.4283895E-01_jprb, &
 -2.4776936E-01_jprb,  3.2819751E-02_jprb, -1.2658812E-03_jprb,  1.0980257E+00_jprb,  1.6805004E+00_jprb,  8.8947965E+00_jprb, &
  2.5939988E-01_jprb, -7.2041861E-01_jprb,  7.2769131E-02_jprb,  1.9150598E-02_jprb, -8.9783678E-02_jprb,  2.9837979E-01_jprb, &
  5.1621308E-01_jprb, -1.2434515E-01_jprb,  2.2147832E-01_jprb,  2.2527534E-02_jprb,  1.2796891E-02_jprb,  9.1235528E-01_jprb, &
  1.5002131E+00_jprb, -1.7300941E-02_jprb, -5.5034623E-01_jprb,  1.3570241E+00_jprb,  2.2887989E-03_jprb, -4.9925923E-04_jprb, &
 -1.5706310E-02_jprb, -7.0882810E-01_jprb, -2.7220385E-01_jprb, -9.3059844E-01_jprb,  1.1897795E+00_jprb, -3.5110995E-03_jprb, &
 -2.8838096E-03_jprb,  3.2593281E-01_jprb, -5.5435555E-01_jprb,  8.7382172E-01_jprb,  1.4812502E-01_jprb, -1.2887409E-01_jprb, &
  1.7814444E-02_jprb,  1.7870056E-02_jprb,  2.2048542E+00_jprb,  3.3338001E+00_jprb, -1.5122702E-01_jprb,  2.9480328E-02_jprb, &
 -2.6435499E-01_jprb,  1.2582424E-02_jprb, -5.4278491E-02_jprb,  3.0179754E-02_jprb, -3.3596852E-01_jprb, -7.3901977E-03_jprb, &
  1.7742340E+00_jprb,  1.0546910E-01_jprb,  1.0318889E-02_jprb, -2.3597687E-03_jprb, -7.3774279E-01_jprb,  2.8350842E-02_jprb, &
 -5.1968914E+00_jprb, -2.2072615E-01_jprb, -6.1835098E-01_jprb, -8.8215870E-04_jprb, -1.8393217E-02_jprb, -6.3641053E-02_jprb, &
 -1.9807606E-01_jprb, -1.0350136E-01_jprb, -2.7094011E-01_jprb,  8.2072819E-01_jprb, -9.6174697E-03_jprb,  1.3744206E-03_jprb, &
 -9.8715825E-02_jprb, -2.7115954E+00_jprb, -9.9020621E-02_jprb, -1.0322486E+00_jprb, -2.4274821E-02_jprb, -6.4753824E-03_jprb, &
 -1.5082613E-03_jprb, -1.9589470E-01_jprb,  6.6195547E-01_jprb, -2.0361767E+00_jprb,  6.2064825E-01_jprb,  1.0239006E-01_jprb, &
 -4.4594669E-02_jprb,  2.0435296E-02_jprb, -5.5920129E-01_jprb,  8.0013844E-02_jprb, -2.3928981E-01_jprb,  4.2506111E+00_jprb, &
  8.4714030E-02_jprb, -4.0794880E-02_jprb,  3.0417362E-03_jprb, -2.3060713E+00_jprb, -4.1114497E-01_jprb,  8.9762861E-04_jprb, &
  1.2687510E-01_jprb,  6.4210203E-01_jprb, -6.5214566E-03_jprb, -2.1601007E-03_jprb, -1.7280536E-01_jprb,  1.9609870E+00_jprb, &
  2.7529535E-01_jprb, -3.6170388E-01_jprb,  7.2112786E-03_jprb,  2.6703889E-02_jprb,  3.9133304E-03_jprb,  1.7390417E+00_jprb, &
 -3.8636709E-01_jprb,  2.1448671E+00_jprb,  2.4970981E-01_jprb,  3.8689279E-01_jprb,  8.3547751E-02_jprb,  2.6255448E-02_jprb, &
  1.1292990E+00_jprb,  2.1076419E+00_jprb, -1.2585349E-01_jprb,  1.2234767E-01_jprb,  7.6365493E-01_jprb, -8.1991312E-03_jprb, &
  1.2017672E-03_jprb, -1.4589281E-01_jprb, -4.2092569E-01_jprb,  1.0392508E-01_jprb, -8.9524945E-02_jprb,  3.4278838E-01_jprb, &
  6.5913294E-03_jprb,  6.0142409E-03_jprb, -4.4300079E-01_jprb, -2.3605519E+00_jprb,  9.1767061E-01_jprb,  1.3155031E-01_jprb, &
  6.3908338E-02_jprb, -6.6384375E-02_jprb, -3.3526097E-02_jprb,  2.5365453E-01_jprb,  4.1522743E-01_jprb,  8.3166180E-02_jprb, &
 -1.1459832E-01_jprb,  4.4557944E-01_jprb,  1.1702505E-03_jprb,  4.8450681E-03_jprb,  2.4562842E-02_jprb, -2.7920209E+00_jprb, &
  2.8232487E-02_jprb, -2.8778147E-01_jprb, -1.3564955E+00_jprb,  1.5725738E-02_jprb,  6.8140455E-03_jprb,  4.9603715E-01_jprb, &
 -6.8216688E-01_jprb,  4.7907167E-01_jprb, -1.7537064E+00_jprb,  5.6709781E+00_jprb, -6.0188577E-02_jprb, -3.1196981E-02_jprb, &
  2.4282497E+00_jprb,  9.1421971E-01_jprb,  3.2472030E+00_jprb,  5.7756426E-01_jprb, -2.3953666E-01_jprb,  3.0806863E-02_jprb, &
  1.3117870E-02_jprb,  2.5567034E+00_jprb,  1.4037983E+00_jprb,  1.4117163E-01_jprb,  3.3931331E-01_jprb,  3.7644001E+00_jprb, &
 -1.6417755E-02_jprb, -1.6291982E-02_jprb, -1.0367724E-01_jprb,  9.5359116E+00_jprb,  5.9436365E+00_jprb,  3.2583613E-01_jprb, &
 -6.1835991E-01_jprb,  6.3803117E-02_jprb,  1.4314761E-02_jprb,  2.2847757E-02_jprb,  3.4108264E-01_jprb,  1.5992651E-01_jprb, &
 -1.7994584E-01_jprb,  4.3447846E-01_jprb, -6.3250609E-02_jprb,  7.9793826E-03_jprb, -1.2341413E+00_jprb, -1.4876583E+00_jprb, &
 -2.9363224E-01_jprb, -2.5880443E-02_jprb,  2.7357307E-01_jprb,  2.0298399E-02_jprb,  1.7152941E-03_jprb,  1.1520119E+00_jprb, &
  1.6544107E+00_jprb, -6.0397770E+00_jprb,  1.0279589E+00_jprb, -6.8207974E-03_jprb,  6.0692461E-02_jprb, -1.3665172E-02_jprb, &
 -2.4633027E-01_jprb, -9.9479295E-01_jprb,  2.1739924E+00_jprb,  3.7534486E-01_jprb, -7.6896243E-01_jprb,  3.9219001E-02_jprb, &
  1.5120752E-02_jprb,  5.2904174E-02_jprb,  2.4239223E-01_jprb, -2.0949715E-01_jprb, -1.4851704E-01_jprb,  1.8768047E-01_jprb, &
 -1.3760515E-01_jprb, -4.3827466E-02_jprb,  8.9940635E-02_jprb, -1.1470065E+00_jprb,  2.4426992E-01_jprb, -5.1378440E-01_jprb, &
  4.0821245E-01_jprb,  7.2050442E-03_jprb,  9.1054881E-03_jprb, -2.6628527E-01_jprb, -1.4994969E+00_jprb, -1.9804572E-03_jprb, &
  5.2425243E-01_jprb, -1.5743700E+00_jprb, -2.6811979E-03_jprb, -3.6739941E-04_jprb, -1.8461967E-01_jprb,  5.1528520E-01_jprb, &
 -2.0373919E-01_jprb, -5.0731456E-01_jprb, -5.6206359E-01_jprb,  3.8039428E-03_jprb,  1.2360624E-02_jprb, -3.8030987E-02_jprb, &
 -3.7917469E-01_jprb,  4.4831709E-02_jprb, -2.9485185E-02_jprb, -4.3494758E-01_jprb, -3.1618972E-02_jprb, -2.4238884E-03_jprb, &
 -4.3010804E-01_jprb,  6.9159413E-01_jprb,  1.1642752E+00_jprb,  5.9654968E-02_jprb, -3.1574179E-01_jprb, -1.4829999E-01_jprb, &
  2.0069061E-01_jprb,  3.6196377E-01_jprb,  9.7201758E-01_jprb,  3.5914657E+00_jprb,  2.1788338E-01_jprb, -6.6906126E-01_jprb, &
  2.1157231E-02_jprb,  7.2160777E-03_jprb, -3.4903701E-02_jprb,  4.0905916E-02_jprb,  8.7875708E-01_jprb,  1.0434406E-01_jprb, &
 -1.5991373E-01_jprb,  2.0636785E-02_jprb,  1.9244743E-02_jprb,  2.2149412E+00_jprb,  3.1417915E+00_jprb, -4.1691029E-01_jprb, &
 -1.1992734E-01_jprb, -1.8207633E-01_jprb,  1.0199683E-01_jprb, -6.6250469E-02_jprb,  7.0477219E-01_jprb,  9.2985345E-02_jprb, &
  5.5769254E-02_jprb, -2.4142753E-01_jprb,  5.8342596E-01_jprb, -1.7699390E-03_jprb,  7.7548130E-04_jprb,  4.5042599E-02_jprb, &
 -2.1235827E+00_jprb,  2.1427991E-02_jprb, -2.9280294E+00_jprb,  7.8994428E+00_jprb,  4.8748413E-04_jprb, -3.6254758E-03_jprb, &
  3.1067939E-01_jprb, -4.5639382E-01_jprb,  3.0453009E-01_jprb,  7.2568442E-01_jprb,  2.2837062E-01_jprb, -2.9408074E-02_jprb, &
 -1.3506292E-02_jprb, -3.0424528E+00_jprb, -3.6166996E-01_jprb,  1.4150755E-01_jprb,  2.2402409E-01_jprb,  1.2852193E+00_jprb, &
 -3.6619906E-02_jprb, -1.5522820E-03_jprb, -3.7891196E-02_jprb,  5.4947069E-01_jprb,  1.4912675E+00_jprb,  1.2013916E-01_jprb, &
 -3.5484044E-01_jprb, -1.3316142E-01_jprb,  1.1117767E-01_jprb,  2.7167487E-01_jprb,  9.4568848E-01_jprb, -4.4805077E+00_jprb, &
 -3.0695915E-01_jprb,  2.3975999E-01_jprb,  2.6152724E-02_jprb,  1.2709816E-02_jprb,  1.7936441E-01_jprb,  1.4198214E-01_jprb, &
  1.4812967E-01_jprb, -8.0719872E-01_jprb,  3.0948956E+00_jprb, -4.1726944E-03_jprb, -7.8223355E-04_jprb,  1.1712041E-01_jprb, &
  4.6070650E-01_jprb,  6.2452553E-01_jprb, -5.6707535E-01_jprb, -2.5040274E-01_jprb,  4.0460534E-02_jprb,  1.1011537E-02_jprb, &
  7.2080943E-01_jprb,  2.8851281E+00_jprb,  9.6880520E-01_jprb,  6.0619615E-01_jprb,  7.6350070E-01_jprb,  5.3240039E-02_jprb, &
  4.0898133E-03_jprb,  1.3573856E-01_jprb,  1.4456979E+00_jprb, -1.1414679E-01_jprb, -5.8073787E-01_jprb, -1.0894564E-02_jprb, &
 -1.8937011E-02_jprb, -3.6301772E-03_jprb, -1.1665106E+00_jprb,  3.1690528E-01_jprb,  4.0554218E-02_jprb, -4.2358954E-02_jprb, &
  3.0129605E-01_jprb,  3.3255194E-03_jprb,  3.5411400E-03_jprb, -1.7518120E-02_jprb, -5.4846647E+00_jprb, -1.6625372E+01_jprb, &
 -1.4622501E-01_jprb,  3.9470457E-01_jprb, -7.6961313E-02_jprb, -3.5586303E-02_jprb, -3.7767436E-02_jprb, -2.0841255E-01_jprb, &
 -8.8348956E-01_jprb, -1.8747583E-03_jprb, -4.8757635E-02_jprb, -1.3748479E-02_jprb,  1.7176851E-02_jprb, -3.2409206E-01_jprb, &
 -3.0875219E-01_jprb,  9.7080876E-02_jprb, -1.4526595E-01_jprb, -2.4930063E-01_jprb,  1.2444368E-01_jprb,  3.8978087E-02_jprb, &
  4.7829022E-03_jprb,  8.8201592E-01_jprb,  2.7971436E-02_jprb, -3.5846265E-01_jprb,  1.2408907E+00_jprb,  1.2454455E-03_jprb, &
  2.6591852E-04_jprb,  3.7646532E-01_jprb, -6.6582303E-01_jprb, -3.2118420E+00_jprb, -8.2005529E-01_jprb, -1.3124842E+00_jprb, &
 -3.0345845E-01_jprb,  2.3540895E-03_jprb, -1.6979239E+00_jprb, -4.6942343E+00_jprb,  2.9449105E-02_jprb, -2.5847671E-01_jprb, &
  6.4465581E-01_jprb, -5.9474755E-03_jprb,  1.6035489E-02_jprb,  2.7972876E-01_jprb,  2.1739975E-01_jprb, -9.0417050E+00_jprb, &
  5.4598167E-02_jprb,  3.0661607E-01_jprb, -1.0620120E-01_jprb, -3.6872650E-02_jprb,  8.3601561E-02_jprb, -5.1719015E-01_jprb, &
  2.3589011E+00_jprb,  5.8289473E-02_jprb, -5.2028314E-01_jprb,  1.1907135E-01_jprb,  4.0710086E-02_jprb, -4.5598321E-02_jprb, &
  5.5453378E-01_jprb, -4.6334351E+00_jprb, -2.7236286E-01_jprb, -5.9172634E-01_jprb, -1.6055487E-04_jprb, -1.8411327E-02_jprb, &
 -9.0282156E-02_jprb, -1.3587436E-01_jprb,  4.9285580E-02_jprb, -5.7404223E+00_jprb,  8.9380075E+00_jprb, -3.9857950E-03_jprb, &
  1.5032942E-03_jprb,  1.5135236E-01_jprb, -8.0159295E-01_jprb,  1.5375947E+00_jprb,  4.0553592E-02_jprb, -2.6658804E-01_jprb, &
  9.9799508E-03_jprb,  1.4015807E-02_jprb,  2.2733468E-01_jprb, -5.9703967E-01_jprb, -6.0503517E-02_jprb,  6.0654278E-01_jprb, &
 -2.0841470E+00_jprb,  1.6256251E-03_jprb,  9.2718282E-04_jprb, -5.2378909E-02_jprb, -8.1125525E-02_jprb,  2.1633171E-01_jprb, &
 -1.8654715E+00_jprb,  5.7318989E+00_jprb, -2.6615986E-02_jprb,  8.8722901E-05_jprb,  1.2835792E+00_jprb, -6.5161172E-01_jprb, &
  6.5372845E-01_jprb, -1.9141577E-01_jprb, -6.6260466E-02_jprb,  3.4861201E-02_jprb,  1.0987446E-02_jprb,  2.1475183E+00_jprb, &
  9.3697753E-01_jprb,  3.0399957E+00_jprb,  3.8236416E-01_jprb, -3.5797355E-01_jprb, -1.0780960E-02_jprb, -1.1692279E-02_jprb, &
 -2.2847213E-01_jprb, -9.1613301E-02_jprb,  9.4293444E+00_jprb,  9.8554563E-02_jprb,  4.0281830E-01_jprb, -4.0576527E-02_jprb, &
  4.1247762E-02_jprb,  2.5342129E-01_jprb,  5.1499208E-01_jprb,  5.5255927E-01_jprb, -5.6296342E-01_jprb,  5.4032307E-01_jprb, &
  6.2675703E-02_jprb,  2.7711623E-02_jprb, -1.3426634E+00_jprb, -8.7670313E-01_jprb, -2.3564122E+00_jprb, -1.1447970E-01_jprb, &
 -7.3469855E-01_jprb,  2.1526674E-01_jprb,  7.2052594E-02_jprb, -8.9075952E-02_jprb, -3.2813806E-01_jprb, -3.5362080E-01_jprb, &
 -1.2954758E-01_jprb, -1.6240606E+00_jprb,  2.1148788E-02_jprb,  1.5808139E-03_jprb, -2.6599960E-01_jprb, -8.0578427E-01_jprb, &
  3.6194864E-02_jprb, -5.3064466E-01_jprb,  8.5358553E-01_jprb, -6.7360252E-03_jprb,  7.4860475E-04_jprb, -1.1340430E-02_jprb, &
 -9.6670016E-01_jprb, -6.8903011E+00_jprb,  1.4939259E+00_jprb, -1.1983082E-01_jprb,  4.2287723E-02_jprb, -3.4447620E-02_jprb, &
 -3.6171301E-01_jprb, -1.0422722E+00_jprb, -9.0393029E-02_jprb,  7.0552982E+00_jprb, -1.8221334E-02_jprb, -2.0932892E-02_jprb, &
 -3.3229899E-03_jprb, -1.4958469E+00_jprb,  1.3030012E+00_jprb,  3.3959820E-02_jprb, -1.6699984E-01_jprb,  2.7472273E+00_jprb, &
  9.5285755E-03_jprb,  1.8165741E-04_jprb,  1.3492964E-01_jprb, -2.1871687E-01_jprb,  1.0358026E+01_jprb, -7.2157356E-03_jprb, &
 -4.2758227E-01_jprb,  8.6869658E-02_jprb,  2.7935838E-02_jprb, -4.9186890E-02_jprb,  3.3489223E-01_jprb,  9.4630835E-02_jprb, &
 -6.3484158E-01_jprb, -1.8724044E-01_jprb, -6.6914748E-03_jprb, -1.2532880E-02_jprb,  2.4215824E-01_jprb,  7.9953086E+00_jprb, &
 -7.6932166E+00_jprb, -2.4350150E-01_jprb,  7.2747125E-01_jprb, -8.1991099E-02_jprb, -2.3057340E-02_jprb,  7.6138640E-02_jprb, &
 -4.1503507E-01_jprb, -1.5398207E+00_jprb, -7.4316277E-02_jprb,  1.3960971E-01_jprb, -1.3919671E-01_jprb, -2.8994547E-02_jprb, &
 -6.6944952E-01_jprb, -1.3421650E+00_jprb, -1.5127416E+00_jprb, -2.2669918E-01_jprb,  1.0241273E+00_jprb,  4.6600672E-02_jprb, &
 -9.6560461E-02_jprb, -2.6314054E+00_jprb, -5.2977737E+00_jprb, -1.0295157E-01_jprb, -3.3260249E-01_jprb, -5.7031521E-01_jprb, &
  4.0467973E-04_jprb,  2.0620752E-03_jprb, -6.2705115E-02_jprb,  3.0474430E-02_jprb, -2.6909807E+00_jprb, -8.0701894E-02_jprb, &
 -9.8291550E-01_jprb,  5.1024521E-02_jprb, -8.1175810E-03_jprb, -8.9602048E-02_jprb, -1.6869999E-01_jprb,  2.5531226E-02_jprb, &
 -8.6561149E-01_jprb,  2.6824531E+00_jprb,  1.5118301E-02_jprb, -3.7793072E-03_jprb,  7.9248630E-01_jprb, -1.0272194E+00_jprb, &
  3.0132789E-03_jprb,  4.8660360E-01_jprb, -1.7552170E+00_jprb, -7.9023638E-03_jprb,  1.9259531E-04_jprb, -5.1361150E-02_jprb, &
  2.7240168E-01_jprb,  1.9307015E-01_jprb,  3.2952760E-01_jprb,  4.6465587E+00_jprb, -1.8084664E-02_jprb, -2.3087786E-02_jprb, &
 -1.2398482E-01_jprb,  1.1627164E+01_jprb, -9.3130358E-01_jprb,  1.2402023E-01_jprb,  2.1843419E-01_jprb,  1.4046520E-01_jprb, &
 -2.2412138E-01_jprb, -4.1547850E-01_jprb, -8.0128690E-01_jprb, -7.7537580E-01_jprb,  3.3705935E-03_jprb,  1.1172207E-01_jprb, &
 -3.6451607E-02_jprb, -1.5155883E-02_jprb, -2.1974677E+00_jprb, -2.0195108E+00_jprb,  1.7480473E+01_jprb,  1.5198164E-01_jprb, &
 -3.3275369E-01_jprb,  6.0960007E-02_jprb,  3.8435609E-02_jprb,  6.5060378E-02_jprb,  1.6908065E-01_jprb,  3.1404407E-01_jprb, &
 -1.7195407E+00_jprb,  5.9401981E+00_jprb,  4.8960134E-04_jprb, -2.5282782E-03_jprb,  1.5704596E-01_jprb,  8.2146673E-01_jprb, &
 -2.6963220E-01_jprb,  1.0316299E-01_jprb, -7.8633342E-03_jprb,  3.9665297E-02_jprb,  9.0251337E-04_jprb,  1.3614607E+00_jprb, &
  5.5240904E-02_jprb,  2.6122879E-02_jprb,  6.4941070E-01_jprb,  3.9142009E-01_jprb, -5.4612892E-03_jprb,  5.0998698E-03_jprb, &
  4.7325251E-01_jprb,  1.9401159E-01_jprb, -2.6262442E+00_jprb, -9.9456206E-02_jprb, -9.1531572E-01_jprb,  6.5265764E-02_jprb, &
  8.8829989E-03_jprb, -6.1107446E-02_jprb, -2.8902380E-01_jprb, -8.9422231E-01_jprb, -3.2824570E-01_jprb,  1.8171322E-01_jprb, &
  2.3476541E-02_jprb, -2.5545885E-02_jprb,  1.2869660E-02_jprb, -6.7547089E-01_jprb, -1.0651669E-01_jprb, -3.4639355E-01_jprb, &
 -1.9249443E+00_jprb,  1.3923016E-02_jprb,  3.7440425E-03_jprb,  6.2531026E-02_jprb, -3.4022750E+00_jprb, -7.5583472E-02_jprb, &
  3.5735262E+00_jprb, -9.2920972E+00_jprb,  9.6633285E-04_jprb,  2.8692954E-03_jprb, -5.8752509E-01_jprb,  8.6848566E-01_jprb /

! bias vector of first layer, vector [no * 1]
! 2 elements in b2
DATA iso_net%b2 /        -3.9442164e+00_jprb,  -2.1429341e+00_jprb /

! weight matrix of first layer, matrix [nh * no]
! (180 * 2) 360 elements in W2
DATA iso_net%W2 /     4.0051335E-01_jprb, -1.8668008E+00_jprb, -2.4449709E+00_jprb,  3.5524317E+00_jprb,  1.7079847E+00_jprb, &
-4.3338332E-01_jprb, -2.0120590E+00_jprb, -1.6663834E+00_jprb,  2.9625577E+00_jprb, -6.1856568E+00_jprb, -5.2609299E+00_jprb, &
-2.3839230E+00_jprb,  8.2927298E-01_jprb,  1.4376675E+00_jprb, -9.9442478E-01_jprb,  1.3558373E+00_jprb, -9.7257268E-01_jprb, &
 9.1198314E-01_jprb,  5.3063500E+00_jprb,  1.6709361E+00_jprb, -9.4159497E-02_jprb,  5.6358221E-01_jprb, -7.4365759E-02_jprb, &
 1.0774980E+00_jprb,  3.4144870E+00_jprb, -4.6889404E-01_jprb, -2.6256608E+00_jprb, -1.3102189E+00_jprb,  1.3227473E+00_jprb, &
-1.1739720E-01_jprb, -2.5768090E+00_jprb,  2.8534561E+00_jprb, -7.1597615E+00_jprb, -6.0735648E-01_jprb,  2.0123086E+00_jprb, &
-1.3939754E+00_jprb, -9.8210578E-01_jprb,  2.3132509E+00_jprb,  3.4095134E+00_jprb, -5.7646066E-01_jprb,  1.6412997E+00_jprb, &
 3.9648273E+00_jprb,  1.5763106E+00_jprb,  8.6477652E-01_jprb, -6.7708789E-01_jprb,  3.0694872E+00_jprb, -1.2862605E+00_jprb, &
-2.4711500E+00_jprb, -2.3673916E-02_jprb,  2.3111727E+00_jprb, -1.1666559E+00_jprb, -3.4948257E+00_jprb, -7.8943401E-01_jprb, &
 3.8086712E+00_jprb,  1.8992240E+00_jprb, -8.7933218E-01_jprb,  3.2413581E+00_jprb,  3.5041918E+00_jprb,  3.0070329E+00_jprb, &
-2.3920377E+00_jprb,  2.9350273E+00_jprb, -3.4652745E-02_jprb,  2.8893640E+00_jprb, -1.4327895E+00_jprb, -1.8702297E+00_jprb, &
-2.9362270E+00_jprb, -2.5564287E-01_jprb, -2.9213633E+00_jprb, -5.5681465E-01_jprb,  2.1842184E+00_jprb, -2.1877956E+00_jprb, &
-2.9563491E-02_jprb, -1.3953482E+00_jprb,  1.4025152E+00_jprb,  1.2642776E+00_jprb,  3.2676204E+00_jprb, -1.7221933E+00_jprb, &
 2.4000951E-01_jprb, -1.2475734E+00_jprb,  2.6407788E+00_jprb, -1.3801703E+00_jprb, -6.5628459E-01_jprb,  2.0977593E+00_jprb, &
-1.6668152E-01_jprb,  5.8359843E-01_jprb, -1.3127595E+00_jprb,  7.4200087E-02_jprb,  2.1051577E+00_jprb, -7.2426263E-01_jprb, &
-1.0959044E+00_jprb, -4.7162552E-01_jprb,  3.4118371E-01_jprb, -1.9949195E+00_jprb, -1.4555299E-01_jprb, -8.2389801E-01_jprb, &
 2.1690292E+00_jprb,  2.7450814E+00_jprb, -8.1057946E-01_jprb,  3.8422445E+00_jprb, -2.2732156E+00_jprb,  5.4573421E-01_jprb, &
 3.2747448E+00_jprb,  3.1128990E+00_jprb, -1.4632424E+00_jprb, -9.9015895E-01_jprb, -3.6746345E+00_jprb,  7.5257562E-01_jprb, &
 3.3992296E+00_jprb, -1.5430161E+00_jprb, -9.3881065E-02_jprb,  1.0501996E-01_jprb,  8.7508654E+00_jprb, -1.0119535E+00_jprb, &
-3.7868059E-01_jprb,  7.2557779E-01_jprb, -5.2832410E+00_jprb,  1.5926276E+00_jprb, -1.1770401E+00_jprb, -2.5987006E+00_jprb, &
-8.0324571E-01_jprb, -6.8335252E-01_jprb, -1.4425510E+00_jprb, -6.3334696E-01_jprb,  2.5225613E+00_jprb, -7.0825419E-01_jprb, &
-4.0343268E-01_jprb,  3.0935324E+00_jprb, -6.5521722E+00_jprb, -1.3148801E+00_jprb, -1.1637123E+00_jprb,  5.8300894E-01_jprb, &
-2.8856674E+00_jprb, -9.3452184E-01_jprb, -7.5804919E-01_jprb,  6.5970632E-01_jprb, -1.0128075E+00_jprb,  4.2750246E+00_jprb, &
 4.2776113E+00_jprb, -8.0805663E-01_jprb, -1.2334613E+00_jprb,  3.0375364E+00_jprb, -1.7695951E-02_jprb,  6.1782596E-01_jprb, &
-5.0020094E+00_jprb,  1.1966881E+00_jprb,  7.8535005E-01_jprb,  1.1862574E+01_jprb,  1.9620570E+00_jprb, -5.5064987E-02_jprb, &
-2.9161469E+00_jprb, -3.6154611E+00_jprb,  2.8174288E+00_jprb,  7.0385122E-01_jprb,  2.2790921E-01_jprb, -3.0873080E-01_jprb, &
 3.2653492E+00_jprb, -3.4617421E+00_jprb,  3.0781994E+00_jprb,  3.9455831E+00_jprb,  8.5267604E-02_jprb, -4.2221678E+00_jprb, &
-2.8673239E+00_jprb, -4.4405557E-02_jprb,  1.5407482E+00_jprb, -9.2303189E-01_jprb,  2.3324287E+00_jprb, -2.0538062E+00_jprb, &
-4.6282526E+00_jprb,  8.4558750E-01_jprb, -3.4403222E+00_jprb, -3.9702581E-01_jprb, -3.0031833E+00_jprb,  4.1715506E+00_jprb, &
-1.6511178E-01_jprb,  1.3453369E+00_jprb, -4.2823682E+00_jprb,  1.6054503E+00_jprb, -2.3218583E-01_jprb,  9.4098870E+00_jprb, &
-8.9086426E+00_jprb,  3.4659814E+00_jprb, -1.6868253E+00_jprb, -1.4857916E+00_jprb, -1.2820893E+00_jprb,  9.7424298E-01_jprb, &
 5.9493796E-02_jprb, -4.7365029E-01_jprb, -1.1064680E+00_jprb,  1.1331303E+00_jprb,  7.9448179E-01_jprb, -2.0262909E+00_jprb, &
-1.3684994E+00_jprb, -3.5802379E+00_jprb,  7.9491502E-01_jprb, -6.2906869E-01_jprb,  2.9467405E+00_jprb,  2.7230751E+00_jprb, &
-4.8486171E-01_jprb,  1.8955042E+00_jprb,  1.7635870E+00_jprb, -1.9145816E-01_jprb, -2.8686357E+00_jprb,  1.5852710E+00_jprb, &
 1.1192260E-01_jprb,  2.5473577E+00_jprb,  1.4475370E+00_jprb, -9.2995167E-02_jprb,  3.2357207E-01_jprb,  1.2713028E+00_jprb, &
-1.6474820E-01_jprb, -1.7731700E+00_jprb, -4.6824645E-01_jprb, -2.2250242E+00_jprb,  9.5490586E-02_jprb, -6.0078334E-01_jprb, &
-7.8766890E-01_jprb, -4.3287525E+00_jprb, -1.1179341E+00_jprb,  9.4209561E-01_jprb, -8.3299613E-01_jprb,  6.9371591E-01_jprb, &
 3.2034936E+00_jprb,  1.3417669E+00_jprb,  9.2474657E-01_jprb,  5.7989613E-01_jprb,  1.1460032E+00_jprb,  5.7563722E-02_jprb, &
-8.8332232E-01_jprb, -2.3822572E-02_jprb,  2.2869219E+00_jprb, -2.6095592E-01_jprb, -1.4361747E+00_jprb, -1.7337292E+00_jprb, &
 7.7196183E-01_jprb, -3.8757090E-01_jprb, -3.9464259E-01_jprb,  7.2202456E-01_jprb,  3.4954400E-01_jprb,  2.1592154E+00_jprb, &
 3.9559730E+00_jprb,  1.2861175E+00_jprb, -3.1170220E-02_jprb, -3.4669003E+00_jprb, -4.8484359E-01_jprb, -9.8833299E-01_jprb, &
 3.4811806E+00_jprb, -1.9011638E+00_jprb, -7.6287357E-01_jprb,  3.0237158E-01_jprb, -1.9901722E-01_jprb, -3.6451027E-01_jprb, &
-4.9558536E-02_jprb,  2.9152486E-02_jprb,  1.0492160E+00_jprb,  1.3630668E+00_jprb, -1.4003403E-01_jprb, -1.3485711E-01_jprb, &
 7.7979248E-02_jprb, -8.4341574E-01_jprb,  7.6531774E-01_jprb, -1.3123982E+00_jprb, -4.0682382E-01_jprb,  4.4636549E-01_jprb, &
 7.0117189E-02_jprb,  2.8055727E-01_jprb, -2.2104825E-01_jprb, -1.1499469E-01_jprb, -1.1989322E-01_jprb, -2.9172967E-01_jprb, &
-2.1102025E+00_jprb, -1.1400687E+00_jprb, -1.5901726E+00_jprb, -9.7429858E-01_jprb,  1.4323956E+00_jprb, -1.0017952E+00_jprb, &
 5.4167462E-01_jprb, -7.4976494E-01_jprb, -2.9233540E-01_jprb,  9.7553178E-01_jprb, -8.4063147E-01_jprb,  1.4426585E-01_jprb, &
 1.6821396E+00_jprb,  1.1770997E+00_jprb, -8.4413142E-01_jprb,  1.9697307E+00_jprb, -1.6966270E+00_jprb,  8.4304108E-01_jprb, &
 1.4071191E+00_jprb,  4.2991387E-01_jprb, -1.8649945E-03_jprb,  4.7934584E-02_jprb,  3.1911443E+00_jprb, -7.5809631E-01_jprb, &
-1.4400161E-01_jprb,  2.0553449E-01_jprb, -1.6014492E+00_jprb, -3.4223521E-01_jprb,  3.5391724E-01_jprb, -7.1795515E-01_jprb, &
 3.2409437E+00_jprb, -1.0619946E+00_jprb,  9.8927685E-01_jprb, -9.0046437E-01_jprb,  1.3619804E+00_jprb, -2.0812317E+00_jprb, &
-5.4215728E-01_jprb,  1.3494155E+00_jprb, -1.7802096E-01_jprb, -6.0060793E-01_jprb, -7.7160925E-01_jprb,  7.3288144E-01_jprb, &
-3.2877224E+00_jprb,  2.5722128E+00_jprb, -3.7674298E-01_jprb,  5.0301896E-01_jprb, -2.5785280E-01_jprb,  1.3768267E+00_jprb, &
 7.0291312E+00_jprb,  7.2057180E-01_jprb, -1.8257160E+00_jprb,  3.2756720E+00_jprb, -1.7520966E-02_jprb,  7.8502008E-01_jprb, &
-1.0111239E+00_jprb,  2.1915733E-01_jprb,  3.9864930E+00_jprb,  1.3903068E-02_jprb,  1.2842439E+00_jprb, -2.8077977E+00_jprb, &
 1.2737483E-01_jprb, -1.5589265E+00_jprb, -5.3101889E-01_jprb,  5.9568399E-01_jprb,  1.0714928E-01_jprb, -3.8397623E-01_jprb, &
 1.7852450E+00_jprb, -8.8076044E-01_jprb,  1.2247539E+00_jprb,  9.7905802E-01_jprb,  1.2942658E-01_jprb,  2.4822568E+00_jprb, &
-9.8259966E-01_jprb, -1.7379550E+00_jprb,  1.0923052E+00_jprb, -9.7212427E-01_jprb,  3.3492789E+00_jprb,  8.5103540E-02_jprb, &
 2.0516115E+00_jprb, -1.3152380E+00_jprb, -1.2395504E+00_jprb, -6.3840442E-01_jprb, -1.7059156E+00_jprb,  6.0866950E+00_jprb, &
 1.1415179E+00_jprb,  2.6117314E-01_jprb,  1.2156630E-02_jprb,  1.8643292E+00_jprb,  1.6979459E+00_jprb,  3.2320472E+00_jprb, &
 1.0162774E+00_jprb /     

!== ti: structure having the settings of the input mapminmax transformation

! 7 elements in ti_gain
DATA iso_net%ti_gain / 2.8591851e-03_jprb, 2.2471910e-02_jprb, 4.0000000e-02_jprb, 6.2500000e-02_jprb, &
                       5.0632911e-02_jprb, 2.3269984e+00_jprb, 2.2429120e+00_jprb /

! 7 elements in ti_xoffset 
DATA iso_net%ti_xoffset / 5.0000000e-01_jprb, 0.0000000e+00_jprb, 0.0000000e+00_jprb, 2.7115000e+02_jprb, &
                          0.0000000e+00_jprb, 1.4029873e-01_jprb, 3.5592941e-03_jprb /

! 1 element in ti_ymin
DATA iso_net%ti_ymin / -1.0000000e+00_jprb /

!== to: structure having the settings of the output mapminmax transformation

! 2 elements in to_gain
DATA iso_net%to_gain / 1.9594736e+00_jprb, 3.2770936e+00_jprb /

! 2 elements in to_xoffset 
DATA iso_net%to_xoffset / -2.8205421e-01_jprb, -1.1063987e-05_jprb /

! 1 element in to_ymin
DATA iso_net%to_ymin / -1.0000000e+00_jprb /

!!===== anisotropic coefficients ====================================

!> Neural network 
TYPE(surfem_aniso_net), PUBLIC, TARGET :: aniso_net

! bias vector of first layer, vector [nh * 1]
! 180 elements in b1
DATA aniso_net%b1 / &
-1.9641305E+00_jprb,  6.0634385E+00_jprb, -2.6917618E+00_jprb,  3.4555850E+00_jprb, -5.9022953E+00_jprb,  9.8907194E+00_jprb, &
 2.5889827E+00_jprb,  1.9383002E+00_jprb,  1.7317644E+00_jprb,  1.1334977E+01_jprb, -7.3444923E+00_jprb,  9.4736195E-01_jprb, &
 7.0209243E+00_jprb,  1.2017111E+00_jprb, -4.3056451E+00_jprb, -2.1026714E+00_jprb, -8.1118486E-01_jprb, -5.3078677E-01_jprb, &
-2.3788259E+00_jprb, -6.6812263E+00_jprb, -3.2967895E+00_jprb,  6.7667306E-01_jprb, -1.6325385E+00_jprb, -6.9872842E+00_jprb, &
 2.3400057E+00_jprb,  1.9331680E+00_jprb, -3.8812285E+00_jprb, -1.0720582E+01_jprb,  6.6562992E+00_jprb,  1.0876691E+00_jprb, &
 7.8879763E+00_jprb, -4.9431161E+00_jprb, -1.1948736E+00_jprb, -1.1166555E-01_jprb,  9.9005466E-01_jprb,  1.7219734E+00_jprb, &
 1.7658274E-01_jprb, -1.0019979E+00_jprb, -1.0358301E+00_jprb,  1.1046754E+00_jprb, -2.3190570E+00_jprb, -1.0182093E+00_jprb, &
-2.0857566E+00_jprb,  2.5291120E+00_jprb,  3.4281048E-01_jprb,  1.0551091E+00_jprb,  6.3554417E-02_jprb,  2.7562264E+00_jprb, &
 6.7626099E+00_jprb,  1.5630068E-01_jprb, -3.1877531E-04_jprb,  4.2802999E-02_jprb,  1.5832587E+00_jprb,  4.2463906E+00_jprb, &
 3.2002853E+00_jprb, -9.3892945E-01_jprb,  8.1996227E-01_jprb, -2.0253547E-01_jprb, -7.0112952E-01_jprb, -1.9847849E-02_jprb, &
 7.9983332E-02_jprb, -1.0260324E+00_jprb, -1.3999837E+00_jprb, -2.8719483E+00_jprb,  9.1926764E-01_jprb,  6.9361864E-01_jprb, &
 3.0518780E+00_jprb,  3.1249964E+00_jprb,  3.2543236E+00_jprb,  2.1205235E-01_jprb, -2.7042040E-01_jprb, -1.3411568E+00_jprb, &
-3.2227829E+00_jprb, -2.5768499E-01_jprb, -6.3554839E-01_jprb,  8.0911607E+00_jprb,  3.0096621E+00_jprb, -6.6176106E-01_jprb, &
 5.1063953E-01_jprb,  1.9520113E+00_jprb,  3.2829173E+00_jprb,  9.9638923E-01_jprb,  1.0915880E+00_jprb,  5.8188203E-01_jprb, &
 6.2728865E+00_jprb,  1.0167847E+00_jprb, -4.9122767E-01_jprb, -8.8297594E+00_jprb,  1.3037281E+00_jprb,  2.9520425E+00_jprb, &
 7.3135492E-01_jprb, -5.6833504E+00_jprb, -2.9146627E-02_jprb,  7.6787372E-01_jprb,  1.6872276E+00_jprb, -1.7008671E+00_jprb, &
-9.7060365E-01_jprb, -7.1107119E+00_jprb,  3.8287934E-01_jprb, -2.9498688E+00_jprb,  5.5186299E+00_jprb,  2.2758984E+00_jprb, &
 2.1669435E-02_jprb,  1.1412920E+00_jprb,  1.3712648E+00_jprb, -3.8569292E+00_jprb,  9.3677190E-01_jprb,  3.1943288E+00_jprb, &
 2.1029758E+00_jprb, -5.4429616E+00_jprb,  1.1624348E+00_jprb,  6.6219675E+00_jprb,  6.4081209E+00_jprb, -1.5770647E+00_jprb, &
 7.3249017E-01_jprb, -3.0644386E+00_jprb,  2.7139675E+00_jprb,  2.9930915E+00_jprb,  1.1010092E+00_jprb,  9.6701796E-02_jprb, &
 8.6045610E-01_jprb, -9.3959018E-01_jprb,  6.5706360E-01_jprb,  5.3323127E+00_jprb,  7.1435944E+00_jprb, -3.2518370E-01_jprb, &
-1.1124667E+01_jprb,  3.2691174E+00_jprb, -3.1589439E+00_jprb, -1.2241433E+00_jprb,  2.3287587E+00_jprb, -1.0353912E+01_jprb, &
 3.3138516E+00_jprb,  1.0695740E+00_jprb,  8.6432872E-01_jprb, -1.0057410E+00_jprb, -4.7141546E+00_jprb, -6.7874471E+00_jprb, &
-1.5477423E+00_jprb,  1.4823472E+01_jprb, -9.7092495E-01_jprb, -3.5744834E+00_jprb, -6.7181225E+00_jprb, -6.5744306E+00_jprb, &
 8.7513815E+00_jprb, -8.1874333E-01_jprb,  3.1413491E+00_jprb,  3.7530157E+00_jprb,  9.7222610E+00_jprb,  2.4274300E+00_jprb, &
 9.6353340E+00_jprb, -4.5699711E-01_jprb,  1.2294318E+01_jprb, -3.1519799E+00_jprb, -3.2713250E+00_jprb, -1.2403267E+00_jprb, &
-7.8047774E+00_jprb, -2.6576320E+00_jprb, -4.6254434E+00_jprb, -9.4023298E-01_jprb,  2.2251087E+01_jprb,  3.9958426E+00_jprb, &
-5.8683599E+00_jprb, -2.2094257E+00_jprb,  3.0892848E+00_jprb, -5.3763005E+00_jprb, -6.2702709E+00_jprb,  2.8773983E+00_jprb, &
 1.9486821E+01_jprb,  3.8565125E+00_jprb, -3.9363093E+00_jprb, -7.7498658E+00_jprb,  4.9322726E+00_jprb,  1.2259897E+01_jprb, &
-8.5269383E-01_jprb,  6.4433446E+00_jprb, -4.0478735E+00_jprb, -7.2798289E+00_jprb, -1.3575911E+01_jprb, -1.1371588E+01_jprb /

! weight matrix of first layer, matrix [nh (columns) * ni (rows)]
! 180 * 7 elements in W1 
DATA aniso_net%W1 / &  
  1.6079875E-01_jprb,  2.2618774E+00_jprb,  5.9596311E-01_jprb, -2.0646864E-02_jprb,  4.5072675E-03_jprb, -1.8807262E-01_jprb, &
  1.4647902E+00_jprb,  3.5222099E-01_jprb, -1.5481654E+00_jprb, -5.1262468E-01_jprb,  1.0379705E-01_jprb,  1.8388250E-02_jprb, &
  1.1623907E-01_jprb,  3.0270454E+00_jprb, -3.0854046E-01_jprb,  1.8590806E+00_jprb,  7.2450309E-01_jprb, -2.3979270E-01_jprb, &
 -3.4447214E-01_jprb, -4.3861193E-01_jprb,  1.3050657E+00_jprb, -8.6080599E-02_jprb, -2.7518169E-02_jprb,  3.3337040E+00_jprb, &
 -2.8270902E-02_jprb,  5.0550023E-03_jprb, -4.5102365E-01_jprb, -2.7272360E-01_jprb, -3.6285667E+00_jprb,  5.4368324E-02_jprb, &
  2.9887411E-01_jprb, -2.6459723E-02_jprb, -2.9365148E-02_jprb, -1.3137931E+00_jprb, -1.0761417E+00_jprb, -3.9004637E-01_jprb, &
 -3.3464486E-01_jprb,  1.0186497E+01_jprb,  3.3043395E-02_jprb, -9.1932172E-03_jprb, -1.3465028E-01_jprb,  5.7650332E-01_jprb, &
  6.1770358E-01_jprb, -1.1067516E+00_jprb, -2.3276503E-01_jprb, -9.3718410E-02_jprb,  1.1961423E-03_jprb,  1.6472286E-01_jprb, &
  6.9975468E-01_jprb, -6.8843285E-02_jprb, -4.6798938E-01_jprb,  1.2724358E+00_jprb, -8.6706634E-03_jprb,  4.6379374E-03_jprb, &
 -3.5773808E-01_jprb, -2.7821031E-01_jprb, -1.0460106E-01_jprb, -1.3723440E+00_jprb,  5.3119622E-01_jprb,  1.1424593E-02_jprb, &
 -6.5447337E-04_jprb,  8.4309762E-02_jprb,  3.4771120E-01_jprb,  9.4985616E+00_jprb, -1.5918728E+00_jprb,  2.5987021E-02_jprb, &
  1.0633786E-02_jprb,  7.7179442E-03_jprb,  5.9388072E-02_jprb, -2.9462088E-02_jprb,  2.6318170E+00_jprb,  9.3913624E-01_jprb, &
 -9.7368051E+00_jprb, -9.4562893E-03_jprb, -6.5438754E-03_jprb,  3.7853291E-03_jprb, -2.5397318E-01_jprb, -2.3773439E-01_jprb, &
  5.9473028E-02_jprb,  1.3209130E+00_jprb, -1.2147994E-02_jprb,  3.0958857E-03_jprb,  8.7052986E-03_jprb, -4.1787334E-01_jprb, &
 -2.6048031E+00_jprb, -9.4515903E-01_jprb,  9.2270213E+00_jprb,  6.0887208E-03_jprb,  7.2098443E-03_jprb,  9.0114600E-03_jprb, &
  1.2021808E-01_jprb, -4.2655960E-02_jprb,  2.4152649E-01_jprb,  1.1870559E+00_jprb, -2.2485030E-02_jprb,  3.2349588E-03_jprb, &
 -5.1531301E-02_jprb, -6.9110936E-01_jprb, -2.1664572E-01_jprb,  1.6022443E+00_jprb, -2.7323285E+00_jprb,  5.2682159E-03_jprb, &
 -6.0266454E-03_jprb, -4.6827645E-02_jprb, -2.3024925E-01_jprb, -1.1810330E-01_jprb,  3.5782705E-01_jprb,  1.0592177E+00_jprb, &
  7.7124278E-02_jprb, -1.7991485E-02_jprb,  6.1113809E-01_jprb,  8.7369182E-01_jprb,  2.6238384E-01_jprb,  3.0054780E-01_jprb, &
  4.0426557E-01_jprb,  4.3603115E-03_jprb, -4.2322534E-03_jprb,  4.2826991E-01_jprb,  2.3234085E-01_jprb,  7.1461494E-01_jprb, &
  2.6614395E-01_jprb, -1.1116459E+00_jprb,  4.6483059E-02_jprb, -4.7730454E-04_jprb,  4.4369170E-01_jprb,  3.9763091E-01_jprb, &
 -1.1316805E-01_jprb,  2.5936695E+00_jprb,  7.6698661E-01_jprb,  2.1843452E-02_jprb, -1.3749144E-02_jprb, -2.8288597E-02_jprb, &
  2.0513294E-01_jprb, -4.3506379E+00_jprb, -9.4513188E-02_jprb,  1.5758635E-01_jprb, -3.4087436E-02_jprb, -9.5802198E-03_jprb, &
  9.9631587E-01_jprb,  2.7008118E-01_jprb,  7.0382651E-01_jprb,  1.9206089E+00_jprb,  8.6843047E-01_jprb,  6.9748745E-01_jprb, &
  5.9918492E-01_jprb, -2.6835494E-01_jprb,  6.9502678E-01_jprb,  1.2258529E+00_jprb, -6.5418384E-01_jprb, -2.8924702E-01_jprb, &
  2.2194884E-02_jprb, -6.2095270E-03_jprb,  3.9680645E-01_jprb,  2.5978740E-01_jprb, -3.0887248E-01_jprb, -5.5523260E-01_jprb, &
  1.2385260E+00_jprb,  2.0968990E-02_jprb, -1.9594990E-02_jprb,  5.0251653E-01_jprb,  3.7260052E-01_jprb, -3.6620292E+00_jprb, &
  1.4240720E-03_jprb,  1.3430822E-01_jprb, -1.9763243E-02_jprb, -3.8110924E-02_jprb, -1.3898374E+00_jprb, -2.0400448E+00_jprb, &
  3.2512402E-02_jprb, -2.0462659E+00_jprb, -7.9456998E-01_jprb,  2.2816905E-01_jprb,  3.4511744E-01_jprb,  5.2639752E-01_jprb, &
 -1.2781350E+00_jprb,  1.1350052E+00_jprb, -8.6509366E-01_jprb,  4.8294059E-01_jprb,  4.9675365E-02_jprb, -1.1463292E-03_jprb, &
  9.1015085E-02_jprb,  3.1555728E-01_jprb,  1.9747621E-01_jprb,  2.6957835E-01_jprb, -3.9324431E+00_jprb,  3.7593467E-02_jprb, &
 -7.2598345E-03_jprb,  3.7833765E-01_jprb,  4.4871732E-01_jprb,  4.1175797E-01_jprb,  3.4511354E-01_jprb, -1.1148213E+01_jprb, &
 -3.2012435E-02_jprb,  1.0195093E-02_jprb,  1.5230390E-01_jprb, -6.7458364E-01_jprb,  1.8397877E+00_jprb,  6.5738279E-01_jprb, &
 -9.7342938E-02_jprb, -7.7705972E-02_jprb,  9.9970274E-02_jprb,  4.3189337E+00_jprb,  1.0235943E+01_jprb,  1.9737500E+00_jprb, &
  3.9923186E-01_jprb, -5.7481941E-01_jprb, -1.4278632E-02_jprb,  3.7226083E-04_jprb,  8.4530794E-02_jprb,  4.0080352E-01_jprb, &
 -6.4644700E-02_jprb, -2.8219701E-01_jprb,  8.1299789E+00_jprb,  3.0065927E-02_jprb, -3.7814952E-03_jprb,  4.2228771E-01_jprb, &
 -2.8787311E-01_jprb, -1.5918794E+00_jprb,  4.9012946E+00_jprb,  3.3618307E-01_jprb, -3.5508458E-02_jprb, -7.6241708E-04_jprb, &
 -1.5023384E-01_jprb,  9.3306876E-01_jprb, -2.5422505E-01_jprb, -1.0721748E-01_jprb,  5.9522114E-01_jprb,  1.6862277E-03_jprb, &
 -5.2645290E-05_jprb,  3.3705034E-01_jprb, -6.5469004E-01_jprb, -1.1198040E-01_jprb,  2.8192171E-01_jprb, -1.8637637E-02_jprb, &
 -8.9886080E-04_jprb,  2.8439296E-03_jprb, -2.2671311E-01_jprb, -8.8189749E-01_jprb, -2.2908993E-01_jprb, -1.7127547E+00_jprb, &
 -5.5002089E-01_jprb, -6.4037952E-02_jprb,  8.5509623E-03_jprb, -8.9465500E-02_jprb, -6.7452282E-01_jprb, -1.0148470E-01_jprb, &
 -2.0990826E+00_jprb, -5.7252930E-01_jprb,  1.6984620E-02_jprb, -2.4166820E-03_jprb,  2.0791146E-01_jprb, -1.1712421E+00_jprb, &
 -3.4758059E-01_jprb, -5.9801831E-01_jprb, -7.1151770E-02_jprb, -4.7359511E-02_jprb, -2.3574675E-03_jprb, -9.7670159E-01_jprb, &
 -9.0150372E-02_jprb,  3.1433590E-01_jprb,  1.7322085E+00_jprb, -3.6608249E-01_jprb, -3.7350106E-02_jprb,  5.0255179E-03_jprb, &
 -2.1117589E-01_jprb,  3.3231080E-01_jprb,  4.9410183E-01_jprb, -5.2753059E-01_jprb,  2.3223385E-01_jprb,  1.0135782E-01_jprb, &
 -2.0056687E-03_jprb,  5.9207518E-01_jprb,  6.9061340E-01_jprb, -1.6416437E-01_jprb, -1.8479596E+00_jprb, -7.8317502E-01_jprb, &
 -3.9250649E-02_jprb,  8.5865597E-03_jprb, -2.0415203E-01_jprb, -8.5867754E-02_jprb, -2.5995129E-01_jprb,  2.7344901E+00_jprb, &
 -1.5906577E-01_jprb, -6.3220574E-03_jprb,  4.1014275E-05_jprb, -5.4942303E-01_jprb,  2.9501305E-01_jprb, -1.7485018E+00_jprb, &
  4.1692016E-01_jprb, -8.3496719E-02_jprb, -6.1155419E-03_jprb, -1.9356531E-02_jprb, -1.9815596E+00_jprb, -1.4933721E-01_jprb, &
 -1.6447586E-01_jprb,  1.5682772E+00_jprb,  2.0103355E-01_jprb, -8.0202789E-02_jprb, -6.6706712E-03_jprb, -1.0235897E+00_jprb, &
  1.6325259E+00_jprb,  2.1112617E-01_jprb, -3.0398529E-01_jprb,  2.5398266E+00_jprb, -7.6298013E-02_jprb,  4.7741383E-03_jprb, &
 -9.2523713E-01_jprb, -9.5724224E-01_jprb, -6.6964851E-01_jprb,  7.2224544E-02_jprb,  1.8133086E+00_jprb,  3.6443040E-02_jprb, &
 -1.1862656E-02_jprb, -1.8971007E-01_jprb, -5.6549421E-01_jprb,  3.1153390E-01_jprb, -8.7217069E-01_jprb, -4.2250202E-01_jprb, &
  1.1567928E-01_jprb, -5.2373420E-04_jprb,  7.2206721E-01_jprb, -3.9936796E-01_jprb,  1.3314477E-02_jprb, -2.4882179E-01_jprb, &
  3.5322208E-02_jprb,  5.8802339E-02_jprb,  3.2360662E-03_jprb, -2.5990725E-01_jprb, -1.2131480E+00_jprb,  6.6448535E-01_jprb, &
 -3.8463100E-01_jprb, -6.4454545E-01_jprb, -2.2140412E-02_jprb,  4.8934508E-03_jprb,  7.4227459E-02_jprb,  5.5076881E-01_jprb, &
  6.1861893E+00_jprb, -1.0429686E+00_jprb, -4.4485456E-01_jprb, -1.1658355E-03_jprb,  4.7340678E-03_jprb,  1.4087819E-01_jprb, &
  4.2677752E-02_jprb,  2.1288956E-01_jprb, -4.2322132E-01_jprb,  5.2667405E-01_jprb,  5.0347901E-02_jprb, -5.2724152E-04_jprb, &
 -1.1831962E-01_jprb,  5.6145512E-01_jprb,  6.9644630E-01_jprb,  5.9816609E-01_jprb, -9.3801611E-01_jprb, -1.7500819E-02_jprb, &
 -5.0187300E-03_jprb, -2.2980368E-01_jprb,  9.4949613E-01_jprb,  2.5243725E-01_jprb, -2.4203677E-01_jprb, -5.7901497E-01_jprb, &
  3.5981000E-02_jprb,  8.2133473E-05_jprb,  2.6699234E-01_jprb,  8.9853015E-01_jprb, -4.3993805E-01_jprb, -8.8496479E-01_jprb, &
 -7.5933640E-01_jprb, -3.5272068E-02_jprb, -4.6216608E-03_jprb, -1.5270172E+00_jprb,  3.6057110E-01_jprb,  1.4627998E-01_jprb, &
 -1.3883531E+00_jprb,  3.0882053E+00_jprb, -4.9536989E-03_jprb,  5.8281051E-03_jprb,  5.4962791E-02_jprb,  2.5453925E-01_jprb, &
  1.6519810E+00_jprb, -7.6930722E-03_jprb, -6.4754072E-01_jprb,  3.0131693E-02_jprb,  6.2731482E-03_jprb,  9.2142592E-01_jprb, &
 -9.5191066E-02_jprb, -1.4135346E+00_jprb, -4.0788537E-01_jprb,  1.1146971E+00_jprb, -2.5350719E-02_jprb,  1.0625266E-03_jprb, &
  2.4353280E-01_jprb, -7.7238217E-01_jprb, -2.2644608E-01_jprb, -1.5925102E+00_jprb, -5.8433845E-01_jprb, -4.5540836E-02_jprb, &
  6.3523611E-03_jprb, -2.0329896E-01_jprb, -3.0489873E-01_jprb,  1.6891196E-01_jprb, -1.3254769E+00_jprb, -6.8700782E-01_jprb, &
 -3.6475152E-02_jprb, -1.0330272E-02_jprb,  4.1739449E-01_jprb, -1.8377843E+00_jprb,  3.1358392E-01_jprb,  1.4968871E+00_jprb, &
 -2.5981205E-01_jprb, -4.2080671E-02_jprb,  6.1793007E-03_jprb, -1.8505454E-01_jprb,  1.9953860E-02_jprb, -3.4228045E-01_jprb, &
 -6.2317902E-02_jprb,  6.4669564E-01_jprb, -2.9002422E-02_jprb,  8.2668899E-04_jprb, -1.7793678E-01_jprb, -5.0671032E-01_jprb, &
  1.2537352E-01_jprb,  7.1101550E-01_jprb,  1.0565185E+00_jprb,  3.1241722E-02_jprb, -8.9982163E-03_jprb, -2.0377190E-01_jprb, &
 -4.7737448E-01_jprb, -8.2188166E-01_jprb,  5.4989705E-02_jprb, -1.1994858E+00_jprb,  8.4837394E-03_jprb, -3.4649026E-03_jprb, &
  6.7188061E-02_jprb,  1.4562030E-01_jprb, -5.7047693E-01_jprb,  8.8102820E-01_jprb, -1.8352287E-01_jprb,  6.0613440E-01_jprb, &
 -1.5635642E-01_jprb,  3.6537211E-01_jprb,  1.2002526E+00_jprb, -4.3333058E+00_jprb, -1.8267243E-01_jprb,  5.6086744E-01_jprb, &
 -1.6144792E-02_jprb,  5.9882989E-04_jprb, -1.2549651E+00_jprb, -1.3567194E+00_jprb,  1.9088658E+00_jprb,  3.9699675E-01_jprb, &
 -1.2425028E+00_jprb, -2.6615628E-02_jprb,  5.2343657E-03_jprb,  3.8868879E-02_jprb,  5.8236447E-01_jprb,  9.1894959E-01_jprb, &
  4.5589957E-01_jprb, -3.1048652E-01_jprb,  6.9992590E-03_jprb,  1.4590729E-03_jprb,  4.4726574E-01_jprb, -1.6747763E-01_jprb, &
  3.3350319E+00_jprb,  4.1423856E-01_jprb, -2.2328622E-02_jprb, -9.4015436E-02_jprb, -9.0557540E-03_jprb, -1.6330006E-02_jprb, &
  9.7183259E-01_jprb, -9.1525090E-02_jprb, -3.9354426E-02_jprb,  3.2775176E+00_jprb,  2.9258834E-03_jprb, -1.4897316E-03_jprb, &
 -4.4617608E-02_jprb, -9.3962129E-02_jprb,  3.3008067E-01_jprb,  9.2493183E+00_jprb,  1.3920845E+01_jprb, -4.7024474E-02_jprb, &
  1.1874635E-02_jprb,  8.1044720E-01_jprb, -1.1693553E+00_jprb,  5.6623050E-03_jprb,  6.7797206E-01_jprb,  1.5649287E-01_jprb, &
  4.8580500E-03_jprb, -1.5225240E-03_jprb, -1.4942454E-01_jprb, -8.7583316E-01_jprb, -1.6485203E-01_jprb,  7.8969649E-01_jprb, &
  3.4673307E-01_jprb,  5.5578158E-03_jprb, -1.9087628E-03_jprb,  1.7779989E-01_jprb, -7.3686726E-01_jprb, -9.6786382E-02_jprb, &
 -1.4952146E+00_jprb, -4.4678739E-01_jprb, -2.0155191E-03_jprb, -2.9307797E-03_jprb,  3.7243457E-01_jprb, -1.3713751E-01_jprb, &
 -3.5495340E-01_jprb, -9.0583286E+00_jprb, -1.4011772E+01_jprb,  5.1885379E-02_jprb, -1.4455505E-02_jprb, -9.6183686E-01_jprb, &
  1.3448964E+00_jprb,  7.2421214E-01_jprb,  9.7908485E-01_jprb,  7.2568280E-01_jprb,  5.9735179E-02_jprb, -1.3473227E-03_jprb, &
  7.9794097E-01_jprb,  1.1123902E-01_jprb, -1.0183905E+00_jprb, -2.7611713E-01_jprb,  2.2446191E-01_jprb,  2.9505459E-02_jprb, &
  9.4470357E-03_jprb, -3.9072049E-01_jprb, -1.1009382E+00_jprb, -2.7722975E+00_jprb, -9.7507820E-01_jprb,  1.0840342E+01_jprb, &
  1.2795248E-02_jprb,  6.4298535E-03_jprb, -1.9742605E-02_jprb,  4.3847440E-01_jprb,  2.9509918E-01_jprb,  5.3657404E-01_jprb, &
  6.0625364E-01_jprb,  7.2662648E-02_jprb, -1.0970171E-02_jprb,  2.4714993E+00_jprb,  2.7183288E+00_jprb, -5.8518622E-01_jprb, &
  8.3419202E-01_jprb,  6.6950443E-01_jprb, -2.2280866E-03_jprb, -4.0782852E-03_jprb, -6.4667607E-02_jprb, -6.3536648E-03_jprb, &
 -4.8352024E-01_jprb, -3.1589389E-01_jprb,  1.3634696E+00_jprb,  2.4006473E-02_jprb, -4.7006776E-03_jprb,  2.1590283E-01_jprb, &
 -8.3513958E-02_jprb,  7.2902612E-01_jprb, -1.4852140E+00_jprb, -3.3273318E-01_jprb, -4.3600398E-02_jprb, -3.6919939E-03_jprb, &
  2.1093352E-01_jprb,  1.6151122E-01_jprb,  3.4538534E+00_jprb, -7.2884067E-01_jprb, -3.4670734E-01_jprb, -1.2534964E-02_jprb, &
 -1.4728539E-02_jprb,  9.3085703E-01_jprb,  1.1670247E+00_jprb,  1.9659030E+00_jprb,  4.2244828E-01_jprb, -1.4221232E+00_jprb, &
 -2.3835994E-02_jprb,  7.7683077E-03_jprb,  1.3006537E-02_jprb,  5.2106202E-01_jprb, -1.0928067E-02_jprb, -1.2042445E+00_jprb, &
  2.3851958E-01_jprb, -1.4245877E-01_jprb,  1.9539081E-02_jprb, -2.3457774E-01_jprb,  6.5719064E-01_jprb, -6.1452125E-01_jprb, &
 -3.5069649E-01_jprb,  1.2840576E+00_jprb,  2.1365621E-02_jprb, -4.1342107E-03_jprb,  2.6900196E-01_jprb, -2.9550839E-01_jprb, &
  5.6717089E+00_jprb, -9.9588079E-01_jprb,  1.2559624E-02_jprb, -2.3061350E-02_jprb, -2.1572142E-03_jprb, -1.1734377E-02_jprb, &
  2.6585966E-01_jprb, -8.7405108E-02_jprb, -1.0633954E+00_jprb,  4.3786716E-01_jprb, -1.3569645E-02_jprb, -1.4197383E-03_jprb, &
 -1.3913697E-01_jprb, -4.7613916E-02_jprb, -1.0228413E+00_jprb, -1.9269256E-01_jprb,  1.5669876E-01_jprb, -1.4283889E-03_jprb, &
  2.7469724E-03_jprb, -1.8849224E-01_jprb,  5.6537562E-02_jprb,  2.9438366E+00_jprb,  1.0262538E+00_jprb, -1.1899882E+01_jprb, &
 -1.4398693E-02_jprb, -6.8767317E-03_jprb,  3.6666529E-02_jprb, -5.7889104E-01_jprb,  1.6856760E-01_jprb, -1.9872358E+00_jprb, &
 -1.2788490E-02_jprb, -3.1412264E-02_jprb,  5.8510152E-03_jprb,  3.6070305E-01_jprb,  4.2627940E-01_jprb,  3.3255310E+00_jprb, &
  1.3533000E-01_jprb, -3.2824009E-02_jprb,  3.0147474E-03_jprb,  4.5471950E-03_jprb,  1.0910965E-01_jprb,  1.4082671E-01_jprb, &
 -9.9663403E-02_jprb, -2.3110319E-02_jprb,  1.2858903E+00_jprb,  4.6771375E-02_jprb, -7.9870529E-03_jprb,  1.2413487E-01_jprb, &
 -2.7720104E-01_jprb, -5.8698135E+00_jprb,  2.2156751E-01_jprb,  3.2193453E-01_jprb,  4.0901304E-03_jprb, -1.1845102E-03_jprb, &
 -1.4098282E-01_jprb, -2.4407969E-01_jprb, -2.2292982E-01_jprb,  5.6937409E-01_jprb, -3.3146530E-01_jprb, -3.7603024E-03_jprb, &
 -1.0021320E-02_jprb, -1.7279446E+00_jprb, -8.5022216E-02_jprb,  3.3293208E-01_jprb, -8.7795495E-01_jprb, -7.0426608E-01_jprb, &
 -1.2144821E-02_jprb,  2.2809429E-03_jprb, -4.3146874E-02_jprb,  3.4653207E-01_jprb,  1.4846534E+00_jprb,  4.5982792E-01_jprb, &
 -3.4226460E-01_jprb, -3.0496983E-02_jprb,  4.7432857E-03_jprb,  1.2658311E-01_jprb,  5.9859299E-02_jprb, -1.3026107E+00_jprb, &
  2.3389392E-01_jprb, -7.9242370E-01_jprb,  1.0888280E-02_jprb,  8.0312714E-03_jprb,  2.5074423E-02_jprb, -5.5077438E-01_jprb, &
 -1.4480577E+00_jprb,  2.0908665E-01_jprb, -1.6386263E-01_jprb, -1.4302659E-02_jprb, -1.5946683E-02_jprb, -1.7324662E+00_jprb, &
 -3.0951188E-01_jprb, -2.7384937E+00_jprb, -5.5727256E-01_jprb,  3.5617868E-01_jprb,  7.0385987E-02_jprb, -1.0846817E-01_jprb, &
 -4.4028502E+00_jprb, -1.0209450E+01_jprb, -1.3516685E-01_jprb, -1.0997933E+00_jprb, -2.1968288E-01_jprb, -9.9249111E-03_jprb, &
 -2.9461392E-03_jprb,  1.0336786E-01_jprb,  1.8323028E-02_jprb, -2.0499075E+00_jprb, -4.0893745E-01_jprb, -5.0919507E-01_jprb, &
 -2.6995807E-02_jprb, -8.5962865E-03_jprb,  1.5069483E-01_jprb, -1.5400682E+00_jprb,  6.2394362E+00_jprb,  1.7463809E-01_jprb, &
 -6.6899636E-01_jprb, -1.7659282E-02_jprb,  4.3106153E-03_jprb,  2.4510703E-01_jprb,  3.9762159E-02_jprb,  1.7649987E+00_jprb, &
 -4.2709524E-01_jprb,  2.4061902E-01_jprb, -2.5569806E-02_jprb, -2.6366725E-03_jprb,  2.2893827E-01_jprb,  4.4948266E-01_jprb, &
  8.1036284E-02_jprb,  6.6120346E-01_jprb,  4.8673647E-01_jprb,  1.1907489E-02_jprb, -5.8074712E-03_jprb, -3.5179059E-01_jprb, &
 -1.1582235E-01_jprb,  1.7721467E+00_jprb,  3.5084223E-01_jprb,  5.9557320E-01_jprb,  5.0787117E-02_jprb,  7.5576582E-03_jprb, &
  4.7945360E-01_jprb,  2.1997527E-01_jprb, -2.5268205E-01_jprb, -4.9820321E-01_jprb,  1.7314307E+00_jprb,  1.2996305E-03_jprb, &
  1.3206292E-03_jprb, -3.2424413E-02_jprb, -5.0259580E-02_jprb, -1.5865574E+00_jprb,  3.9485138E+00_jprb, -3.7467063E-01_jprb, &
 -2.8979900E-02_jprb, -2.0884769E-02_jprb, -1.9068319E+00_jprb, -2.9071358E-01_jprb,  9.7229121E-01_jprb, -1.1895750E-01_jprb, &
  1.3447467E+00_jprb, -4.2408020E-03_jprb,  6.4732986E-03_jprb,  1.8717926E-02_jprb, -2.4255445E-01_jprb,  3.8278900E-01_jprb, &
  8.8663398E+00_jprb,  1.4145540E+01_jprb, -5.7097344E-02_jprb,  1.7392441E-02_jprb,  1.1346957E+00_jprb, -1.5447234E+00_jprb, &
  2.5534870E+00_jprb, -1.0905359E-01_jprb, -7.9613403E-01_jprb, -2.4742381E-02_jprb,  8.4734159E-04_jprb, -1.1796278E-02_jprb, &
  4.9315402E-01_jprb, -5.2674090E+00_jprb, -1.5139754E+00_jprb,  3.1891312E-01_jprb,  1.3701594E-02_jprb,  2.7153047E-03_jprb, &
  5.7470576E-01_jprb, -7.3468957E-01_jprb,  5.8025389E-01_jprb, -1.2728848E-01_jprb,  8.5571359E-01_jprb, -1.0271110E-02_jprb, &
  1.9709157E-03_jprb, -1.1776417E-01_jprb,  6.6009907E-02_jprb,  1.9364057E+00_jprb,  6.7924785E-01_jprb, -1.1314317E-01_jprb, &
 -7.7587520E-02_jprb,  1.0192949E-01_jprb,  4.2664652E+00_jprb,  1.0127927E+01_jprb,  6.9223721E+00_jprb, -3.1539108E-01_jprb, &
 -2.8607853E-01_jprb, -1.5998290E-02_jprb,  1.2387676E-03_jprb,  1.2033192E-01_jprb, -2.2097383E-01_jprb, -1.1967652E+00_jprb, &
  4.9941025E-01_jprb, -2.3160287E-01_jprb, -3.8932675E-02_jprb,  3.4739356E-03_jprb, -8.7088632E-03_jprb, -5.6666844E-01_jprb, &
 -4.6231634E-01_jprb, -4.3034870E-01_jprb,  1.7605224E+00_jprb,  3.0322699E-02_jprb, -4.1439774E-03_jprb,  2.2025263E-01_jprb, &
  1.6742413E-01_jprb, -1.4447185E+00_jprb,  1.4540939E+00_jprb,  7.0417687E-02_jprb, -8.1699924E-02_jprb,  1.2814880E-01_jprb, &
 -2.0638852E+00_jprb, -2.7173949E+00_jprb,  3.4342892E+00_jprb,  3.0283526E-01_jprb, -5.9863862E-01_jprb, -6.4934600E-03_jprb, &
  2.5323418E-03_jprb,  6.6371087E-02_jprb,  3.5204736E-01_jprb,  3.0496813E+00_jprb,  5.1594462E-01_jprb,  1.1926172E-01_jprb, &
 -1.0892746E-01_jprb, -1.0828120E-02_jprb, -3.1402381E-02_jprb,  9.5904865E-01_jprb,  1.5772022E-02_jprb, -8.1938270E-01_jprb, &
  2.1777591E-01_jprb, -1.5856359E-01_jprb,  2.1256497E-02_jprb, -4.5044686E-01_jprb,  8.6282283E-01_jprb, -2.2610083E-01_jprb, &
 -8.4832642E-01_jprb, -4.9911966E-01_jprb, -2.3988030E-02_jprb, -5.8769662E-03_jprb, -9.3947968E-02_jprb, -6.3707563E-01_jprb, &
  5.7800940E-01_jprb, -7.7956164E-01_jprb, -8.4645287E-01_jprb,  2.9740187E-02_jprb,  2.5479206E-03_jprb,  2.2560261E-01_jprb, &
 -1.8559523E-01_jprb, -1.1008306E+00_jprb, -4.3957866E-02_jprb, -3.5222550E-01_jprb, -3.7303091E-02_jprb, -1.5940095E-02_jprb, &
 -1.5577762E+00_jprb, -2.7699878E-01_jprb,  2.2763547E-01_jprb,  4.5150927E-01_jprb, -6.9825639E-01_jprb, -4.7077345E-03_jprb, &
  3.6301860E-03_jprb,  3.3856346E-01_jprb,  3.0022334E-01_jprb,  5.3525187E+00_jprb, -3.4477142E-01_jprb, -2.7404177E-01_jprb, &
  1.3902881E-02_jprb,  3.8404229E-03_jprb,  1.9234626E-01_jprb,  2.6193302E-01_jprb,  2.7892314E+00_jprb,  4.8193162E-01_jprb, &
 -4.0441038E-01_jprb, -6.5838007E-02_jprb,  1.0649285E-01_jprb,  4.4059475E+00_jprb,  1.0144882E+01_jprb,  3.2897177E-02_jprb, &
 -6.2988709E-01_jprb, -2.7318104E-01_jprb,  9.5424971E-03_jprb,  1.4671605E-04_jprb,  1.5481598E-01_jprb,  8.1925584E-01_jprb, &
 -1.0569768E+01_jprb,  1.6043315E-01_jprb,  2.0451487E-01_jprb, -1.2223687E-02_jprb, -3.8529953E-03_jprb,  2.5398907E-02_jprb, &
 -5.2476400E-01_jprb, -9.9602624E-02_jprb, -8.1220423E-02_jprb,  2.8939389E+00_jprb,  1.2940634E-02_jprb, -2.9982005E-03_jprb, &
  8.9178129E-02_jprb, -7.0570392E-02_jprb, -3.0643094E+00_jprb, -5.8961105E-01_jprb, -2.5589000E-01_jprb,  1.2458181E-01_jprb, &
  1.1419512E-02_jprb,  3.5341193E-03_jprb, -8.5045035E-01_jprb, -1.9339978E+00_jprb,  1.2645542E-01_jprb,  2.3992058E-01_jprb, &
 -6.6835922E-01_jprb, -6.6548848E-01_jprb,  2.3842168E-01_jprb,  4.9846459E-01_jprb,  1.5434109E+00_jprb, -9.7134284E-01_jprb, &
  3.3137239E-01_jprb, -1.6962902E-02_jprb,  4.0802164E-02_jprb,  2.2271881E+00_jprb,  1.9473687E+00_jprb, -1.0322639E+01_jprb, &
 -4.8774697E-02_jprb,  4.2511552E-01_jprb,  9.7775667E-03_jprb, -3.7667640E-03_jprb, -9.4116233E-02_jprb,  4.3730664E-02_jprb, &
  3.4553538E+00_jprb,  4.0162394E-02_jprb, -1.6609587E-01_jprb,  6.1417247E-02_jprb,  1.1008723E-02_jprb,  3.9299927E-02_jprb, &
  9.2102024E-01_jprb,  4.3618904E-01_jprb, -9.4792316E-01_jprb,  3.4663745E-01_jprb,  3.4845290E-01_jprb, -1.2887969E-02_jprb, &
  4.8065402E-01_jprb,  1.8419046E+00_jprb,  9.2612563E-02_jprb,  9.1183021E-02_jprb,  1.2179289E+00_jprb,  5.3212237E-02_jprb, &
 -9.0914049E-03_jprb,  2.8403662E-02_jprb, -4.3177264E-01_jprb, -5.1821853E-01_jprb,  4.4530956E-01_jprb,  8.8846662E-01_jprb, &
 -7.2052685E-02_jprb,  3.9385289E-03_jprb,  2.9150314E-02_jprb,  5.7697333E-02_jprb, -3.8567120E+00_jprb,  1.0730456E+00_jprb, &
  3.1759397E-01_jprb, -4.1220797E-02_jprb, -1.6684425E-02_jprb, -3.8361612E-01_jprb,  3.8583490E-01_jprb, -6.2492546E+00_jprb, &
  6.1277093E-01_jprb,  3.5643044E-01_jprb, -3.2197213E-02_jprb, -1.8424504E-03_jprb,  1.4937898E-02_jprb,  3.3211943E-01_jprb, &
  1.4921436E-01_jprb,  1.1380888E+00_jprb, -1.5616567E+00_jprb, -4.2749709E-02_jprb, -1.1946092E-02_jprb, -2.8427960E-01_jprb, &
 -7.7427071E-01_jprb,  1.3416895E+01_jprb, -1.9162935E+00_jprb, -3.3570045E-01_jprb,  2.7045638E-02_jprb,  1.1929815E-02_jprb, &
  4.9987911E-02_jprb, -2.9342764E-01_jprb, -1.4019670E+00_jprb, -5.4688428E-01_jprb,  1.3408320E+00_jprb, -5.8598368E-02_jprb, &
  6.4501436E-03_jprb,  2.4450647E-01_jprb, -9.8151589E-01_jprb, -1.5017836E+00_jprb,  1.1424357E+00_jprb, -1.5579048E+00_jprb, &
  2.1382967E-03_jprb,  2.0894064E-04_jprb, -7.2896361E-02_jprb, -7.3575748E-02_jprb, -6.2943524E+00_jprb,  5.0122477E-02_jprb, &
 -3.2760482E-01_jprb,  7.4826551E-04_jprb, -3.3954832E-03_jprb, -1.8500614E-02_jprb,  1.9750065E-01_jprb, -1.5644496E+00_jprb, &
 -7.4290039E-01_jprb, -8.0836335E-03_jprb,  7.5966204E-02_jprb, -8.8661101E-02_jprb, -4.2611580E+00_jprb, -1.0287330E+01_jprb, &
  5.1258397E+00_jprb, -1.0862163E+00_jprb,  2.7119217E+00_jprb,  3.5596349E-02_jprb,  1.0841311E-02_jprb,  3.9582735E-01_jprb, &
  7.8304477E-01_jprb, -2.9233144E-01_jprb, -2.4946049E-01_jprb, -1.0454179E+00_jprb, -5.7545646E-04_jprb,  1.6380075E-04_jprb, &
 -4.8729744E-01_jprb,  9.7973341E-01_jprb,  3.3856059E+00_jprb,  1.9145871E-01_jprb,  2.0940848E-01_jprb, -1.1100401E-02_jprb, &
  1.9470627E-03_jprb,  9.9879627E-02_jprb,  6.7819845E-02_jprb,  1.6676830E-01_jprb,  1.2058519E+00_jprb, -1.0627885E+00_jprb, &
 -1.1312255E-01_jprb,  2.1944208E-02_jprb,  3.6820707E-01_jprb, -1.7703403E+00_jprb,  8.9669544E+00_jprb, -4.5562427E-02_jprb, &
 -4.0135293E-01_jprb,  1.2178928E-02_jprb,  2.7061957E-03_jprb,  2.2541171E-02_jprb,  3.7095955E-01_jprb,  1.6361637E+00_jprb, &
 -4.2427618E-01_jprb,  4.4177180E-01_jprb, -2.4881177E-02_jprb, -1.5678178E-03_jprb,  3.6737746E-01_jprb,  4.7940038E-01_jprb, &
  9.5194850E+00_jprb,  4.2774126E-02_jprb, -4.2485446E-01_jprb, -2.9645451E-03_jprb,  2.4661474E-03_jprb,  2.1284154E-01_jprb, &
 -1.2709351E-01_jprb,  4.2595515E-01_jprb,  6.3755694E-01_jprb, -1.7752514E-01_jprb,  2.7111160E-02_jprb, -3.5012063E-03_jprb, &
  7.2754277E-01_jprb,  4.5222062E-02_jprb,  1.1056273E+01_jprb, -1.7717747E-01_jprb,  3.8207519E-01_jprb,  5.5663036E-03_jprb, &
  7.8328456E-03_jprb, -3.7461562E-02_jprb,  5.9239614E-01_jprb, -7.5423283E-03_jprb,  1.1732833E+00_jprb, -2.1491057E+00_jprb, &
  4.8931552E-02_jprb, -5.0474600E-03_jprb,  6.2398855E-01_jprb,  7.0518404E-01_jprb, -3.1177752E+00_jprb, -7.6700009E-02_jprb, &
  6.0933012E-01_jprb,  2.1598587E-02_jprb,  8.4656421E-04_jprb,  6.8363340E-01_jprb, -3.1453347E-01_jprb, -1.0952597E+00_jprb, &
  6.0963799E-01_jprb, -3.8229664E-01_jprb, -7.5195960E-02_jprb,  1.4945549E-02_jprb, -1.4568539E-01_jprb, -8.8492540E-01_jprb, &
 -6.7409826E+00_jprb,  5.4373365E-01_jprb,  5.1803289E-01_jprb, -4.6911247E-02_jprb,  4.2342947E-03_jprb,  7.0806030E-03_jprb, &
  6.4789039E-01_jprb, -5.1624595E-01_jprb,  2.0271336E+00_jprb, -3.7693348E-01_jprb, -1.4697737E-02_jprb,  2.5650676E-03_jprb, &
 -2.0502828E-01_jprb, -2.6750354E-01_jprb, -1.6807476E+00_jprb,  1.5460528E+00_jprb,  4.0632445E-01_jprb,  2.8651050E-03_jprb, &
 -1.7065101E-02_jprb, -2.4688866E+00_jprb, -1.0651239E+00_jprb,  1.2868672E-01_jprb, -8.3886060E-01_jprb,  6.9614495E-01_jprb, &
 -3.7124955E-03_jprb,  1.9950471E-02_jprb,  4.8280077E-01_jprb,  1.1619930E+00_jprb,  2.1485584E+01_jprb, -5.5808095E-01_jprb, &
 -3.0405553E-01_jprb, -2.6295814E-03_jprb,  6.6761982E-03_jprb, -1.3230198E-02_jprb, -1.2218652E-03_jprb,  4.0049501E+00_jprb, &
 -4.0893770E-01_jprb, -6.1664502E-01_jprb, -1.4222051E-03_jprb,  1.8716488E-03_jprb,  6.8285690E-02_jprb,  4.7165916E-01_jprb, &
 -6.1436018E+00_jprb,  2.1324722E-01_jprb,  4.8455649E-01_jprb,  1.0263915E-02_jprb, -2.5088297E-04_jprb, -1.0086594E-01_jprb, &
  1.4415640E-01_jprb, -1.1035755E+00_jprb,  1.8685458E-02_jprb, -1.1987114E+00_jprb, -1.0426665E-02_jprb,  3.1721047E-04_jprb, &
 -4.9120001E-02_jprb,  3.3831895E-01_jprb,  1.3184448E+00_jprb,  1.1103819E-02_jprb,  1.9585301E+00_jprb, -1.4692981E-02_jprb, &
 -5.5344003E-04_jprb, -2.3123311E-01_jprb, -1.3907855E-01_jprb, -5.0468373E+00_jprb, -1.4711494E+00_jprb,  2.8534176E-01_jprb, &
  1.7639291E-02_jprb,  2.3135365E-03_jprb,  2.6186543E-01_jprb, -3.8059194E-01_jprb, -7.7612762E-01_jprb,  8.1410334E-02_jprb, &
 -5.2924027E+00_jprb, -8.6117847E-03_jprb, -9.6643934E-03_jprb, -6.9793380E-01_jprb,  7.7651371E-02_jprb,  1.7924276E+00_jprb, &
  3.4609152E-01_jprb, -2.2386655E-01_jprb, -1.2609197E-01_jprb,  2.2311891E-02_jprb,  2.4199463E-01_jprb, -3.2626452E-01_jprb, &
  1.8132694E+01_jprb, -4.5146183E-01_jprb, -3.3475570E-01_jprb, -3.9589763E-03_jprb,  3.3711924E-03_jprb, -1.4863433E-03_jprb, &
 -4.1266582E-02_jprb,  6.7874162E-01_jprb, -2.1136114E+00_jprb, -8.1855552E-01_jprb, -2.8540591E-01_jprb,  9.8956857E-02_jprb, &
 -5.1936941E-01_jprb,  1.7833453E+00_jprb, -1.3901899E-01_jprb,  2.0408873E+00_jprb,  3.6183858E-01_jprb,  1.7812821E-02_jprb, &
 -4.0562379E-03_jprb, -4.4874367E-01_jprb, -1.0155032E+00_jprb, -6.6442184E+00_jprb, -9.0677391E-02_jprb,  4.4708766E-01_jprb, &
 -2.2770800E-02_jprb, -3.6177709E-04_jprb, -5.7697233E-01_jprb,  3.0174568E-01_jprb,  2.1741024E+00_jprb,  1.1606564E+00_jprb, &
  1.6326010E-01_jprb,  1.7620961E-02_jprb, -3.4166061E-02_jprb,  2.1748681E+00_jprb,  3.9327744E+00_jprb,  2.3133801E-01_jprb, &
 -1.8476848E+00_jprb,  1.1248976E+01_jprb, -6.6580815E-03_jprb,  5.8145610E-03_jprb,  2.3411861E-01_jprb,  1.1481861E-01_jprb, &
 -3.0645908E-01_jprb,  1.7839102E+00_jprb,  5.1395075E-01_jprb,  1.6980614E-01_jprb, -7.8907915E-04_jprb, -1.2943711E+00_jprb, &
  2.6529923E+00_jprb,  5.7738807E+00_jprb, -1.0525285E+00_jprb, -3.5962750E-01_jprb,  1.8503179E-03_jprb,  4.6964609E-03_jprb, &
  1.3983248E-01_jprb,  1.0083350E-01_jprb, -3.1907853E+00_jprb,  1.2994350E+00_jprb,  5.5953439E-01_jprb,  4.3119008E-02_jprb, &
  9.1952340E-05_jprb,  5.6738849E-01_jprb,  3.5811843E-01_jprb, -6.6549833E+00_jprb,  7.7352017E-02_jprb, -4.5651096E-01_jprb, &
  1.2274670E-02_jprb, -4.7222544E-03_jprb, -3.5837316E-02_jprb,  6.3638227E-01_jprb, -1.2403861E+01_jprb,  3.1261199E-01_jprb, &
 -2.5438453E-01_jprb,  6.8744326E-03_jprb, -7.1437247E-03_jprb, -6.1324986E-03_jprb, -9.1617615E-02_jprb, -2.6574799E-01_jprb, &
  1.6075732E+00_jprb, -1.0216848E+01_jprb,  2.3156481E-03_jprb, -4.2162684E-03_jprb, -1.3920002E-01_jprb, -2.2077556E-01_jprb /

! bias vector of first layer, vector [no * 1]
! 8 elements in b2
DATA aniso_net%b2 / -7.9389554e-01_jprb, -1.6530464e+00_jprb, 1.2170673e+00_jprb,  7.7658185e-01_jprb, &
                    -1.2672478e+00_jprb,  2.7329894e-01_jprb, 2.9603450e+00_jprb, -4.1883279e+00_jprb  /

! weight matrix of first layer, matrix [nh * no]
! 180 * 8 elements in W2
DATA aniso_net%W2 / &
 4.4859959E-01_jprb,  2.8223066E+00_jprb,  4.0482956E-02_jprb,  6.9912169E-02_jprb, -7.4536136E-02_jprb,  4.3291102E-03_jprb, &
 3.4419222E-01_jprb, -3.6709312E-01_jprb,  7.6206071E-01_jprb, -9.1072440E+00_jprb, -9.4066498E-02_jprb,  1.4529595E-01_jprb, &
-4.5819107E-02_jprb, -2.6382870E-01_jprb, -4.1208972E-01_jprb, -9.7060558E-02_jprb,  7.1764224E-01_jprb, -1.8719611E-01_jprb, &
 5.6760623E-01_jprb,  4.1845140E+00_jprb,  6.3506404E-03_jprb,  2.4362293E-01_jprb, -4.6957685E-02_jprb, -5.0472708E-01_jprb, &
-1.1580961E-03_jprb, -5.9554233E-01_jprb, -6.3566296E-03_jprb,  2.4752111E-03_jprb,  9.6077315E-01_jprb, -7.1742477E-01_jprb, &
-2.9049034E-02_jprb, -6.3915845E-01_jprb, -1.4609030E-01_jprb,  1.0038827E+00_jprb,  1.2295745E+00_jprb,  1.1616559E+00_jprb, &
-7.7916158E-01_jprb,  1.3421996E+00_jprb,  1.0179922E-02_jprb, -3.4462121E-02_jprb,  9.2488883E-01_jprb,  1.0581681E+00_jprb, &
 1.4420446E+00_jprb, -1.8922641E-02_jprb, -2.5541064E-02_jprb,  4.7695137E-01_jprb, -2.1147473E-01_jprb, -2.9813380E+00_jprb, &
 2.7884462E+00_jprb, -2.1413481E-01_jprb,  1.7312873E-01_jprb, -1.5926585E-01_jprb,  1.4065148E-01_jprb, -3.4423136E-01_jprb, &
 5.4766021E-02_jprb, -4.1853211E-01_jprb, -1.2052211E+00_jprb,  5.5363283E-02_jprb, -2.1805417E+00_jprb, -1.5249846E+00_jprb, &
 9.4318370E-02_jprb, -1.2843690E+00_jprb,  2.1143474E-02_jprb, -2.2319172E-02_jprb,  2.7946473E-01_jprb, -1.2521918E+00_jprb, &
-6.5496034E-01_jprb, -1.8938612E-01_jprb,  3.1558388E-02_jprb, -1.4506692E+00_jprb,  1.0754242E-01_jprb, -1.4243059E-01_jprb, &
 6.6021469E-02_jprb,  4.2364129E-02_jprb,  1.9551465E-01_jprb, -8.2335775E-02_jprb,  9.2240256E-03_jprb,  7.0863092E-02_jprb, &
-5.8460816E-01_jprb,  5.8049987E-02_jprb, -5.0906801E-02_jprb, -2.2555567E-01_jprb, -3.4380123E-01_jprb,  1.3733513E-01_jprb, &
-2.0327495E+00_jprb,  4.3353306E-01_jprb, -3.6263492E-01_jprb, -3.3749763E-02_jprb,  5.6896203E-01_jprb,  3.9164291E-01_jprb, &
-8.4013131E-02_jprb,  2.6077481E+00_jprb, -1.6207791E-01_jprb, -5.0525727E-01_jprb,  1.4339459E+00_jprb,  2.3199040E-02_jprb, &
-1.8410989E+00_jprb, -3.5576581E-01_jprb, -2.2679005E+00_jprb,  2.6151035E-01_jprb, -3.2893746E-01_jprb, -3.7038837E-01_jprb, &
-1.2648387E-01_jprb,  9.3144113E-02_jprb, -1.5828015E-01_jprb,  1.4862753E-01_jprb, -4.8395144E-01_jprb,  3.5241337E-02_jprb, &
 8.1471115E-01_jprb, -1.7989250E+00_jprb, -1.2391650E+00_jprb, -6.2210837E-01_jprb,  5.3266183E-01_jprb,  2.7883541E-01_jprb, &
 2.2904268E-01_jprb, -2.0547484E-02_jprb, -9.7349800E-01_jprb,  5.7335716E-01_jprb,  4.1763661E-01_jprb,  2.5054154E-01_jprb, &
 9.1282856E-02_jprb,  7.5684156E-01_jprb, -7.4808255E-01_jprb,  1.1803589E-02_jprb, -3.6428217E-01_jprb, -8.5162682E-01_jprb, &
 2.6755175E+00_jprb,  2.1941089E-01_jprb,  1.3087809E-01_jprb,  1.3848450E-03_jprb,  5.7558655E-02_jprb,  6.2006357E+00_jprb, &
-3.5423275E-01_jprb,  3.0550403E-02_jprb,  1.7272999E-01_jprb,  9.9050395E-02_jprb,  2.4385990E+00_jprb,  3.8012501E-01_jprb, &
 1.8701211E-02_jprb,  4.3859632E+00_jprb,  1.5921024E-01_jprb,  1.9065461E-01_jprb, -6.0199661E+00_jprb,  3.4738889E-01_jprb, &
 7.6250719E-02_jprb, -1.0137672E-01_jprb,  6.7396838E-01_jprb,  1.9608029E-01_jprb,  8.5224224E-01_jprb,  4.5031147E-01_jprb, &
 6.6080776E+00_jprb, -1.2156399E+00_jprb, -1.9742957E+00_jprb, -3.3935201E-02_jprb, -1.9968384E-04_jprb,  3.4625465E-02_jprb, &
 4.9618914E-01_jprb, -9.0441498E-01_jprb,  9.4334845E-02_jprb,  3.4512551E-02_jprb, -3.7576859E+00_jprb,  1.6319972E+00_jprb, &
 5.5617806E-01_jprb,  9.9099640E-02_jprb,  4.9869579E-01_jprb,  1.9756745E+00_jprb, -4.4094402E-02_jprb, -1.9744616E-01_jprb, &
 5.4056657E+00_jprb,  4.9967483E-02_jprb, -4.9261805E-01_jprb, -8.3716778E-01_jprb,  2.8599871E-02_jprb,  6.9523882E-02_jprb, &
-1.1571101E-01_jprb,  3.1551422E+00_jprb, -1.8776460E-01_jprb,  4.3839655E+00_jprb, -8.7321812E+00_jprb,  1.0447575E-01_jprb, &
-1.6078390E+00_jprb,  9.2029655E-01_jprb, -7.5062398E-02_jprb, -5.2246315E-02_jprb, -8.8482743E-01_jprb,  1.9558597E-02_jprb, &
 8.3185139E-01_jprb, -1.0258824E-02_jprb, -2.0181898E-01_jprb, -7.7213422E-01_jprb,  1.4395750E-01_jprb, -3.2952800E-01_jprb, &
 5.9904947E-02_jprb,  2.8924845E-01_jprb,  1.5848230E-01_jprb, -2.1202967E-01_jprb,  9.0003354E-02_jprb, -1.0622328E-01_jprb, &
-8.4896976E-02_jprb,  4.1589686E-01_jprb,  4.1626947E-03_jprb,  2.1292052E-01_jprb, -5.5308707E-03_jprb,  6.7441792E-01_jprb, &
-5.1714664E-02_jprb, -3.0455141E-01_jprb, -1.1916887E-02_jprb,  1.7538881E-02_jprb, -1.8069001E-01_jprb, -9.4065899E-02_jprb, &
 2.6932302E-02_jprb, -9.7404402E-03_jprb,  3.4477905E-01_jprb,  4.7276853E-01_jprb, -6.3111929E-01_jprb, -2.0312542E+00_jprb, &
-1.1402331E-01_jprb,  2.0693920E-01_jprb,  1.0893809E-01_jprb, -1.2473281E-01_jprb, -4.6353449E-01_jprb,  4.4626404E-02_jprb, &
-1.1096911E-02_jprb, -6.3697182E-06_jprb, -9.1575422E-03_jprb,  1.4134330E-01_jprb, -5.2092260E-02_jprb, -1.5019193E+00_jprb, &
-2.8668857E+00_jprb, -1.8804857E-01_jprb,  1.6304726E-03_jprb, -2.5505838E-01_jprb,  3.4579924E-02_jprb,  1.1470591E-01_jprb, &
-2.9302956E-01_jprb, -4.8391206E-02_jprb,  5.1675446E-01_jprb, -9.8902073E-02_jprb, -2.4016267E-01_jprb, -2.7097559E-01_jprb, &
 8.4562918E-02_jprb,  3.8782829E-01_jprb,  5.7508455E-03_jprb,  2.4083206E-03_jprb, -5.4644297E-02_jprb,  3.0446923E-01_jprb, &
-7.7795732E-01_jprb,  7.6699168E-02_jprb,  1.8057545E-02_jprb,  3.5419875E-01_jprb,  1.9782665E-01_jprb, -1.8172036E-01_jprb, &
 3.3386690E-02_jprb,  5.5192740E-02_jprb,  5.2445912E-03_jprb,  1.5024763E-01_jprb, -3.8806192E-02_jprb, -7.6541173E-02_jprb, &
 3.3410020E-01_jprb, -1.0509918E+00_jprb,  1.9329794E-02_jprb,  1.9425021E-02_jprb,  1.2845824E-02_jprb, -1.0427465E-01_jprb, &
-7.0509644E-02_jprb,  3.0324748E-01_jprb,  2.7651224E-01_jprb,  6.6459015E-02_jprb,  2.0785180E-01_jprb,  6.0966185E-01_jprb, &
-1.4139250E-01_jprb, -4.3336346E-01_jprb,  6.5592542E-02_jprb,  6.6968600E-03_jprb,  3.7149355E-01_jprb, -1.3035553E-01_jprb, &
-1.6981536E-01_jprb,  9.2409103E-02_jprb, -4.5199473E-01_jprb,  9.7359698E-02_jprb, -8.7271248E-02_jprb,  1.6523602E-01_jprb, &
-4.7467723E-01_jprb, -8.3611736E-02_jprb,  1.5218986E-01_jprb,  4.9504699E-02_jprb,  1.9270637E-01_jprb,  1.4925888E-02_jprb, &
-1.0594934E-01_jprb, -3.3338333E-01_jprb,  2.6081065E-01_jprb,  1.2603156E-01_jprb, -2.3500338E-02_jprb, -5.2031609E-01_jprb, &
-1.1194131E-01_jprb, -2.0491208E-03_jprb, -8.2344120E-02_jprb,  1.6188966E+00_jprb,  2.3585660E-02_jprb,  4.4880339E-01_jprb, &
 3.3156548E-02_jprb,  1.4184302E-01_jprb, -2.8455628E-01_jprb, -1.0645731E+00_jprb,  8.5979357E-02_jprb,  9.2646858E-01_jprb, &
 7.0809351E-01_jprb, -1.8300333E-01_jprb,  7.5146608E-01_jprb,  9.5324221E-04_jprb,  2.2555211E-02_jprb,  4.8144803E-01_jprb, &
 2.6334702E-01_jprb,  1.2694432E-02_jprb,  1.1181519E-01_jprb,  3.0207331E-02_jprb,  3.1829210E-01_jprb,  1.0055247E-01_jprb, &
 5.2388029E-03_jprb,  3.6599894E-01_jprb,  5.8245799E-02_jprb,  8.2006936E-02_jprb, -4.9196032E-01_jprb, -6.2984114E-02_jprb, &
-1.1150492E-02_jprb,  4.3123716E-02_jprb, -5.4220069E-01_jprb, -1.5405213E-01_jprb,  1.6367514E+00_jprb, -2.1422605E-01_jprb, &
 4.3789259E-01_jprb,  4.6704546E-02_jprb, -3.5367862E-01_jprb, -1.1736987E-02_jprb,  4.7009247E-02_jprb,  1.7070330E-01_jprb, &
-3.9850663E-01_jprb,  8.1322668E-02_jprb, -1.0418382E-02_jprb,  1.5904179E-02_jprb, -1.2014929E+00_jprb, -2.9332860E-01_jprb, &
-2.5986642E-01_jprb, -2.9515123E-01_jprb, -1.4057544E-01_jprb,  3.7618521E-01_jprb,  2.2149217E-02_jprb, -4.7019404E-02_jprb, &
 2.6075644E+00_jprb, -2.9045953E-03_jprb,  2.8820522E-01_jprb,  7.8588313E-01_jprb,  2.3420699E-02_jprb, -3.8757459E-02_jprb, &
-1.1642744E-01_jprb,  3.8973345E+00_jprb, -5.0161837E-02_jprb,  6.4132681E-01_jprb, -1.4679000E+00_jprb, -5.4329757E-02_jprb, &
 1.3913050E+00_jprb,  2.5554502E-01_jprb,  2.3341782E-01_jprb,  4.3311067E-02_jprb, -1.6517203E-02_jprb, -1.2661819E-02_jprb, &
-1.5282787E-01_jprb,  2.4044829E-01_jprb, -1.9024182E-01_jprb, -1.8834697E+00_jprb, -2.7189183E-01_jprb,  2.0369499E-02_jprb, &
-1.1898793E-01_jprb,  1.3338586E-01_jprb,  1.4421315E-01_jprb, -3.9761860E-01_jprb, -4.4383021E-01_jprb,  1.1152781E-01_jprb, &
 3.4890648E-01_jprb,  1.4791199E+00_jprb, -9.0323058E-03_jprb,  6.0239746E-02_jprb, -2.4240699E-02_jprb, -6.8538777E-02_jprb, &
 1.0287288E-01_jprb, -3.6340488E-01_jprb,  3.0616496E-03_jprb, -1.3827566E-02_jprb, -1.8769514E-01_jprb,  4.3191361E-01_jprb, &
-2.7013418E-02_jprb, -6.0607900E-02_jprb, -1.4406452E+00_jprb,  1.0763897E-01_jprb, -3.4484802E-02_jprb,  1.5452227E+00_jprb, &
 4.9172440E-01_jprb,  4.6080759E-01_jprb,  2.6912879E-01_jprb, -1.4220075E-01_jprb,  1.4773013E-02_jprb, -6.5637890E-02_jprb, &
 7.6715154E-01_jprb, -3.0271624E-02_jprb, -5.5590710E-03_jprb,  7.3253741E-02_jprb, -9.5816347E-02_jprb, -7.4851202E-01_jprb, &
 1.1127621E-01_jprb, -3.7191460E-01_jprb,  5.0516238E-02_jprb,  3.7163015E-01_jprb, -7.4629913E-03_jprb,  7.4808686E-02_jprb, &
-4.8577531E-01_jprb, -3.0636014E-02_jprb,  1.1542839E-01_jprb, -6.2942249E-02_jprb, -4.4405618E-01_jprb,  8.2086325E-01_jprb, &
 2.3477348E-01_jprb, -4.3982648E-01_jprb,  5.4955114E-03_jprb, -3.7512150E-02_jprb, -5.6487621E-02_jprb,  3.0994030E-02_jprb, &
 3.1549048E-01_jprb, -9.1590303E-02_jprb,  2.1867303E-03_jprb,  3.9348909E-01_jprb, -1.6453761E-01_jprb,  3.3956151E-02_jprb, &
 1.9994654E-02_jprb, -7.0480707E-02_jprb,  3.0182528E-01_jprb, -2.7011379E-01_jprb,  7.7588316E-03_jprb, -7.1845236E-01_jprb, &
-4.0293016E-01_jprb,  8.7837102E-02_jprb, -8.8366528E-02_jprb,  9.9550051E-02_jprb, -1.4507226E-01_jprb,  2.5631487E-01_jprb, &
-6.1708097E-01_jprb, -3.6633421E-01_jprb, -5.7323377E-02_jprb, -1.1650805E-01_jprb, -1.3683936E-02_jprb, -2.0308880E+00_jprb, &
 5.3534470E-02_jprb, -2.5512544E-01_jprb, -4.8603432E-02_jprb, -7.8075020E-01_jprb, -1.9885653E+00_jprb, -2.7004452E-02_jprb, &
 2.0018232E-01_jprb, -2.7188339E-02_jprb,  8.4937753E-01_jprb, -2.7864117E-01_jprb,  4.8818731E-02_jprb, -3.7423990E-01_jprb, &
-2.0438579E-01_jprb,  1.6006033E-01_jprb, -5.4947514E-02_jprb,  2.1464813E-03_jprb, -1.8538934E-01_jprb,  1.8368469E-02_jprb, &
 7.7578477E-02_jprb, -6.4476129E-04_jprb, -3.0628669E-01_jprb,  7.4981472E-02_jprb, -8.7632304E-01_jprb,  1.5489848E-01_jprb, &
 1.2945158E-01_jprb, -1.2570770E-02_jprb,  6.0247436E-01_jprb, -1.2863085E+00_jprb,  1.1574314E-01_jprb, -2.2087976E-01_jprb, &
 1.1985682E-01_jprb, -1.7892066E-01_jprb,  7.6540466E-01_jprb, -1.3859674E-01_jprb, -1.7272802E-02_jprb,  3.4194251E-01_jprb, &
 3.1478368E-01_jprb,  1.0382686E-01_jprb, -8.5963633E-01_jprb,  7.1542715E-04_jprb,  2.4055425E-03_jprb, -3.0805576E-01_jprb, &
-7.6507371E-02_jprb, -1.4160560E-02_jprb, -1.0071741E-01_jprb,  4.2699031E-01_jprb,  7.9631408E-02_jprb, -1.6660133E+00_jprb, &
-2.3986257E-02_jprb,  8.8154703E-01_jprb, -1.4080886E-02_jprb, -4.0887982E-02_jprb, -1.5345795E+00_jprb, -1.1157097E-01_jprb, &
 4.3380203E-02_jprb, -2.6444288E-01_jprb,  1.4961022E+00_jprb, -9.2665367E-01_jprb, -2.7363473E+00_jprb,  3.6683249E-01_jprb, &
-5.9653465E-01_jprb,  4.1956684E-01_jprb,  8.5587402E-02_jprb, -5.6988603E-02_jprb, -1.8823488E-01_jprb, -1.2711785E-01_jprb, &
 2.1984369E+00_jprb, -5.0319907E-01_jprb,  3.5302751E-02_jprb, -4.4919871E-02_jprb, -7.9019459E-01_jprb,  2.9024430E-01_jprb, &
-6.4934354E-01_jprb,  3.7055084E-01_jprb,  2.2596814E-01_jprb, -2.1560387E-02_jprb, -1.3422798E-02_jprb,  4.5108569E-01_jprb, &
 2.9271183E+00_jprb,  1.5879095E-02_jprb, -4.0952236E-01_jprb, -3.3176095E+00_jprb, -1.5063379E-02_jprb, -1.0298582E-03_jprb, &
 4.8374160E-02_jprb,  6.5030694E-01_jprb,  2.8762866E-01_jprb,  8.5506029E-01_jprb, -8.2457466E-01_jprb, -2.5480933E-03_jprb, &
 3.6001111E-01_jprb,  4.2678821E-01_jprb, -3.1801909E-01_jprb, -1.9536704E-02_jprb,  7.1093998E-02_jprb,  1.7482733E-02_jprb, &
 6.4936350E-01_jprb,  2.1975109E-01_jprb, -3.6859265E-01_jprb,  1.6860542E+00_jprb, -2.1656778E-02_jprb, -2.4961925E-01_jprb, &
-2.3642927E-02_jprb,  3.4526025E-01_jprb,  2.1531834E-01_jprb, -9.7765379E-01_jprb, -3.1410641E-01_jprb,  3.1416108E-02_jprb, &
 2.0770769E-01_jprb, -4.1238301E-01_jprb, -1.4261739E-02_jprb, -5.2956247E-02_jprb,  2.1276132E-02_jprb, -9.7331833E-02_jprb, &
-3.3693691E-02_jprb,  5.8467162E-01_jprb, -9.3877842E-02_jprb,  5.8109066E-03_jprb, -1.1071666E-01_jprb,  4.9561174E-01_jprb, &
-1.0732947E-02_jprb,  2.4061525E-01_jprb, -6.7574913E-01_jprb, -5.4334795E-01_jprb, -4.7919151E-01_jprb, -1.8704158E-01_jprb, &
 3.3547886E-01_jprb, -6.1173867E-01_jprb, -1.1685929E-02_jprb,  2.6691023E-01_jprb,  2.6293102E-02_jprb, -1.9567929E-01_jprb, &
-1.7079512E-01_jprb,  7.5225197E-03_jprb,  2.8678251E-02_jprb,  8.0624420E-02_jprb,  6.8286503E-02_jprb, -2.2743302E-01_jprb, &
-1.1222872E+00_jprb, -2.7346350E-01_jprb,  8.2628320E-02_jprb, -2.2204647E-01_jprb, -4.2724836E-02_jprb,  9.2749607E-02_jprb, &
-3.0119383E-01_jprb,  1.4717586E-02_jprb,  4.7868693E-01_jprb, -4.5533873E-02_jprb,  9.2012920E-01_jprb,  6.3595422E-01_jprb, &
 4.2818925E-01_jprb,  1.6869475E-01_jprb,  2.8713109E-02_jprb,  2.7753527E-03_jprb, -1.0746538E-01_jprb, -6.4435578E-02_jprb, &
 1.4819279E-01_jprb, -1.4977711E-01_jprb,  2.5657357E-02_jprb, -1.8995190E-01_jprb, -2.7453579E-01_jprb, -5.6961173E-02_jprb, &
 3.4904743E-02_jprb, -1.2786091E-02_jprb, -1.3973940E-01_jprb,  5.2147132E-02_jprb,  3.1488882E-02_jprb, -1.3940611E-01_jprb, &
-6.0386814E-02_jprb, -6.2645979E-01_jprb,  2.0444153E-02_jprb, -8.9125524E-03_jprb, -4.3695473E-02_jprb,  7.2241515E-02_jprb, &
 2.7067389E-01_jprb,  3.2703752E-01_jprb, -1.3451186E-01_jprb,  4.9580370E-02_jprb, -4.5126907E-03_jprb, -1.2871867E-01_jprb, &
 1.4303622E-01_jprb, -7.5716216E-02_jprb, -6.3649546E-02_jprb,  2.4458090E-01_jprb, -3.8347491E-01_jprb,  7.5229029E-03_jprb, &
 2.6077015E-01_jprb,  5.8455513E-02_jprb,  6.9630582E-01_jprb, -2.4537714E-01_jprb,  3.4784620E-01_jprb,  3.4205150E-02_jprb, &
-7.3578396E-01_jprb, -1.7924453E-01_jprb, -1.0706486E-02_jprb, -8.2269173E-02_jprb,  1.0030927E-01_jprb,  9.5654444E-03_jprb, &
-2.6609673E-01_jprb, -1.3857576E-01_jprb,  7.9908923E-02_jprb, -2.8805362E-01_jprb, -8.2390334E-01_jprb,  5.5930343E-01_jprb, &
-4.9715293E-02_jprb,  3.6335733E-02_jprb,  4.2505378E-01_jprb, -1.4473619E-01_jprb, -2.1175644E-01_jprb, -9.9326002E-02_jprb, &
-8.7239391E-03_jprb, -3.7648207E-02_jprb,  4.1876647E-01_jprb,  6.6780517E-01_jprb,  2.5456612E-01_jprb, -4.6087316E-01_jprb, &
 6.7583350E-01_jprb,  7.9813756E-02_jprb,  1.0967540E-01_jprb, -1.8381670E-02_jprb, -6.7145631E-02_jprb, -1.0557886E+00_jprb, &
-4.9807523E-01_jprb,  1.3816038E-02_jprb,  3.7548723E-02_jprb,  4.1819679E-01_jprb, -6.6629675E-02_jprb, -6.6850735E-01_jprb, &
 6.0487589E-02_jprb, -1.7022166E-01_jprb, -1.2182799E-01_jprb,  9.0347354E-02_jprb,  5.5730087E-01_jprb, -2.2155943E-01_jprb, &
-1.1097164E-02_jprb,  5.1674004E-02_jprb, -3.9152652E-02_jprb, -7.7146785E-01_jprb,  1.0073932E-01_jprb,  2.8261658E-01_jprb, &
-7.1511806E-01_jprb,  5.3212574E-01_jprb, -3.6578809E-02_jprb,  4.4895418E-02_jprb,  1.3918759E-01_jprb,  3.1110425E-02_jprb, &
 5.2848658E-01_jprb, -3.1927556E-01_jprb,  2.3300588E-02_jprb,  6.8770770E-02_jprb, -3.7887972E-01_jprb,  1.3436440E-01_jprb, &
-8.8375172E-01_jprb, -1.3086007E-01_jprb,  8.4649907E-02_jprb,  2.7234190E-01_jprb, -4.3411940E-03_jprb, -6.0978514E-01_jprb, &
-1.2180735E+00_jprb, -8.6318821E-02_jprb,  4.2675441E-01_jprb, -4.7338873E-02_jprb, -1.0028958E-02_jprb,  9.4878467E-02_jprb, &
-7.2885637E-02_jprb, -6.4719313E-01_jprb, -8.4685109E-02_jprb, -2.2826137E-01_jprb, -5.7344535E-01_jprb,  1.1208573E-01_jprb, &
-2.5692224E-01_jprb,  1.5125199E+00_jprb,  8.0575213E-02_jprb, -1.6551629E+00_jprb,  3.1161446E-01_jprb,  1.8264800E+00_jprb, &
 6.0387086E-01_jprb, -3.4014195E-01_jprb, -1.8076474E+00_jprb, -2.2045395E+00_jprb, -1.9345955E+00_jprb,  2.0289087E+00_jprb, &
-6.2705823E-01_jprb, -2.0150464E+00_jprb, -1.9858645E+00_jprb, -2.7103905E-01_jprb, -7.1785085E-01_jprb, -2.6810021E-01_jprb, &
-1.9532234E-01_jprb,  1.5488287E+00_jprb, -2.5903210E-02_jprb, -2.0915012E-01_jprb,  4.1901859E-02_jprb, -5.6361463E-01_jprb, &
 5.1088006E-02_jprb, -1.0672438E+00_jprb, -7.1009806E-01_jprb,  1.3656456E+00_jprb, -1.8461511E+00_jprb, -1.0559280E+00_jprb, &
 2.8153243E-01_jprb,  2.9571038E-01_jprb,  1.9127826E+00_jprb, -3.3792373E-01_jprb, -1.6487095E+00_jprb,  8.9288347E-02_jprb, &
 3.8626664E-01_jprb,  6.4899199E-01_jprb,  1.5368874E-02_jprb, -4.6983484E-01_jprb,  7.7518207E-01_jprb, -5.0578446E-01_jprb, &
 1.3858873E+00_jprb,  1.5463883E-01_jprb,  2.2296385E-02_jprb,  4.5024642E-01_jprb, -9.0830633E-02_jprb, -7.9061143E-01_jprb, &
 3.1862625E+00_jprb, -3.9747900E-01_jprb,  5.1634941E-01_jprb,  1.3487753E-01_jprb,  4.2037343E-02_jprb, -1.9017827E+00_jprb, &
-1.8862159E+00_jprb, -1.3849402E+00_jprb,  2.2187227E+00_jprb, -6.7799153E-02_jprb,  9.4345518E-03_jprb, -7.1834575E-01_jprb, &
 9.5469588E-01_jprb,  4.4799622E-01_jprb,  1.0496194E-02_jprb, -7.1303580E-02_jprb,  9.6979715E-01_jprb, -2.5024181E-02_jprb, &
-1.6753796E+00_jprb,  2.6918309E+00_jprb, -2.5014077E+00_jprb,  4.2500810E-01_jprb, -1.2724043E+00_jprb, -2.7736589E-02_jprb, &
-4.7133990E+00_jprb, -1.5460214E-01_jprb,  2.5116631E-02_jprb, -2.4178884E+00_jprb,  4.1371843E-02_jprb, -1.5864897E+00_jprb, &
-4.0786715E-02_jprb, -7.2829677E-01_jprb, -8.0132520E-02_jprb, -8.5011579E-01_jprb, -2.7434567E-01_jprb, -3.5720859E-01_jprb, &
-1.4227762E+00_jprb,  1.2592639E+00_jprb,  6.9852818E-01_jprb, -1.1169920E+00_jprb,  6.5205363E-01_jprb,  3.0644056E-01_jprb, &
-1.5402368E+00_jprb, -2.6669423E+00_jprb, -1.9486222E-01_jprb, -2.6189143E+00_jprb,  1.5501575E+00_jprb, -2.2477646E-01_jprb, &
 1.0923901E+00_jprb,  1.0729543E+00_jprb, -1.6686655E-01_jprb,  1.6925443E-02_jprb, -5.9043101E-01_jprb, -2.5079180E+00_jprb, &
-1.1581953E+00_jprb, -2.1296759E-01_jprb, -3.0564101E-01_jprb, -4.0643122E-02_jprb,  4.7601026E-01_jprb, -2.2129389E+00_jprb, &
 1.2854432E+00_jprb, -8.2902737E-01_jprb, -1.6735483E+00_jprb,  1.6707779E+00_jprb, -4.5211560E-01_jprb, -1.2148039E+00_jprb, &
 2.3763770E-01_jprb, -9.0192816E-03_jprb, -1.2452570E+00_jprb,  2.6483510E+00_jprb,  1.1634415E-01_jprb, -2.7894301E-01_jprb, &
-8.2403459E-01_jprb, -4.5590351E-01_jprb, -6.4555135E-01_jprb, -2.5587565E+00_jprb,  9.0318602E-01_jprb,  9.2590753E-02_jprb, &
 4.7648231E-01_jprb, -4.1785142E+00_jprb,  1.1421816E+00_jprb,  6.5282989E-04_jprb,  1.3260414E-02_jprb,  2.2673800E+00_jprb, &
-2.6744852E-02_jprb, -1.3613834E-02_jprb,  1.3842803E+00_jprb,  1.1971031E-01_jprb,  1.4714367E+00_jprb, -2.2001097E+00_jprb, &
-7.3942916E-02_jprb,  1.3590389E+00_jprb,  5.4555789E-01_jprb, -2.0707342E-01_jprb,  1.3061223E+00_jprb, -3.6657591E-01_jprb, &
 2.1585661E-01_jprb, -9.3341347E-01_jprb, -1.0087318E-02_jprb,  1.4606519E-01_jprb,  2.1702319E+00_jprb,  2.1362041E+00_jprb, &
 3.4040008E-01_jprb,  4.4949115E-01_jprb, -3.2672313E+00_jprb, -3.0468373E-01_jprb, -9.4215852E-01_jprb,  1.9318774E-01_jprb, &
 2.1511598E+00_jprb, -2.3460231E+00_jprb,  1.0012196E-01_jprb,  1.8804541E-02_jprb, -1.8651721E+00_jprb,  2.1126784E+00_jprb, &
 1.1603765E+00_jprb, -3.1935876E+00_jprb, -1.2283491E+00_jprb,  8.5411981E-01_jprb,  1.3561147E-01_jprb,  1.5230133E-01_jprb, &
 3.0981923E-01_jprb, -1.4451134E-02_jprb, -1.6729691E+00_jprb, -3.6270868E+00_jprb, -5.5229419E-02_jprb,  2.0263410E+00_jprb, &
 1.3523907E-02_jprb, -1.7415885E-01_jprb,  9.8232003E-01_jprb, -2.6047559E+00_jprb, -6.0255880E+00_jprb,  2.6311645E+00_jprb, &
-1.0178589E+00_jprb,  2.0304214E-01_jprb,  4.1192353E-01_jprb, -8.4691517E-01_jprb, -1.6497784E+00_jprb,  1.0637069E+00_jprb, &
 6.7378525E-02_jprb, -1.5312484E-01_jprb, -1.2114742E-01_jprb, -2.8457813E-01_jprb, -1.5719709E+00_jprb,  2.0276191E+00_jprb, &
-6.2713036E-01_jprb, -1.7351759E+00_jprb, -1.5608109E+00_jprb,  1.6106646E-01_jprb,  1.5816429E-01_jprb,  1.9443840E-01_jprb, &
 2.5119572E-01_jprb, -1.1257213E-01_jprb, -1.2222921E-02_jprb, -1.8812043E-01_jprb, -8.6707870E-02_jprb,  1.2271510E+00_jprb, &
 1.8666324E-01_jprb,  2.3744931E-01_jprb, -6.9319324E-01_jprb,  6.9669835E-01_jprb,  7.7222864E-01_jprb, -2.4162777E-02_jprb, &
-1.5909673E-01_jprb,  1.6569251E-01_jprb, -9.3506730E-01_jprb,  4.9677121E-01_jprb, -1.8218251E-01_jprb, -1.2671604E+00_jprb, &
 1.0232278E-01_jprb, -1.5976154E-02_jprb, -8.1115056E-03_jprb, -4.3076828E-01_jprb,  4.2107766E-01_jprb, -1.0989360E-01_jprb, &
 8.9368367E-01_jprb,  1.2775301E-01_jprb, -1.5592684E-02_jprb,  1.8473056E-01_jprb, -1.7526740E-01_jprb, -4.1306839E-02_jprb, &
-1.5006874E+00_jprb, -4.7304432E-02_jprb,  1.0896354E-01_jprb,  8.0791089E-01_jprb, -6.4487036E-02_jprb, -1.4705799E+00_jprb, &
-7.4696242E-01_jprb,  2.2217487E-01_jprb,  6.5582069E-01_jprb,  5.5189395E-02_jprb,  1.4182230E-01_jprb,  1.1710699E+00_jprb, &
 3.3750001E-01_jprb,  1.1979621E+00_jprb,  9.0584763E-03_jprb,  3.6656275E-03_jprb,  2.6238016E-01_jprb, -2.4194241E-02_jprb, &
 1.3264095E+00_jprb,  3.2595990E-01_jprb, -1.1755738E+00_jprb,  7.5924466E-01_jprb, -6.1733595E-01_jprb,  8.7607882E-03_jprb, &
-2.2433776E+00_jprb, -5.2869355E-02_jprb, -1.1187196E-01_jprb, -1.5418570E+00_jprb,  4.6688317E-03_jprb,  2.9147520E-01_jprb, &
-8.9041057E-01_jprb, -2.3441351E-01_jprb,  1.4632118E-03_jprb, -1.6391210E-01_jprb,  8.5746791E-02_jprb,  3.5223881E-01_jprb, &
-1.5691227E-01_jprb, -1.5714658E-01_jprb,  8.8988684E-02_jprb, -6.0521269E-01_jprb,  4.0557343E-02_jprb,  1.5480291E-01_jprb, &
-1.3446446E+00_jprb, -1.3234933E+00_jprb, -1.8791071E-02_jprb, -5.8021916E-01_jprb,  1.2052445E-01_jprb, -1.4618855E-01_jprb, &
 2.6736133E-01_jprb, -8.4468846E-01_jprb, -1.6087864E-02_jprb, -1.2652253E-01_jprb,  5.4580389E-02_jprb, -2.2615703E-01_jprb, &
-7.2208123E-01_jprb, -7.6635391E-02_jprb, -4.1417422E-01_jprb, -1.3678523E-02_jprb,  6.0300646E-01_jprb, -1.0678104E+00_jprb, &
 2.0435047E-01_jprb,  3.2513737E-01_jprb,  7.5863968E-01_jprb, -1.0834605E+00_jprb, -7.7435914E-02_jprb,  3.8998181E-01_jprb, &
 3.2387861E-01_jprb, -9.1728715E-03_jprb, -5.6322875E-02_jprb, -2.2844879E+00_jprb, -7.9039545E-02_jprb, -1.3625847E-01_jprb, &
 5.1830474E-01_jprb, -1.5745443E-01_jprb,  2.8149116E-01_jprb, -1.0602608E+00_jprb, -6.5258222E-01_jprb,  5.5333785E-01_jprb, &
 5.4168506E-01_jprb, -4.1524734E-01_jprb, -9.0121801E-01_jprb,  7.5596803E-04_jprb,  7.5004318E-03_jprb, -6.6017770E-01_jprb, &
-2.7181873E-01_jprb,  2.1283929E-02_jprb,  1.2839883E+00_jprb,  2.3830709E-01_jprb,  7.4540271E-01_jprb, -2.5006911E+00_jprb, &
 6.9269247E-02_jprb,  5.3190631E-01_jprb, -2.1453568E-01_jprb, -2.3761204E-01_jprb,  1.6789695E+00_jprb, -1.2263047E-01_jprb, &
 4.3418025E-02_jprb, -7.3602433E-02_jprb, -2.2932702E-01_jprb, -1.3059041E-02_jprb, -1.2778069E+00_jprb,  2.6003490E-01_jprb, &
-2.1065197E-01_jprb,  1.1092627E-01_jprb, -2.1431179E-02_jprb, -1.5074173E-01_jprb, -5.3175051E-02_jprb, -5.9542766E-02_jprb, &
 2.1283368E+00_jprb, -8.4029174E-01_jprb,  3.8357756E-02_jprb, -6.2985813E-02_jprb, -1.0010094E+00_jprb,  2.3289235E-01_jprb, &
 5.6829933E-01_jprb, -1.0114262E+00_jprb, -5.4098880E-01_jprb, -3.6090461E-01_jprb, -9.5892004E-02_jprb, -4.1762745E-01_jprb, &
 2.4483744E+00_jprb,  1.6352260E-02_jprb, -6.6315873E-01_jprb,  4.5140575E-01_jprb,  1.7037494E-02_jprb,  1.3609704E+00_jprb, &
 2.0157600E-02_jprb,  2.0300973E+00_jprb,  8.1666539E-02_jprb, -1.5673218E+00_jprb, -6.3772001E-01_jprb,  1.8977429E+00_jprb, &
 2.9041951E-01_jprb,  3.0282848E-02_jprb, -2.7639940E-01_jprb, -3.0412163E-01_jprb, -1.9116702E+00_jprb, -1.0946505E-01_jprb, &
 8.1379080E-01_jprb,  1.2954400E+00_jprb, -5.1277403E-01_jprb,  1.6352420E+00_jprb, -6.0259795E+00_jprb,  1.5160836E+00_jprb, &
-2.6998833E+00_jprb, -1.0090368E+00_jprb, -2.2343756E+00_jprb,  1.9805387E-01_jprb,  9.4680407E-01_jprb,  5.4962146E-01_jprb, &
 1.3292576E-01_jprb, -2.8020518E-01_jprb,  3.9540986E-02_jprb, -3.9303471E-01_jprb, -1.4784787E-01_jprb,  1.5498600E+00_jprb, &
-2.1238867E-01_jprb,  1.0918320E+00_jprb, -3.2074810E-01_jprb, -2.1278551E-01_jprb,  1.0046815E+00_jprb,  1.0630906E+00_jprb, &
-2.3935158E-01_jprb,  1.4384003E-01_jprb, -1.9503640E+00_jprb,  2.9533115E-01_jprb, -8.0276594E-01_jprb,  4.1243078E-01_jprb, &
-4.1277520E-01_jprb, -6.4123919E-01_jprb, -1.1549054E-01_jprb, -2.5558993E-01_jprb, -2.1400587E-01_jprb, -9.5714142E-02_jprb, &
-1.1089254E+00_jprb,  9.6803009E-02_jprb,  9.5059759E-02_jprb,  6.9177742E-02_jprb,  3.2047889E-01_jprb,  1.3053623E+00_jprb, &
 4.9401038E-01_jprb,  1.0453220E-02_jprb, -8.0681393E-02_jprb,  7.3637276E-01_jprb, -5.9703740E-02_jprb, -1.9748803E+00_jprb, &
-2.0447377E-01_jprb,  2.2181216E-01_jprb,  7.5623860E-01_jprb, -2.2389677E-02_jprb,  7.7687976E-01_jprb,  1.9942678E+00_jprb, &
-3.2443634E-01_jprb,  2.0874809E+00_jprb, -1.1621371E-03_jprb,  5.9854027E-02_jprb, -1.1381179E+00_jprb, -4.6230907E-01_jprb, &
 1.4560423E+00_jprb, -3.2265692E-01_jprb, -1.4408204E+00_jprb,  7.6655645E-02_jprb, -3.7226905E-01_jprb, -3.6861789E-02_jprb, &
-2.7430744E+00_jprb,  9.3081457E-02_jprb, -4.7120548E-01_jprb, -5.3013955E+00_jprb, -3.8581400E-02_jprb,  2.3796107E+00_jprb, &
-2.6125626E+00_jprb, -1.0238252E+00_jprb,  1.0613432E-01_jprb,  8.3376466E-01_jprb, -9.3122504E-01_jprb,  1.2015333E+00_jprb, &
 8.1967618E-01_jprb, -6.1322148E-01_jprb, -7.9668075E-01_jprb, -1.9973995E+00_jprb,  3.0941705E-01_jprb,  1.9073165E+00_jprb, &
 1.1602423E+00_jprb, -3.7518506E-01_jprb,  8.6725398E-02_jprb,  4.3669736E-02_jprb, -1.5387107E+00_jprb, -2.2300659E-01_jprb, &
 1.0442598E-01_jprb, -8.8247695E-01_jprb,  8.1013209E-01_jprb, -6.0246773E-02_jprb,  5.1883485E-01_jprb, -6.8442591E-01_jprb, &
 1.8562667E-01_jprb, -9.0610109E-02_jprb, -1.3259335E+00_jprb, -2.6274360E-02_jprb,  7.8372559E-01_jprb, -1.3027352E+00_jprb, &
-8.1302623E-01_jprb,  4.7444156E-01_jprb,  2.9372566E+00_jprb, -1.2160146E+00_jprb,  1.5889218E+00_jprb,  1.3080084E+00_jprb, &
 7.9115830E-01_jprb,  1.4506684E-02_jprb,  6.4966780E-01_jprb, -1.8300539E+00_jprb,  8.5563242E-01_jprb,  6.1180666E-01_jprb, &
 2.2580805E+00_jprb, -3.9691690E-03_jprb,  7.3202982E-01_jprb,  3.9394410E-01_jprb, -7.0534222E-01_jprb,  9.7230607E-02_jprb, &
-3.2992015E+00_jprb,  1.1409947E+00_jprb, -6.2501497E-01_jprb, -2.3363790E-04_jprb, -2.1324168E-02_jprb, -5.6555976E-01_jprb, &
-2.3404408E-01_jprb,  3.3223901E-03_jprb, -6.9187025E-01_jprb,  7.0187434E-01_jprb,  3.9740902E-01_jprb, -5.3295603E-01_jprb, &
 2.1322500E-02_jprb, -1.9755864E-01_jprb, -1.5363311E-01_jprb, -8.6530396E-01_jprb,  3.1915747E+00_jprb, -3.4895361E-02_jprb, &
 1.0982985E-01_jprb,  2.6136235E-01_jprb, -1.7666795E+00_jprb,  6.4129366E-03_jprb, -5.6415723E+00_jprb,  2.0304666E-01_jprb, &
 1.4966792E+00_jprb, -2.8261403E-01_jprb, -8.0268150E-01_jprb,  2.5069764E-01_jprb,  1.1971231E-01_jprb, -1.4010004E-01_jprb, &
 4.4190540E-01_jprb, -2.1707099E-01_jprb, -3.2290309E-03_jprb,  1.6545951E-02_jprb,  1.1839735E+00_jprb, -9.0527418E-01_jprb, &
 1.8511228E+00_jprb,  4.4234130E-01_jprb, -4.1240619E-01_jprb, -4.5363886E-01_jprb, -1.8507337E-01_jprb, -4.3721233E-01_jprb, &
-1.5303662E+00_jprb,  5.6736665E-03_jprb,  1.0266359E+00_jprb,  1.4870760E+00_jprb,  3.7965559E-02_jprb,  1.3018426E+00_jprb, &
 6.4059970E-02_jprb, -2.3869117E+00_jprb, -5.7513609E-02_jprb, -1.0374138E+00_jprb,  1.5212467E+00_jprb,  1.7959294E+00_jprb, &
 1.8406706E-01_jprb,  9.6179974E-02_jprb,  2.9923748E-01_jprb,  3.2559449E-02_jprb,  1.4086012E+00_jprb, -2.4369824E-01_jprb, &
 2.3827844E-02_jprb,  2.6987321E-01_jprb,  2.9828749E-01_jprb, -1.4581791E+00_jprb,  1.1445773E+00_jprb,  4.8001090E-02_jprb, &
 4.8728036E-01_jprb, -1.0988055E-01_jprb, -1.4496207E-02_jprb,  1.2328160E-02_jprb,  1.6071725E-01_jprb, -6.5511239E-02_jprb, &
 2.1521590E-01_jprb, -8.7500754E-02_jprb,  3.9437257E-02_jprb,  1.6053971E-02_jprb,  3.1599478E-02_jprb, -1.0554204E+00_jprb, &
 4.0005799E-02_jprb, -2.6161481E-01_jprb,  7.1961776E-02_jprb, -1.9811868E-01_jprb, -9.5265536E-02_jprb, -7.8659622E-02_jprb, &
 6.9575801E-02_jprb, -1.5004516E-01_jprb,  3.3551043E-02_jprb, -2.2154992E-02_jprb,  5.9450836E-02_jprb,  4.8348751E-01_jprb, &
-6.4878823E-02_jprb,  4.3494442E-01_jprb, -1.1185000E-02_jprb, -1.0170250E-01_jprb, -1.0438476E-01_jprb, -2.3664579E-02_jprb, &
-1.8983858E-01_jprb, -9.6928878E-03_jprb, -7.2938964E-03_jprb,  1.3886746E-01_jprb,  1.4591753E-01_jprb, -4.0708087E+00_jprb, &
-4.9058969E-01_jprb, -7.6822978E-02_jprb, -3.1951276E-02_jprb,  1.2566221E-01_jprb,  1.2139909E-02_jprb,  4.2177220E-02_jprb, &
 2.2495978E+00_jprb, -5.8604879E-02_jprb, -5.5026848E-02_jprb, -1.0021463E-01_jprb, -3.4555191E-01_jprb, -1.7470769E-01_jprb, &
-3.9932435E-02_jprb,  1.8143048E-01_jprb, -3.9578721E-03_jprb,  3.8603038E-02_jprb,  8.4419513E-02_jprb, -4.9810005E-02_jprb, &
 1.0669516E+00_jprb,  1.9285162E-01_jprb, -1.7415507E-02_jprb, -2.7957004E-01_jprb, -1.3836452E-02_jprb, -4.4457030E-02_jprb, &
-3.7942613E-02_jprb,  4.8273428E-02_jprb, -2.6338978E-02_jprb,  1.0877460E+00_jprb, -1.0411855E-02_jprb, -6.9195356E-02_jprb, &
 3.0296699E-01_jprb,  2.4880948E-01_jprb,  7.5701726E-02_jprb, -8.3652985E-02_jprb, -5.5701692E-01_jprb, -1.2951810E-01_jprb, &
-3.1654383E-02_jprb,  2.9917855E-01_jprb, -6.4261028E-02_jprb,  4.3391894E-01_jprb, -1.0920119E-02_jprb,  5.7497164E-01_jprb, &
-2.7208064E-01_jprb,  2.2926530E+00_jprb,  2.3629086E-02_jprb, -6.6584539E-02_jprb,  5.8471379E-02_jprb, -6.8394923E-02_jprb, &
 1.7772778E-03_jprb,  1.0587923E-01_jprb, -3.5519648E-01_jprb, -4.9366386E-01_jprb,  7.3305892E-03_jprb, -2.3849785E-01_jprb, &
-1.2073742E-01_jprb,  5.5600680E-02_jprb,  5.4123179E-02_jprb, -8.4694023E-03_jprb,  1.1505568E-01_jprb, -2.0119830E-02_jprb, &
 5.4296163E-02_jprb, -2.3934309E-01_jprb, -3.3217530E-02_jprb,  1.1621584E-01_jprb,  1.3872963E+00_jprb, -1.6740549E-02_jprb, &
-1.0027480E-01_jprb, -1.0922583E-02_jprb, -1.6454686E-01_jprb, -2.2103687E+00_jprb,  5.0524205E-01_jprb,  3.4619141E-01_jprb, &
-1.9724205E-01_jprb,  2.1208096E-02_jprb,  1.7671213E-02_jprb,  9.5471442E-01_jprb,  8.1040790E-02_jprb, -3.2413877E-01_jprb, &
 1.4339394E+00_jprb, -3.9999602E-01_jprb, -1.2011187E+00_jprb,  4.6543561E-03_jprb, -1.7319917E-02_jprb,  5.1753090E-01_jprb, &
-3.8790816E-01_jprb,  8.9212204E-03_jprb,  1.5651793E-01_jprb, -2.1418356E-02_jprb, -3.1773439E-02_jprb, -9.6587007E-01_jprb, &
 1.2800219E-02_jprb,  6.8967785E-01_jprb,  4.9629784E-02_jprb,  6.5557423E-02_jprb,  4.8581565E-01_jprb, -2.0741530E-03_jprb, &
-4.1746564E-02_jprb, -2.2575711E-02_jprb, -3.3067519E-01_jprb, -2.2599649E-02_jprb,  3.1588667E-01_jprb, -1.6985556E-02_jprb, &
 1.6628496E+00_jprb,  6.0465774E-02_jprb,  2.4014808E+00_jprb,  4.0469393E-02_jprb,  1.3584988E-01_jprb, -3.8456627E-02_jprb, &
 3.7167867E+00_jprb,  2.7505514E-01_jprb, -4.1011010E-02_jprb,  2.0888218E-02_jprb, -3.2710588E+00_jprb,  1.4413220E-01_jprb, &
 5.5513182E-01_jprb, -4.7993179E-02_jprb, -5.3328000E-02_jprb,  2.3697276E-01_jprb,  6.2064911E-02_jprb, -1.4667561E-01_jprb, &
 1.4233136E+01_jprb,  8.3083487E-03_jprb,  3.3526710E-01_jprb,  1.6534487E+00_jprb, -1.0197704E-02_jprb, -8.5443607E-02_jprb, &
 3.7784463E-02_jprb,  1.0385661E+00_jprb, -1.3692545E-01_jprb, -6.6563395E-02_jprb,  2.4313127E+00_jprb, -1.1467764E-01_jprb /

!== ti: structure having the settings of the input mapminmax transformation

! 7 elements in ti_gain
DATA aniso_net%ti_gain / 2.8591851e-03_jprb, 2.2471910e-02_jprb, 4.0000000e-02_jprb, 6.2500000e-02_jprb, &
                         5.0632911e-02_jprb, 2.3269984e+00_jprb, 2.2429120e+00_jprb /

! 7 elements in ti_xoffset 
DATA aniso_net%ti_xoffset / 5.0000000e-01_jprb, 0.0000000e+00_jprb, 0.0000000e+00_jprb, 2.7115000e+02_jprb, &
                            0.0000000e+00_jprb, 1.4029873e-01_jprb, 3.5592941e-03_jprb /

! 1 element in ti_ymin
DATA aniso_net%ti_ymin / -1.0000000e+00_jprb /

!== to: structure having the settings of the output mapminmax transformation

! 8 elements in to_gain
DATA aniso_net%to_gain / 3.0580307e+01_jprb, 4.1162395e+01_jprb, 5.3868692e+01_jprb, 1.1524030e+03_jprb, &
                         3.7019172e+01_jprb, 4.0606951e+01_jprb, 2.2551726e+01_jprb, 5.5613649e+01_jprb /

! 8 elements in to_xoffset 
DATA aniso_net%to_xoffset / -3.3371218e-02_jprb, -5.7938410e-03_jprb, -3.6754564e-02_jprb, -1.7318159e-03_jprb, &
                            -2.8659080e-02_jprb, -2.5723239e-02_jprb, -3.6839388e-02_jprb, -1.0555830e-04_jprb /

! 1 element in to_ymin
DATA aniso_net%to_ymin / -1.0000000e+00_jprb /

END MODULE rttov_surfem_ocean_coef_mod
