! Description:
!> @file
!!   Contains data for the FASTEM-1,2,3 MW sea surface emissivity models
!
!> @brief
!!   Contains data for the FASTEM-1,2,3 MW sea surface emissivity models
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
MODULE mod_rttov_fastem3_coef

  USE parkind1, ONLY : jprb, jpim

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fastem3_coef, freqfixed

  REAL(KIND=jprb), PARAMETER :: freqfixed(4) = (/7.0_jprb, 10.0_jprb, 19.0_jprb, 37.0_jprb/)

  REAL(KIND=jprb), PARAMETER :: fastem3_coef_1(125) = (/ &
  0.175350E+02_jprb, -0.617670E+00_jprb,  0.894800E-02_jprb,  0.318420E+01_jprb,  0.191890E-01_jprb, &
 -0.108730E-01_jprb,  0.258180E-03_jprb,  0.683960E+02_jprb, -0.406430E+00_jprb,  0.228320E-01_jprb, &
 -0.530610E-03_jprb,  0.476290E+01_jprb,  0.154100E+00_jprb, -0.337170E-01_jprb,  0.844280E-03_jprb, &
  0.782870E+02_jprb, -0.434630E-02_jprb,  0.531250E+01_jprb, -0.114770E-01_jprb,  0.314160E+01_jprb, &
 -0.100000E+01_jprb,  0.195000E-04_jprb,  0.255000E+01_jprb, -0.182390E+01_jprb, -0.434790E-02_jprb, &
  0.646320E-04_jprb,  0.278640E+01_jprb,  0.878460E-02_jprb, -0.102670E-03_jprb, -0.101890E+01_jprb, &
 -0.426820E-02_jprb,  0.396520E-04_jprb,  0.730720E-01_jprb,  0.261790E-02_jprb, -0.950500E-05_jprb, &
  0.295330E-03_jprb,  0.443690E-05_jprb, -0.140160E-07_jprb, -0.717940E-01_jprb, -0.267870E-02_jprb, &
  0.949560E-05_jprb, -0.334690E+00_jprb,  0.951660E-02_jprb,  0.964400E-05_jprb,  0.470780E+00_jprb, &
 -0.148970E-01_jprb, -0.987460E-05_jprb, -0.142750E+00_jprb,  0.565380E-02_jprb,  0.118850E-05_jprb, &
 -0.137840E+00_jprb, -0.216950E-02_jprb,  0.793130E-05_jprb,  0.237840E-04_jprb,  0.869500E-06_jprb, &
  0.282490E-08_jprb,  0.138790E+00_jprb,  0.209470E-02_jprb, -0.797900E-05_jprb, -0.637180E+01_jprb, &
  0.253920E-01_jprb,  0.357570E-04_jprb,  0.942930E+01_jprb, -0.332840E-01_jprb, -0.647720E-04_jprb, &
 -0.329280E+01_jprb,  0.965450E-02_jprb,  0.281590E-04_jprb,  0.252680E+00_jprb,  0.343870E-02_jprb, &
 -0.156360E-04_jprb, -0.156670E-03_jprb,  0.139490E-04_jprb, -0.407630E-07_jprb, -0.141320E+00_jprb, &
 -0.356560E-02_jprb,  0.142870E-04_jprb, -0.240700E+01_jprb, -0.563890E-01_jprb,  0.325230E-03_jprb, &
  0.296010E+01_jprb,  0.704680E-01_jprb, -0.426440E-03_jprb, -0.751250E+00_jprb, -0.191930E-01_jprb, &
  0.125940E-03_jprb, -0.288250E+00_jprb, -0.102650E-02_jprb,  0.226700E-05_jprb, -0.119070E-02_jprb, &
 -0.263170E-04_jprb,  0.114600E-06_jprb,  0.406300E+00_jprb,  0.200030E-02_jprb, -0.781640E-05_jprb, &
 -0.675700E-01_jprb,  0.214600E+00_jprb, -0.363000E-02_jprb,  0.636730E+01_jprb,  0.900610E+00_jprb, &
 -0.524880E+00_jprb, -0.370920E+01_jprb, -0.143310E+01_jprb,  0.397450E+00_jprb,  0.823100E-01_jprb, &
 -0.255980E+00_jprb,  0.552000E-02_jprb,  0.208000E+01_jprb,  0.244920E+01_jprb, -0.456420E+00_jprb, &
 -0.224900E-01_jprb,  0.616900E-01_jprb, -0.344000E-02_jprb, -0.507570E+01_jprb, -0.360670E+01_jprb, &
  0.118750E+01_jprb,  0.124950E+00_jprb,  0.121270E+00_jprb,  0.714000E-02_jprb,  0.736620E+01_jprb, &
 -0.114060E+00_jprb, -0.272910E+00_jprb, -0.504350E+01_jprb, -0.336450E+00_jprb,  0.161260E+00_jprb /)

  REAL(KIND=jprb), PARAMETER :: fastem3_coef_2(125) = (/ &
 -0.154290E+00_jprb, -0.141070E+00_jprb, -0.809000E-02_jprb,  0.395290E+01_jprb,  0.958580E+00_jprb, &
 -0.159080E+00_jprb,  0.368500E-01_jprb,  0.307100E-01_jprb,  0.810000E-03_jprb, -0.619960E+01_jprb, &
 -0.172580E+01_jprb,  0.641360E+00_jprb,  0.100000E+01_jprb,  0.200000E-01_jprb,  0.300000E+00_jprb, &
 -0.585336E-04_jprb,  0.141135E-03_jprb,  0.341558E-05_jprb,  0.163655E-07_jprb,  0.184676E-03_jprb, &
 -0.956046E-04_jprb,  0.544262E-05_jprb, -0.121126E-06_jprb, -0.453125E-04_jprb,  0.154844E-04_jprb, &
 -0.812972E-06_jprb,  0.160754E-07_jprb, -0.590621E-05_jprb,  0.721022E-04_jprb,  0.331280E-05_jprb, &
 -0.116781E-08_jprb,  0.790312E-03_jprb, -0.345218E-03_jprb,  0.146029E-04_jprb, -0.312496E-06_jprb, &
  0.312823E-04_jprb, -0.138377E-04_jprb,  0.225909E-07_jprb,  0.229783E-08_jprb,  0.550691E-04_jprb, &
 -0.106330E-03_jprb,  0.453266E-06_jprb,  0.909021E-08_jprb,  0.694857E-03_jprb, -0.286702E-03_jprb, &
  0.944863E-05_jprb, -0.217880E-06_jprb,  0.612423E-04_jprb, -0.315418E-04_jprb,  0.107982E-05_jprb, &
 -0.250954E-07_jprb,  0.212483E-04_jprb, -0.598084E-05_jprb, -0.609132E-06_jprb,  0.159695E-07_jprb, &
 -0.613516E-03_jprb,  0.263937E-03_jprb, -0.762252E-05_jprb,  0.180749E-06_jprb,  0.339556E-05_jprb, &
 -0.235102E-05_jprb,  0.489815E-06_jprb, -0.139383E-07_jprb, -0.408839E-04_jprb,  0.115903E-03_jprb, &
  0.133087E-04_jprb, -0.232691E-06_jprb,  0.212243E-03_jprb, -0.106434E-03_jprb,  0.465887E-05_jprb, &
 -0.461192E-07_jprb, -0.621221E-04_jprb,  0.177511E-04_jprb,  0.799048E-07_jprb, -0.576266E-08_jprb, &
  0.381961E-04_jprb,  0.585374E-04_jprb,  0.756621E-05_jprb, -0.110985E-06_jprb,  0.129621E-02_jprb, &
 -0.554118E-03_jprb,  0.173120E-04_jprb, -0.276256E-06_jprb,  0.515273E-04_jprb, -0.207117E-04_jprb, &
  0.487256E-07_jprb, -0.235559E-08_jprb,  0.152246E-03_jprb, -0.159825E-03_jprb, -0.500807E-05_jprb, &
  0.126266E-06_jprb,  0.102881E-02_jprb, -0.414324E-03_jprb,  0.923915E-05_jprb, -0.150247E-06_jprb, &
  0.888053E-04_jprb, -0.392334E-04_jprb, -0.688354E-06_jprb,  0.293177E-07_jprb,  0.504310E-04_jprb, &
 -0.191818E-04_jprb,  0.490998E-06_jprb, -0.158696E-07_jprb, -0.615485E-03_jprb,  0.257073E-03_jprb, &
 -0.467360E-05_jprb,  0.676351E-07_jprb,  0.247840E-05_jprb, -0.154153E-05_jprb,  0.333460E-06_jprb, &
 -0.784914E-08_jprb, -0.621877E-04_jprb,  0.124143E-03_jprb,  0.170023E-04_jprb, -0.368643E-06_jprb, &
  0.101425E-03_jprb, -0.630114E-04_jprb,  0.435736E-05_jprb, -0.101644E-06_jprb, -0.796174E-04_jprb, &
  0.265038E-04_jprb,  0.537454E-07_jprb, -0.145468E-07_jprb,  0.442053E-04_jprb,  0.459572E-04_jprb /)

  REAL(KIND=jprb), PARAMETER :: fastem3_coef_3(125) = (/ &
  0.971810E-05_jprb, -0.170817E-06_jprb,  0.133357E-02_jprb, -0.557281E-03_jprb,  0.142888E-04_jprb, &
 -0.171095E-06_jprb,  0.531728E-04_jprb, -0.217787E-04_jprb,  0.245581E-06_jprb, -0.101500E-07_jprb, &
  0.115092E-03_jprb, -0.140989E-03_jprb, -0.103311E-04_jprb,  0.255829E-06_jprb,  0.109355E-02_jprb, &
 -0.439655E-03_jprb,  0.847483E-05_jprb, -0.100246E-06_jprb,  0.148653E-03_jprb, -0.629077E-04_jprb, &
  0.114331E-05_jprb, -0.215387E-07_jprb,  0.132775E-04_jprb, -0.418720E-05_jprb, -0.105548E-05_jprb, &
  0.289720E-07_jprb, -0.572372E-03_jprb,  0.234306E-03_jprb, -0.264171E-05_jprb,  0.348850E-08_jprb, &
  0.155316E-04_jprb, -0.823120E-05_jprb,  0.106889E-05_jprb, -0.298319E-07_jprb,  0.863755E-05_jprb, &
  0.595888E-04_jprb,  0.254421E-04_jprb, -0.413468E-06_jprb,  0.227688E-03_jprb, -0.113986E-03_jprb, &
  0.725093E-05_jprb, -0.127069E-06_jprb, -0.337521E-04_jprb,  0.437196E-05_jprb,  0.256526E-05_jprb, &
 -0.749534E-07_jprb,  0.714562E-04_jprb,  0.201237E-04_jprb,  0.138150E-04_jprb, -0.188276E-06_jprb, &
  0.130476E-02_jprb, -0.520253E-03_jprb,  0.463495E-05_jprb,  0.320702E-07_jprb,  0.391178E-04_jprb, &
 -0.155648E-04_jprb, -0.461218E-06_jprb, -0.361295E-08_jprb,  0.111352E-03_jprb, -0.122809E-03_jprb, &
 -0.186779E-04_jprb,  0.325278E-06_jprb,  0.100263E-02_jprb, -0.378765E-03_jprb,  0.246420E-06_jprb, &
  0.317558E-07_jprb,  0.127648E-03_jprb, -0.498883E-04_jprb, -0.867037E-06_jprb,  0.147253E-07_jprb, &
  0.107976E-04_jprb, -0.382161E-05_jprb, -0.949457E-06_jprb,  0.158475E-07_jprb, -0.478420E-03_jprb, &
  0.182279E-03_jprb,  0.881766E-06_jprb, -0.359735E-07_jprb,  0.186481E-04_jprb, -0.547143E-05_jprb, &
  0.498428E-06_jprb, -0.826455E-08_jprb, -0.585336E-04_jprb,  0.141135E-03_jprb,  0.341558E-05_jprb, &
  0.163655E-07_jprb,  0.184676E-03_jprb, -0.956046E-04_jprb,  0.544262E-05_jprb, -0.121126E-06_jprb, &
 -0.453125E-04_jprb,  0.154844E-04_jprb, -0.812972E-06_jprb,  0.160754E-07_jprb, -0.590621E-05_jprb, &
  0.721022E-04_jprb,  0.331280E-05_jprb, -0.116781E-08_jprb,  0.790312E-03_jprb, -0.345218E-03_jprb, &
  0.146029E-04_jprb, -0.312496E-06_jprb,  0.312823E-04_jprb, -0.138377E-04_jprb,  0.225909E-07_jprb, &
  0.229783E-08_jprb,  0.550691E-04_jprb, -0.106330E-03_jprb,  0.453266E-06_jprb,  0.909021E-08_jprb, &
  0.694857E-03_jprb, -0.286702E-03_jprb,  0.944863E-05_jprb, -0.217880E-06_jprb,  0.612423E-04_jprb, &
 -0.315418E-04_jprb,  0.107982E-05_jprb, -0.250954E-07_jprb,  0.212483E-04_jprb, -0.598084E-05_jprb, &
 -0.609132E-06_jprb,  0.159695E-07_jprb, -0.613516E-03_jprb,  0.263937E-03_jprb, -0.762252E-05_jprb /)

  REAL(KIND=jprb), PARAMETER :: fastem3_coef_4(125) = (/ &
  0.180749E-06_jprb,  0.339556E-05_jprb, -0.235102E-05_jprb,  0.489815E-06_jprb, -0.139383E-07_jprb, &
 -0.408839E-04_jprb,  0.115903E-03_jprb,  0.133087E-04_jprb, -0.232691E-06_jprb,  0.212243E-03_jprb, &
 -0.106434E-03_jprb,  0.465887E-05_jprb, -0.461192E-07_jprb, -0.621221E-04_jprb,  0.177511E-04_jprb, &
  0.799048E-07_jprb, -0.576266E-08_jprb,  0.381961E-04_jprb,  0.585374E-04_jprb,  0.756621E-05_jprb, &
 -0.110985E-06_jprb,  0.129621E-02_jprb, -0.554118E-03_jprb,  0.173120E-04_jprb, -0.276256E-06_jprb, &
  0.515273E-04_jprb, -0.207117E-04_jprb,  0.487256E-07_jprb, -0.235559E-08_jprb,  0.152246E-03_jprb, &
 -0.159825E-03_jprb, -0.500807E-05_jprb,  0.126266E-06_jprb,  0.102881E-02_jprb, -0.414324E-03_jprb, &
  0.923915E-05_jprb, -0.150247E-06_jprb,  0.888053E-04_jprb, -0.392334E-04_jprb, -0.688354E-06_jprb, &
  0.293177E-07_jprb,  0.504310E-04_jprb, -0.191818E-04_jprb,  0.490998E-06_jprb, -0.158696E-07_jprb, &
 -0.615485E-03_jprb,  0.257073E-03_jprb, -0.467360E-05_jprb,  0.676351E-07_jprb,  0.247840E-05_jprb, &
 -0.154153E-05_jprb,  0.333460E-06_jprb, -0.784914E-08_jprb, -0.621877E-04_jprb,  0.124143E-03_jprb, &
  0.170023E-04_jprb, -0.368643E-06_jprb,  0.101425E-03_jprb, -0.630114E-04_jprb,  0.435736E-05_jprb, &
 -0.101644E-06_jprb, -0.796174E-04_jprb,  0.265038E-04_jprb,  0.537454E-07_jprb, -0.145468E-07_jprb, &
  0.442053E-04_jprb,  0.459572E-04_jprb,  0.971810E-05_jprb, -0.170817E-06_jprb,  0.133357E-02_jprb, &
 -0.557281E-03_jprb,  0.142888E-04_jprb, -0.171095E-06_jprb,  0.531728E-04_jprb, -0.217787E-04_jprb, &
  0.245581E-06_jprb, -0.101500E-07_jprb,  0.115092E-03_jprb, -0.140989E-03_jprb, -0.103311E-04_jprb, &
  0.255829E-06_jprb,  0.109355E-02_jprb, -0.439655E-03_jprb,  0.847483E-05_jprb, -0.100246E-06_jprb, &
  0.148653E-03_jprb, -0.629077E-04_jprb,  0.114331E-05_jprb, -0.215387E-07_jprb,  0.132775E-04_jprb, &
 -0.418720E-05_jprb, -0.105548E-05_jprb,  0.289720E-07_jprb, -0.572372E-03_jprb,  0.234306E-03_jprb, &
 -0.264171E-05_jprb,  0.348850E-08_jprb,  0.155316E-04_jprb, -0.823120E-05_jprb,  0.106889E-05_jprb, &
 -0.298319E-07_jprb,  0.863755E-05_jprb,  0.595888E-04_jprb,  0.254421E-04_jprb, -0.413468E-06_jprb, &
  0.227688E-03_jprb, -0.113986E-03_jprb,  0.725093E-05_jprb, -0.127069E-06_jprb, -0.337521E-04_jprb, &
  0.437196E-05_jprb,  0.256526E-05_jprb, -0.749534E-07_jprb,  0.714562E-04_jprb,  0.201237E-04_jprb, &
  0.138150E-04_jprb, -0.188276E-06_jprb,  0.130476E-02_jprb, -0.520253E-03_jprb,  0.463495E-05_jprb, &
  0.320702E-07_jprb,  0.391178E-04_jprb, -0.155648E-04_jprb, -0.461218E-06_jprb, -0.361295E-08_jprb /)

  REAL(KIND=jprb), PARAMETER :: fastem3_coef_5(24) = (/ &
  0.111352E-03_jprb, -0.122809E-03_jprb, -0.186779E-04_jprb,  0.325278E-06_jprb,  0.100263E-02_jprb, &
 -0.378765E-03_jprb,  0.246420E-06_jprb,  0.317558E-07_jprb,  0.127648E-03_jprb, -0.498883E-04_jprb, &
 -0.867037E-06_jprb,  0.147253E-07_jprb,  0.107976E-04_jprb, -0.382161E-05_jprb, -0.949457E-06_jprb, &
  0.158475E-07_jprb, -0.478420E-03_jprb,  0.182279E-03_jprb,  0.881766E-06_jprb, -0.359735E-07_jprb, &
  0.186481E-04_jprb, -0.547143E-05_jprb,  0.498428E-06_jprb, -0.826455E-08_jprb /)

  INTEGER(KIND=jpim) :: i

  !> FASTEM-1,2,3 coefficients
  REAL(KIND=jprb), TARGET :: fastem3_coef(524) = (/ &
                (fastem3_coef_1 (i), i = 1, 125),   &
                (fastem3_coef_2 (i), i = 1, 125),   &
                (fastem3_coef_3 (i), i = 1, 125),   &
                (fastem3_coef_4 (i), i = 1, 125),   &
                (fastem3_coef_5 (i), i = 1,  24) /)

END MODULE mod_rttov_fastem3_coef
