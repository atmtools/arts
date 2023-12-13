!
  SUBROUTINE rttov_fastem5( fastem_version, &  ! Input
                              Frequency   , &  ! Input
                              Zenith_Angle, &  ! Input
                              Temperature , &  ! Input
                              Salinity    , &  ! Input
                              Wind_Speed  , &  ! Input
                              Emissivity  , &  ! Output
                              Reflectivity, &  ! Output
                              Transmittance,&  ! Input, may not be used
                              Rel_Azimuth  )   ! Input, may not be used
! Description:
!> @file
!!   Compute FASTEM-4,5,6 emissivity and reflectance for a single channel
!
!> @brief
!!   Compute FASTEM-4,5,6 emissivity and reflectance for a single channel
!!
!! @details
!!   References for FASTEM are given in the user guide.
!!
!! @param[in]     fastem_version        FASTEM version to compute (4, 5 or 6)
!! @param[in]     frequency             channel frequency (GHz)
!! @param[in]     zenith_angle          profile zenith angle (degrees)
!! @param[in]     temperature           profile skin temperature (K)
!! @param[in]     salinity              profile salinity (practical salinity units)
!! @param[in]     wind_speed            profile wind speed (m/s)
!! @param[out]    emissivity            calculated emissivity (4 Stokes components)
!! @param[out]    reflectivity          calculated reflectivity (4 Stokes components)
!! @param[in]     transmittance         surface-to-space transmittance
!! @param[in]     rel_azimuth           relative azimuth angle
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
!INTF_OFF
    USE mod_rttov_fastem5_coef, ONLY : FresnelVariables_type, PermittivityVariables_type,&
        ZERO, ONE, TWO, PI, DEGREES_TO_RADIANS, transmittance_limit_lower,&
        transmittance_limit_upper, e0_4, e0_5, min_f, max_f, min_wind, max_wind, A_COEF, Lcoef4, Lcoef5,&
        Scoef, t_c4, t_c5, b_coef, FR_COEFF, x, y, coef_mk_azi
    USE mod_rttov_fastem3_coef, ONLY : fastem3_coef, freqfixed
    USE parkind1, ONLY : jprb
!INTF_ON
    USE mod_rttov_fastem5_coef, ONLY: fp
    USE parkind1, ONLY : jpim, jplm
    IMPLICIT NONE
    ! Arguments
    INTEGER(jpim),   INTENT(IN)     :: fastem_version
    REAL(fp),        INTENT(IN)     :: Frequency
    REAL(fp),        INTENT(IN)     :: Zenith_Angle
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Salinity
    REAL(fp),        INTENT(IN)     :: Wind_Speed
    REAL(fp),        INTENT(OUT)    :: Emissivity(4), Reflectivity(4)
    REAL(fp), INTENT(IN)  :: Transmittance
    REAL(fp), INTENT(IN)  :: Rel_Azimuth

!INTF_END

  !local variable

    REAL(fp) :: e0
    REAL(fp) :: cos_z, Foam_Cover
    REAL(fp) :: scor, small_corr, Azimuth_Emi(4),RV_Fresnel,RH_Fresnel
    REAL(fp) :: Ev,Eh,RvL,RhL,RvS,RhS
    REAL(fp) :: zreflmod_v,zreflmod_h,zrough_v,zrough_h
    INTEGER :: i, j, L, m

    REAL(fp) :: Foam_Rv,Foam_Rh

    ! Local variables for the permittivity model
    REAL(fp) :: einf, sigma25
    REAL(fp) :: tau1, tau2, es, e1
    REAL(fp) :: perm_Real, perm_imag
    TYPE(PermittivityVariables_type) :: iVar
    COMPLEX( fp ) :: Permittivity

    ! Local variables for Fresnel reflectance
    COMPLEX(fp) :: zRv ! Vertical
    COMPLEX(fp) :: zRh ! Horizontal
    TYPE(FresnelVariables_type) :: frVar

    ! Local variables for small-scale
    REAL(fp) :: windspeed, freq_S
    ! Local variables for large-scale
    REAL(fp) :: seczen, zc(12)

    ! Local variables for including transmittance
    REAL(fp) :: variance,varm,opdpsfc,zx(9)
    ! Foam reflectance
    REAL(fp) :: Fh, Foam_ref
    ! Local variables for azimuth angle
    REAL(fp) :: ac, sc, fre_c, phi, wind10
    ! Local arrays to hold FASTEM-4/5 coefs
    REAL(fp) :: Lcoef(size(Lcoef5)), t_c(size(t_c5))

    INTEGER  :: ifreq,ipol
    REAL(fp) :: azimuth_component(2,6) ! pol,  freq
    REAL(fp),parameter :: xs11=2
    REAL(fp),parameter :: xs12=2
    REAL(fp),parameter :: xs21=1
    REAL(fp),parameter :: xs22=4
    REAL(fp),parameter :: theta_ref=55.2d0
    REAL(fp) :: theta

    ! FASTEM-3 azimuth calculation
    INTEGER(jpim) :: istokes, ich, i_freq
    REAL(jprb)    :: u19, freq_ghz, dfreq
    REAL(jprb)    :: tbfixed(4,4,3)     ! Surface brightness temperature azimuthal variation terms for 37, 19, 10, 7 GHz
    REAL(jprb)    :: efixed(4,4,3)      ! Emissivity azimuthal variation terms for 7, 10, 19, 37 GHz
    REAL(jprb)    :: einterpolated(4,3) ! Emissivity azimuthal variation terms for interpolated to required frequency
    REAL(jprb)    :: a1e, a2e, a3e      ! coefficients used in azimuthal emissivity model
    REAL(jprb), POINTER :: c(:)! pointer to FASTEM-3 coefs

    REAL(fp),dimension(6) :: A1v,A1h,A2v,A2h
    REAL(fp),dimension(6) :: A1s1,A1s2,A2s1,A2s2
    REAL(fp),dimension(6) :: A2s2_theta0, A1s1_theta,A2s1_theta,A1s2_theta,A2s2_theta
    REAL(fp),dimension(6) :: A1v_theta, A1h_theta,A2v_theta,A2h_theta
    REAL(fp) :: fratio

    IF (fastem_version == 4) THEN
      e0 = e0_4
      Lcoef = Lcoef4
      t_c = t_c4
    ELSE
      e0 = e0_5
      Lcoef = Lcoef5
      t_c = t_c5
    ENDIF
    cos_z = cos( Zenith_Angle*DEGREES_TO_RADIANS )

  ! Permittivity Calculation
  ! ------------------------
    !1.2 calculate permittivity using double-debye formula
    !-----------------------------------------------------
    !Set values for temperature polynomials (convert from kelvin to celsius)
    iVar%t = Temperature - 273.15_fp
    iVar%t_sq = iVar%t * iVar%t     !quadratic
    iVar%t_cu = iVar%t_sq * iVar%t  !cubic
    iVar%S = Salinity
    !-----------------------------------------------------
    !1.2 Pure or fresh water
    !-----------------------------------------------------
    einf = A_COEF(0) + A_COEF(1)*iVar%t
    es   = A_COEF(2) + A_COEF(3)*iVar%t  + A_COEF(4)*iVar%t_sq + A_COEF(5)*iVar%t_cu
    e1   = A_COEF(9) + A_COEF(10)*iVar%t + A_COEF(11)*iVar%t_sq
    tau1 = A_COEF(15) + A_COEF(16)*iVar%t + A_COEF(17)*iVar%t_sq + A_COEF(18)*iVar%t_cu
    tau2 = A_COEF(22) + A_COEF(23)*iVar%t + A_COEF(24)*iVar%t_sq + A_COEF(25)*iVar%t_cu

    iVar%es_k = es
    iVar%e1_k = e1
    iVar%tau1_k = tau1
    iVar%tau2_k = tau2
    perm_imag = ZERO

    IF( iVar%S > ZERO ) THEN
      iVar%delta = 25.0_fp - iVar%t
      iVar%beta  = A_COEF(29) +A_COEF(30)*iVar%delta +A_COEF(31)*iVar%delta**2  &
            + iVar%S*(A_COEF(32) +A_COEF(33)*iVar%delta +A_COEF(34)*iVar%delta**2)
      sigma25 = iVar%S*(A_COEF(35) +A_COEF(36)*iVar%S +A_COEF(37)*iVar%S**2  &
              +A_COEF(38)*iVar%S**3)
      iVar%sigma = sigma25*exp(-iVar%delta*iVar%beta)

      iVar%ces = ONE + iVar%S*(A_COEF(6) + A_COEF(7)*iVar%S + A_COEF(8)*iVar%t )
      iVar%ce1 = ONE + iVar%S*(A_COEF(12) + A_COEF(13)*iVar%S +A_COEF(14)*iVar%t )
      iVar%ctau1 = ONE + iVar%S*(A_COEF(19) +A_COEF(20)*iVar%t + A_COEF(21)*iVar%t_sq)
      iVar%ctau2 = ONE + iVar%S*(A_COEF(26) + A_COEF(27)*iVar%t + A_COEF(28)*iVar%S**2 )
      es = iVar%es_k * iVar%ces
      e1 = iVar%e1_k * iVar%ce1
      tau1 = iVar%tau1_k * iVar%ctau1
      tau2 = iVar%tau2_k * iVar%ctau2
      perm_imag = -iVar%sigma/(TWO*PI*e0*Frequency)
    END IF
    !Define two relaxation frequencies, f1 and f2
    iVar%f1 = Frequency*tau1
    iVar%f2 = Frequency*tau2
    iVar%del1 = es - e1
    iVar%del2 = e1 - einf

    perm_Real = einf + iVar%del1/(ONE + iVar%f1**2) + iVar%del2/(ONE + iVar%f2**2)
    perm_imag = -perm_imag + iVar%del1*iVar%f1/(ONE + iVar%f1**2)  &
              + iVar%del2*iVar%f2/(ONE + iVar%f2**2)
    Permittivity = Cmplx(perm_Real,-perm_imag,fp)

  ! Compute Fresnel reflectance code, adopted from Masahiro Kazumori, JMA
  !
    ! Compute the complex reflectivity components
    frVar%z1 = SQRT(permittivity - ONE + (cos_z*cos_z))
    frVar%z2 = permittivity * cos_z
    zRh = (cos_z  -frVar%z1) / (cos_z  +frVar%z1)
    zRv = (frVar%z2-frVar%z1) / (frVar%z2+frVar%z1)

    ! The square of the vertical abs value
    frVar%rzRv = REAL(zRv,fp)
    frVar%izRv = AIMAG(zRv)
    Rv_Fresnel = frVar%rzRv**2 + frVar%izRv**2

    ! The square of the horizontal abs value
    frVar%rzRh = REAL(zRh,fp)
    frVar%izRh = AIMAG(zRh)
    Rh_Fresnel = frVar%rzRh**2 + frVar%izRh**2


  ! Apply small-scale correction
  ! --------------------------------
  ! Note from Steve English: 'windspeed' is restricted to be between min_wind and
  ! max_wind here. After this section the unrestricted 'wind_speed' is used. This
  ! is done intentionally.
    windspeed = Wind_Speed
    IF( windspeed < min_wind ) windspeed = min_wind
    IF( windspeed > max_wind ) windspeed = max_wind

    freq_S = Frequency
    IF( freq_S < min_f ) freq_S = min_f
    IF( freq_S > max_f ) freq_S = max_f

    scor = Scoef(1) *windspeed*freq_S +Scoef(2) *windspeed*freq_S**2 &
           + Scoef(3) *windspeed**2* freq_S +Scoef(4) *windspeed**2* freq_S**2 &
           + Scoef(5) *windspeed**2 /freq_S +Scoef(6) *windspeed**2 /freq_S**2 &
           + Scoef(7) *windspeed + Scoef(8) *windspeed**2

    small_corr = exp(-scor*cos_z*cos_z )
    RvS = Rv_Fresnel * small_corr
    RhS = Rh_Fresnel * small_corr

  ! Large Scale Correction Calculation
  ! ----------------------------------
    seczen = ONE/cos_z
    ! compute fitting coefficients for a given frequency
    DO j = 1, 12
      zc(j) = Lcoef(j*3-2) + Lcoef(j*3-1)*frequency + Lcoef(j*3)*frequency**2
    END DO

    RvL = zc(1) + zc(2)*seczen + zc(3)*seczen**2 + zc(4)*Wind_Speed &
      + zc(5)*Wind_Speed**2 + zc(6)*Wind_Speed*seczen
    RhL = zc(7) + zc(8)*seczen + zc(9)*seczen**2 + zc(10)*Wind_Speed &
      + zc(11)*Wind_Speed**2 + zc(12)*Wind_Speed*seczen
    
    ! change foam coverage back to FASTEM1,2,3 for FASTEM5
    IF (fastem_version == 4) THEN
      ! Compute foam coverage after Tang, 1974
      Foam_Cover = 7.75E-06_fp * Wind_Speed ** 3.231_fp
    ELSE
      ! Monahan et al., 1986 without surface stability term
      Foam_Cover = 1.95E-05_fp * Wind_Speed ** 2.55_fp
    END IF
    
  ! The foam vertical and horizontal reflectanc codes, adopted from Masahiro Kazumori, JMA
  ! ----------------------------------
    Foam_Rv = FR_COEFF(1)
    Fh = ONE + Zenith_Angle*(FR_COEFF(2) +  Zenith_Angle*(FR_COEFF(3)  &
       + Zenith_Angle*FR_COEFF(4)))
    Foam_Rh = ONE + FR_COEFF(5)*Fh

    ! Added frequency dependence derived from Stogry model
    Foam_ref = 0.4_fp * exp(-0.05_fp*Frequency )
    Foam_Rv = Foam_Rv * Foam_ref
    Foam_Rh = Foam_Rh * Foam_ref

    Ev = (ONE-Foam_Cover)*(ONE - RvS + RvL) + Foam_Cover*(ONE-Foam_Rv)
    Eh = (ONE-Foam_Cover)*(ONE - RhS + RhL) + Foam_Cover*(ONE-Foam_Rh)

    Emissivity(1) = Ev
    Emissivity(2) = Eh

    zreflmod_v = ONE
    zreflmod_h = ONE

  ! correction for anisotropic downward radiation, adopted from the FASTEM3
  ! ----------------------------------
    IF( Transmittance > transmittance_limit_lower .and. Transmittance < transmittance_limit_upper) THEN
        !Using the Cox and Munk model to compute slope variance
        variance = 0.00512_fp * Wind_Speed + 0.0030_fp
        varm     = variance * t_c(43)
        variance = varm * ( t_c(44) * Frequency + t_c(45) )
        IF ( variance >= varm ) variance = varm
        IF ( variance <= ZERO  ) variance = ZERO
        !Compute surface to space optical depth
        opdpsfc = -log(Transmittance ) * cos_z

        !Define nine predictors for the effective angle calculation
        zx(1) = ONE
        zx(2) = variance
        zx(4) = ONE / cos_z
        zx(3) = zx(2) * zx(4)
        zx(5) = zx(3) * zx(3)
        zx(6) = zx(4) * zx(4)
        zx(7) = zx(2) * zx(2)
        zx(8) = log(opdpsfc)
        zx(9) = zx(8) * zx(8)

        zrough_v = ONE
        zrough_h = ONE
        DO i = 1, 7
           j = i-1
          !Switched h to v Deblonde SSMIS june 7, 2001
          zrough_h = zrough_h + zx(i) *(t_c(1+j*3) + zx(8)*t_c(2+j*3) + zx(9)*t_c(3+j*3) )
          zrough_v = zrough_v + zx(i) *(t_c(22+j*3)+ zx(8)*t_c(23+j*3)+ zx(9)*t_c(24+j*3))
        END DO
        zreflmod_v = (ONE-Transmittance ** zrough_v)/(ONE-Transmittance )
        zreflmod_h = (ONE-Transmittance ** zrough_h)/(ONE-Transmittance )

    END IF

  ! azimuthal component
  ! --------------------------------
    Azimuth_Emi = ZERO

    IF( abs(Rel_Azimuth) <= 360.0_fp ) THEN
      if(fastem_version == 6) then  !M.Kazumori azimuth model function
        phi = Rel_Azimuth * DEGREES_TO_RADIANS
        wind10 = Wind_Speed
        theta = Zenith_Angle

        ! freq.
        DO ifreq = 1,6
          IF(wind10>18.0_fp) THEN
            ipol=1
            A1v(ifreq) = coef_mk_azi(1,ifreq,ipol) * ( exp(-coef_mk_azi(5,ifreq,ipol) * 18.0_fp * 18.0_fp ) - ONE ) * &
                      ( coef_mk_azi(2,ifreq,ipol) * 18.0_fp + coef_mk_azi(3,ifreq,ipol) * 18.0_fp * 18.0_fp + &
                        coef_mk_azi(4,ifreq,ipol) * 18.0_fp * 18.0_fp * 18.0_fp )
            A2v(ifreq) = coef_mk_azi(6,ifreq,ipol) * 18.0_fp

            ipol=2
            A1h(ifreq) = coef_mk_azi(1,ifreq,ipol) * 18.0_fp
            A2h(ifreq) = coef_mk_azi(2,ifreq,ipol) * ( exp(-coef_mk_azi(6,ifreq,ipol) * 18.0_fp * 18.0_fp ) - ONE ) * &
                      ( coef_mk_azi(3,ifreq,ipol) * 18.0_fp + coef_mk_azi(4,ifreq,ipol) * 18.0_fp * 18.0_fp + &
                        coef_mk_azi(5,ifreq,ipol) * 18.0_fp * 18.0_fp * 18.0_fp )
          ELSE
            ipol=1
            A1v(ifreq) = coef_mk_azi(1,ifreq,ipol) * ( exp(-coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 ) - ONE ) * &
                      ( coef_mk_azi(2,ifreq,ipol) * wind10 + coef_mk_azi(3,ifreq,ipol) * wind10 * wind10 + &
                        coef_mk_azi(4,ifreq,ipol) * wind10 * wind10 * wind10 )
            A2v(ifreq) = coef_mk_azi(6,ifreq,ipol) * wind10

            ipol=2
            A1h(ifreq) = coef_mk_azi(1,ifreq,ipol) * wind10
            A2h(ifreq) = coef_mk_azi(2,ifreq,ipol) * ( exp(-coef_mk_azi(6,ifreq,ipol) * wind10 * wind10 ) - ONE ) * &
                      ( coef_mk_azi(3,ifreq,ipol) * wind10 + coef_mk_azi(4,ifreq,ipol) * wind10 * wind10 + &
                        coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 * wind10 )
          END IF

          A1s1(ifreq) = (A1v(ifreq) + A1h(ifreq))/TWO
          A1s2(ifreq) =  A1v(ifreq) - A1h(ifreq)
          A2s1(ifreq) = (A2v(ifreq) + A2h(ifreq))/TWO
          A2s2(ifreq) =  A2v(ifreq) - A2h(ifreq)

          IF(Frequency>37.0_fp)THEN
            IF(wind10>15.0_fp)THEN
              A2s2_theta0(ifreq) = (15.0_fp*15.0_fp - 15.0_fp*15.0_fp*15.0_fp/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/37.0_fp) )
            ELSE
              A2s2_theta0(ifreq) = (wind10*wind10 - wind10*wind10*wind10/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/37.0_fp) )
            END IF
          ELSE
            IF(wind10>15.0_fp)THEN
              A2s2_theta0(ifreq) = (15.0_fp*15.0_fp - 15.0_fp*15.0_fp*15.0_fp/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/Frequency) )
            ELSE
              A2s2_theta0(ifreq) = (wind10*wind10 - wind10*wind10*wind10/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/Frequency) )
            END IF
          END IF

          A1s1_theta(ifreq)= A1s1(ifreq)*((theta/theta_ref)**xs11)
          A2s1_theta(ifreq)= A2s1(ifreq)*((theta/theta_ref)**xs12)
          A1s2_theta(ifreq)= A1s2(ifreq)*((theta/theta_ref)**xs21)
          A2s2_theta(ifreq)= A2s2_theta0(ifreq) + (A2s2(ifreq) - A2s2_theta0(ifreq))*((theta/theta_ref)**xs22)

          A1v_theta(ifreq) = 0.5d0*(2.d0*A1s1_theta(ifreq) + A1s2_theta(ifreq))
          A1h_theta(ifreq) = 0.5d0*(2.d0*A1s1_theta(ifreq) - A1s2_theta(ifreq))
          A2v_theta(ifreq) = 0.5d0*(2.d0*A2s1_theta(ifreq) + A2s2_theta(ifreq))
          A2h_theta(ifreq) = 0.5d0*(2.d0*A2s1_theta(ifreq) - A2s2_theta(ifreq))

          azimuth_component(1,ifreq) = A1v_theta(ifreq) * cos (real(1,8)*phi) + A2v_theta(ifreq) * cos(real(2,8)*phi)
          azimuth_component(2,ifreq) = A1h_theta(ifreq) * cos (real(1,8)*phi) + A2h_theta(ifreq) * cos(real(2,8)*phi)

        END DO

        IF( Frequency >= 1.4_fp .and. Frequency < 6.925_fp ) THEN
          Azimuth_Emi(1) = azimuth_component(1,1)
          Azimuth_Emi(2) = azimuth_component(2,1)
        ELSE IF( Frequency >= 6.925_fp .and. Frequency < 10.65_fp ) THEN
          fratio = ONE-(Frequency - 6.925_fp)/(10.65_fp - 6.925_fp)
          Azimuth_Emi(1) = azimuth_component(1,1)*fratio + (ONE-fratio)*azimuth_component(1,2)
          Azimuth_Emi(2) = azimuth_component(2,1)*fratio + (ONE-fratio)*azimuth_component(2,2)
        ELSE IF( Frequency > 10.65_fp .and. Frequency <= 18.7_fp ) THEN
          fratio = ONE-(Frequency - 10.65_fp)/(18.7_fp - 10.65_fp)
          Azimuth_Emi(1) = azimuth_component(1,2)*fratio + (ONE-fratio)*azimuth_component(1,3)
          Azimuth_Emi(2) = azimuth_component(2,2)*fratio + (ONE-fratio)*azimuth_component(2,3)
        ELSE IF( Frequency > 18.7_fp .and. Frequency <= 23.8_fp ) THEN
          fratio = ONE-(Frequency - 18.7_fp)/(23.8_fp - 18.7_fp)
          Azimuth_Emi(1) = azimuth_component(1,3)*fratio + (ONE-fratio)*azimuth_component(1,4)
          Azimuth_Emi(2) = azimuth_component(2,3)*fratio + (ONE-fratio)*azimuth_component(2,4)
        ELSE IF( Frequency > 23.8_fp .and. Frequency <= 36.5_fp ) THEN
          fratio = ONE-(Frequency - 23.8_fp)/(36.5_fp - 23.8_fp)
          Azimuth_Emi(1) = azimuth_component(1,4)*fratio + (ONE-fratio)*azimuth_component(1,5)
          Azimuth_Emi(2) = azimuth_component(2,4)*fratio + (ONE-fratio)*azimuth_component(2,5)
        ELSE IF( Frequency > 36.5_fp .and. Frequency <= 89.0_fp ) THEN
          fratio = ONE-(Frequency - 36.5_fp)/(89.0_fp - 36.5_fp)
          Azimuth_Emi(1) = azimuth_component(1,5)*fratio + (ONE-fratio)*azimuth_component(1,6)
          Azimuth_Emi(2) = azimuth_component(2,5)*fratio + (ONE-fratio)*azimuth_component(2,6)
        ELSE IF( Frequency > 89.0_fp .and. Frequency <= 200.0_fp ) THEN
          Azimuth_Emi(1) = azimuth_component(1,6)
          Azimuth_Emi(2) = azimuth_component(2,6)
        END IF
        ! Azimuthal component from Fuzhong Weng (NOAA/NESDIS) based on work by Dr. Gene Poe (NRL)
        ! FASTEM-3 azimuthal calculation for Stokes 3/4 only
        ! Assume 19m wind = 10m wind for now (fix later).
        u19 = wind10
        freq_ghz = Frequency
        c => fastem3_coef
        DO ich = 0, 15
          a1e = c(141 + ich * 12) + u19 * (c(142 + ich * 12) + u19 * (c(143 + ich * 12) + u19 * c(144 + ich * 12)))
          a2e = c(145 + ich * 12) + u19 * (c(146 + ich * 12) + u19 * (c(147 + ich * 12) + u19 * c(148 + ich * 12)))
          a3e = c(149 + ich * 12) + u19 * (c(150 + ich * 12) + u19 * (c(151 + ich * 12) + u19 * c(152 + ich * 12)))
          i_freq = int(ich / 4_jpim) + 1                      ! 37, 19, 10, 7 GHz
          istokes                     = mod(ich, 4_jpim) + 1
          tbfixed(istokes, i_freq, 1) = a1e                  !* prof%skin%t
          tbfixed(istokes, i_freq, 2) = a2e                  !* prof%skin%t
          tbfixed(istokes, i_freq, 3) = a3e                  !* prof%skin%t
        ENDDO
        DO M = 1, 3
          DO istokes = 3, 4
            efixed(1, istokes, M) = tbfixed(istokes, 4, M)!/prof%skin%t  ! 7  GHz
            efixed(2, istokes, M) = tbfixed(istokes, 3, M)!/prof%skin%t  ! 10  GHz
            efixed(3, istokes, M) = tbfixed(istokes, 2, M)!/prof%skin%t  ! 19  GHz
            efixed(4, istokes, M) = tbfixed(istokes, 1, M)!/prof%skin%t  ! 37  GHz
          ENDDO
          ! Interpolate results to required frequency based on 7, 10, 19, 37 GHz
          IF (freq_ghz .LE. freqfixed(1)) THEN
            DO istokes = 3, 4
              einterpolated(istokes, M) = efixed(1, istokes, M)
            ENDDO
          ELSE IF (freq_ghz .GE. freqfixed(4)) THEN
            DO istokes = 3, 4
              einterpolated(istokes, M) = efixed(4, istokes, M)
            ENDDO
          ELSE
            IF (freq_ghz .LT. freqfixed(2)                                 ) ifreq = 2
            IF (freq_ghz .LT. freqfixed(3) .AND. freq_ghz .GE. freqfixed(2)) ifreq = 3
            IF (freq_ghz .GE. freqfixed(3)                                 ) ifreq = 4
            dfreq = (freq_ghz - freqfixed(ifreq - 1)) / (freqfixed(ifreq) - freqfixed(ifreq - 1))
            DO istokes = 3, 4
              einterpolated(istokes, M) =      &
                & efixed(ifreq - 1, istokes, M) + dfreq * (efixed(ifreq, istokes, M) - efixed(ifreq - 1, istokes, M))
            ENDDO
          ENDIF
        ENDDO
        DO istokes = 3, 4
          DO M = 1, 3
            Azimuth_Emi(istokes) = Azimuth_Emi(istokes) +      &
              & einterpolated(istokes, M) * sin(m * phi) * (1.0_JPRB - cos_z) / (1.0_JPRB - 0.6018_JPRB)
          ENDDO
        ENDDO

      else ! FASTEM-4/5    M.Liu      azimuth model function
        Fre_C = ZERO
        IF( Frequency >= min_f .or. Frequency <= max_f ) THEN
          DO i = 1, 8
            IF( Frequency >= x(i) .and. Frequency < x(i+1) ) THEN
              Fre_C = y(i) + (y(i+1)-y(i))/(x(i+1)-x(i))*(Frequency-x(i))
            END IF
          END DO
        END IF

        phi = Rel_Azimuth * DEGREES_TO_RADIANS
        wind10 = WInd_Speed
        DO m = 1, 3
          L = 10*(m-1)
          ac = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(1) = Azimuth_Emi(1) + ac*cos(m*phi)

          L = 10*(m-1) + 30
          ac = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(2) = Azimuth_Emi(2) + ac*cos(m*phi)

          L = 10*(m-1) + 60
          sc = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(3) = Azimuth_Emi(3) + sc*sin(m*phi)

          L = 10*(m-1) + 90
          sc = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(4) = Azimuth_Emi(4) + sc*sin(m*phi)
        END DO

        Azimuth_Emi = Azimuth_Emi * Fre_C
      endif

    END IF

    Emissivity(1) = Emissivity(1) + Azimuth_Emi(1)
    Emissivity(2) = Emissivity(2) + Azimuth_Emi(2)
    Emissivity(3) = Azimuth_Emi(3)
    Emissivity(4) = Azimuth_Emi(4)

    Reflectivity(1) = zreflmod_v * (ONE-Emissivity(1))
    Reflectivity(2) = zreflmod_h * (ONE-Emissivity(2))
    IF (fastem_version == 6) THEN
      Reflectivity(3) = -0.5_jprb * (zreflmod_v + zreflmod_h) * Emissivity(3)
      Reflectivity(4) = -0.5_jprb * (zreflmod_v + zreflmod_h) * Emissivity(4)
    ELSE
      ! FASTEM-4/5 set Stokes 3/4 refl to zero: this is not recommended
      Reflectivity(3:4) = ZERO
    ENDIF

   RETURN

  END SUBROUTINE rttov_fastem5
!
