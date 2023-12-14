! Description:
!> @file
!!   Subroutines for SURFEM-Ocean MW sea surface emissivity model
!
!> @brief
!!   Subroutines for SURFEM-Ocean MW sea surface emissivity model
!!
!! @details
!!   This contains the code which implements SURFEM-Ocean for the direct model.
!!
!!   SURFEM-Ocean is a neural network-based emissivity model applicable between 500 MHz
!!   and 700GHz
!!
!!   It is recommmended to use SURFEM-Ocean for channels below 1.4 GHz and above 200GHz. 
!!
!!   Reference:
!!   Kilic, L., Prigent, C., Jimenez. C. (2022)
!!   Development of the SURface Fast Emissivity Model for Ocean (SURFEM-Ocean)
!!   NWP SAF associate scientist mission NWP_AS20_01
!!
!!   Meissner, T., & Wentz, F. J. (2004). The complex dielectric constant of pure
!!   and sea water from microwave satellite observations. IEEE Transactions on
!!   Geoscience and Remote Sensing, 42(9), 18361849.
!!   https://doi.org/10.1109/TGRS.2004.831888
!!
!!   Meissner, T., & Wentz, F. J. (2012). The Emissivity of the Ocean Surface
!!   Between 6 and 90 GHz Over a Large Range of Wind Speeds and Earth Incidence
!!   Angles. IEEE Transactions on Geoscience and Remote Sensing, 50(8), 30043026.
!!   https://doi.org/10.1109/TGRS.2011.2179662
!!
!!   Meissner, T., Wentz, F. J., & Ricciardulli, L. (2014). The emission and
!!   scattering of L-band microwave radiation from rough ocean surfaces and wind
!!   speed measurements from the Aquarius sensor. Journal of Geophysical Research C:
!!   Oceans, 119(9), 64996522. https://doi.org/10.1002/2014JC009837
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
!
MODULE rttov_surfem_ocean_mod

  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : pi

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rttov_surfem_ocean
  PUBLIC :: epsilon, e_neutral, net_surfem

CONTAINS

  !> Return H-pol , V-pol and S3 and S4 emissivities for given frequency, zenith angle, 
  !! Tskin, salinity, wind speed and reletive wind direction
  !! @param[in]     f_GHz                 frequency (GHz)
  !! @param[in]     theta                 zenith angle (degrees)
  !! @param[in]     windspeed             wind speed (m/s)
  !! @param[in]     SST                   skin temperature (K)
  !! @param[in]     SSS                   salinity (practical salinity unit)
  !! @param[in]     phi                   relative wind direction (degrees)
  !! @param[in]     transmittance         surface-to-space transmittance
  !! @param[out]    emissivity            calculated emissivity (4 Stokes components)
  !! @param[out]    reflectivity          calculated reflectivity (4 Stokes components)
  SUBROUTINE rttov_surfem_ocean(f_GHz, theta, windspeed, SST, SSS, phi, transmittance, emissivity, reflectivity)

    USE rttov_surfem_ocean_coef_mod, ONLY : &
        surfem_nn_min_sst,       &
        surfem_nn_max_sst,       &
        surfem_nn_min_sss,       &
        surfem_nn_max_sss,       &
        surfem_nn_min_windspeed, &
        surfem_nn_max_windspeed, &
        surfem_nn_min_freq,      &
        surfem_nn_max_freq     
    USE mod_rttov_fastem5_coef, ONLY : &
        t_c5, transmittance_limit_lower, transmittance_limit_upper, DEGREES_TO_RADIANS
    IMPLICIT NONE
    REAL(jprb), INTENT(IN)  :: f_GHz, theta, windspeed, SST, SSS, phi, transmittance
    REAL(jprb), INTENT(OUT) :: emissivity(4), reflectivity(4)

    ! Local copies of input variables clipped to neural network training limits
    REAL(jprb) :: SSTl, SSSl, windspeedl, f_GHzl
    INTEGER :: iz, jz
    COMPLEX(jprb) :: dielec     ! complex dielectric constant of sea water
    REAL(jprb)    :: e_Nv       ! vertical component of e Neutral
    REAL(jprb)    :: e_Nh       ! horizontal component of e Neutral
    REAL(jprb)    :: e_v0       ! vertical component of e isotropic
    REAL(jprb)    :: e_h0       ! horizontal component of e isotropic
    REAL(jprb)    :: e_iso(2)   ! e isotropic
    REAL(jprb)    :: e_aniso(8) ! e anisotropic

    REAL(jprb) :: zreflmod_v,zreflmod_h,zrough_v,zrough_h
    ! Local variables for including transmittance
    REAL(jprb) :: cos_z, variance,varm,opdpsfc,zx(9)
    ! calculate dielectric constant of sea water
    CALL epsilon(f_GHz, SST, SSS, dielec)

    ! calculate flat sea emissivity (Neutral)
    CALL e_neutral(dielec, theta, e_Nv, e_Nh)

    ! Ensure NN training limits are not exceeded
    f_GHzl     = MAX(MIN(f_GHz, surfem_nn_max_freq), surfem_nn_min_freq)
    SSTl       = MAX(MIN(SST, surfem_nn_max_sst), surfem_nn_min_sst)
    SSSl       = MAX(MIN(SSS, surfem_nn_max_sss), surfem_nn_min_sss)
    windspeedl = MAX(MIN(windspeed, surfem_nn_max_windspeed), surfem_nn_min_windspeed)

    ! calculate isotropic emissivity via neural network
    CALL net_surfem(f_GHzl, theta, windspeedl, SSTl, SSSl, e_Nv, e_Nh, e_iso)
    e_v0 = e_iso(1)
    e_h0 = e_iso(2)

    ! calculate anisotropic emissivity via neural network
    CALL net_surfem(f_GHzl, theta, windspeedl, SSTl, SSSl, e_Nv, e_Nh, e_aniso)

    ! calculate output emissitivities
    CALL e_total(e_Nv, e_Nh, e_v0, e_h0, e_aniso, phi, emissivity(1), emissivity(2), emissivity(3), emissivity(4))

    cos_z = cos( theta*DEGREES_TO_RADIANS )

    ! correction term for reflectivity same as FASTEM 5/6 - rearrange later ------------------
    IF( transmittance > transmittance_limit_lower .and. transmittance < transmittance_limit_upper) THEN
      !Using the Cox and Munk model to compute slope variance
      variance = 0.00512_JPRB * windspeed + 0.0030_JPRB
      varm     = variance * t_c5(43)
      variance = varm * ( t_c5(44) * f_GHz + t_c5(45) )
      IF ( variance >= varm ) variance = varm
      IF ( variance <= 0.0_JPRB ) variance = 0.0_JPRB
      !Compute surface to space optical depth
      opdpsfc = -log(transmittance ) * cos_z ! note cos_z is recalculated in fastem5 sub so bit diff
      !Define nine predictors for the effective angle calculation
      zx(1) = 1.0_JPRB
      zx(2) = variance
      zx(4) = 1.0_JPRB / cos_z
      zx(3) = zx(2) * zx(4)
      zx(5) = zx(3) * zx(3)
      zx(6) = zx(4) * zx(4)
      zx(7) = zx(2) * zx(2)
      zx(8) = log(opdpsfc)
      zx(9) = zx(8) * zx(8)

      zrough_v = 1.0_JPRB
      zrough_h = 1.0_JPRB
      DO iz = 1, 7
        jz = iz-1
        !Switched h to v Deblonde SSMIS june 7, 2001
        zrough_h = zrough_h + zx(iz) *(t_c5(1+jz*3) + zx(8)*t_c5(2+jz*3) + zx(9)*t_c5(3+jz*3) )
        zrough_v = zrough_v + zx(iz) *(t_c5(22+jz*3)+ zx(8)*t_c5(23+jz*3)+ zx(9)*t_c5(24+jz*3))
      END DO
      zreflmod_v = (1.0_JPRB-transmittance ** zrough_v)/(1.0_JPRB-transmittance )
      zreflmod_h = (1.0_JPRB-transmittance ** zrough_h)/(1.0_JPRB-transmittance )
    ELSE
      zreflmod_v = 1._jprb
      zreflmod_h = 1._jprb
    END IF
    ! end correction term -----------------------------------------------------------

    reflectivity(1) = zreflmod_v * (1.0_JPRB - emissivity(1))
    reflectivity(2) = zreflmod_h * (1.0_JPRB - emissivity(2))
    reflectivity(3) = -0.5_JPRB * (zreflmod_v + zreflmod_h) * emissivity(3)
    reflectivity(4) = -0.5_JPRB * (zreflmod_v + zreflmod_h) * emissivity(4)

  END SUBROUTINE rttov_surfem_ocean

  SUBROUTINE epsilon(f_GHz, SSTk, SSSi, dielec)

    USE rttov_surfem_ocean_coef_mod, ONLY : &
        f0, es_coef, ai_coef,       &
        sigma35_coef, R15_coef,     &
        alpha0_coef, alpha1_coef,   &
        bi_new1_coef, bi_new2_coef, &
        c1_coef, bi_coef,           &
        surfem_dielec_min_sst,      &
        surfem_dielec_max_sst,      &
        surfem_dielec_min_sss,      &
        surfem_dielec_max_sss

    USE rttov_const, ONLY : t0

    ! Computes the dielectric constant of sea water according to the model
    ! by Meissner and Wentz (2004) including corrections by Meissner and Wentz
    ! 2012 and Meissner et al. 2014.
    ! Translated from matlab to fortran by E. Turner (2022) 
    !
    !       Inputs:
    !               SSTk : Sea surface temperature in kelvin
    !               SSSi : Sea surface salinity in psu
    !               f_GHz    : frequency in Hz
    !
    !       Outputs:
    !               dielec : complex dielectric constant of sea water
    !
    REAL(jprb),       INTENT(IN)    :: f_GHz    !frequency in GHz
    REAL(jprb),       INTENT(IN)    :: SSTk     !Sea surface temperature in Kelvin
    REAL(jprb),       INTENT(IN)    :: SSSi     !Sea surface salinity in psu
    COMPLEX(jprb),    INTENT(OUT)   :: dielec   !complex dielectric constant of sea water

    REAL(jprb) :: SST, SSS, SST2, SST3, SST4, SSS2
    REAL(jprb) :: es, e1, einf, nu1, nu2 ! Fresh water parameters
    REAL(jprb) :: sigma35, R15, alpha0, alpha1, RTR15, sigma  ! Salt water
    REAL(jprb) :: es_s, e1_s, einf_s, nu1_s, nu2_s ! Salt water parameters
    REAL(jprb) :: c1, c2, c3, c4 ! conversion coef from fresh to salt water
    COMPLEX(jprb), PARAMETER :: i = (0, 1)   ! sqrt(-1) 

    ! clip input variables if they are outside of reasonable limits
    SST = MAX(MIN(SSTk, surfem_dielec_max_sst), surfem_dielec_min_sst) - t0 ! K -> deg C
    SSS = MAX(MIN(SSSi, surfem_dielec_max_sss), surfem_dielec_min_sss)

    ! combined local variables ??
    SST2 = SST * SST
    SST3 = SST2 * SST
    SST4 = SST3 * SST
    SSS2 = SSS * SSS

    ! Fresh water parameters
    ! es:  Static dielectric constant for pure water by Stogryn et al. 1995
    es = (es_coef(1) - es_coef(2) * SST) / (es_coef(3) + SST)
    e1 = ai_coef(1) + ai_coef(2) * SST + ai_coef(3) * SST2
    nu1 = (45.0_jprb + SST) / ( ai_coef(4) + (ai_coef(5) * SST) + ai_coef(6)*SST2 )
    einf = ai_coef(7) + (ai_coef(8) * SST)
    nu2 = (45.0_jprb + SST) / ( ai_coef(9) + (ai_coef(10) * SST) + (ai_coef(11)*SST2) )

    ! Salt water
    ! Conductivity [s/m] by Stogryn et al. 1995
    sigma35 = sigma35_coef(1) + (sigma35_coef(2) * SST) + (sigma35_coef(3) * SST2) - &
          (sigma35_coef(4) * SST3) + (sigma35_coef(5)*SST4)
    R15 = SSS * ( R15_coef(1) + (R15_coef(2) * SSS) + (R15_coef(3) * SSS2) ) / & 
          ( R15_coef(4) + (R15_coef(5) * SSS) + SSS2 )
    alpha0 = ( alpha0_coef(1) + (alpha0_coef(2) * SSS) - (alpha0_coef(3) * SSS2) ) / &
          ( alpha0_coef(4) + (alpha0_coef(5) * SSS) + SSS2 )
    alpha1 =  alpha1_coef(1) - (alpha1_coef(2) * SSS) + alpha1_coef(3) * SSS2
    RTR15 = 1.0_jprb + ( SST - 15.0_jprb ) * alpha0 / ( alpha1 + SST )
    sigma = sigma35 * R15 * RTR15

    ! Salt water parameters 
    es_s = es * EXP ( (bi_new1_coef(1) * SSS) + (bi_new1_coef(2) * SSS2) + (bi_new1_coef(3) * SSS * SST ) )

    ! conversion coefs from fresh to salt water

    IF ( SST <= 30._jprb ) THEN
      c1 = 1.0_jprb + SSS * ( bi_new2_coef(1) + (bi_new2_coef(2) * SST) + (bi_new2_coef(3) * SST2) + &
            (bi_new2_coef(4) * SST3) + (bi_new2_coef(5) * SST4) )
    ELSE
      c1 = 1.0_jprb + SSS * ( c1_coef(1) + c1_coef(2) * (SST - 30.0_jprb) )
    ENDIF

    nu1_s = nu1*c1

    c2  = EXP((bi_coef(7)*SSS) + (bi_coef(8)*SSS2) + (bi_coef(9)*SSS*SST)) 

    e1_s = e1*c2

    !  c3 = 1.0_jprb + SSS *_coef(bi_coef(10) + bi_coef(11)*SST) ! Expression for c3  amended in Meissner et al. 2016
    c3 = 1.0_jprb + SSS * ( bi_coef(10) + (0.5_jprb * bi_coef(11)) * (SST + 30.0_jprb) )

    nu2_s = nu2 * c3

    c4 = 1.0_jprb  + SSS * ( bi_coef(12) + (bi_coef(13) * SST) )

    einf_s = einf * c4  

    ! dielectric returned. Original complex form.
    ! dielec = ( es_s - e1_s ) / ( 1.0_jprb - i * ( f_GHz/nu1_s ) ) + &
    !          ( e1_s - einf_s ) / ( 1.0_jprb - i * ( f_GHz/nu2_s ) ) + &
    !          einf_s + i*sigma*f0/f_GHz 

    ! multiply by complex conjugate (1.0_jprb + i * (f_GHz/nu1_s)) to speed up calculation 
    dielec = (es_s - e1_s) * (1._jprb + i * f_GHz/nu1_s) / (1._jprb + (f_GHz/nu1_s)**2) + &
              (e1_s - einf_s) * ( 1._jprb + i * f_GHz/nu2_s) / (1._jprb + (f_GHz/nu2_s)**2) + &
              einf_s + &
              i*sigma*f0/f_GHz

  END SUBROUTINE epsilon

  ! Estimation of e neutral (flat sea no wind) from dielectric constant
  SUBROUTINE e_neutral(dielec, theta, e_Nv, e_Nh)

    COMPLEX(jprb),  INTENT(IN)  :: dielec  !complex dielectric constant of sea water
    REAL(jprb),     INTENT(IN)  :: theta     !incidence angle in degrees
    REAL(jprb),     INTENT(OUT) :: e_Nv   ! vertical component of e Neutral
    REAL(jprb),     INTENT(OUT) :: e_Nh   ! horizontal component of e Neutral

    COMPLEX(jprb) :: f_v ! vertical fresnel reflection coefficient unsquared
    COMPLEX(jprb) :: f_h ! horizontal fresnel reflection coefficient unsquared
    REAL(jprb) :: fres_v ! vertical fresnel reflection coefficient
    REAL(jprb) :: fres_h ! horizontal fresnel reflection coefficient
    REAL(jprb)    :: cos_theta
    REAL(jprb)    :: sin_theta
    COMPLEX(jprb) :: sqrt_d

    cos_theta = COS(theta*pi/180._jprb) 
    sin_theta = SIN(theta*pi/180._jprb)

    sqrt_d = SQRT(dielec - sin_theta * sin_theta)

    f_v = (dielec * cos_theta - sqrt_d) / &
          (dielec * cos_theta + sqrt_d) 

    f_h = (cos_theta - sqrt_d) / &
          (cos_theta + sqrt_d) 

    fres_v = REAL(f_v)**2 + AIMAG(f_v)**2
    fres_h = REAL(f_h)**2 + AIMAG(f_h)**2

    e_Nv = 1.0_jprb - fres_v
    e_Nh = 1.0_jprb - fres_h

  END SUBROUTINE e_neutral

  ! NEURAL NETWORK PROCEDURES 
  ! X is a vector of [f_GHz_vec,theta_vec,windspeed_vec,SST_vec,SSS_vec,evn,ehn]
  ! 1. Input X transformed with mapminmax = Process matrices by mapping row minimum and maximum &
  !    values to [-1 1] => aux = iymi + (X - ioff).* igai [y = (ymax-ymin)*(x-xmin)/(xmax-xmin) + ymin]
  ! 2. Propagate through the first layer => x1  = W1 * aux + b1
  ! 3. Transform with tansig(tanh) => aux = 2 ./ (1+exp(-2 * x1 ) ) -1 (activation function)
  ! 4. Propagate through output layer => aux = W2 * aux + b2 (needed?)
  ! 5. Output Y transformed back with mapminmax => Y = ((aux - oymi ) ./ ogai ) + ooff

  !> Return emissivity components for given frequency, zenith angle, T-skin
  !! salinity and wind speed
  !! @param[in]     f_GHz       frequency (GHz)
  !! @param[in]     theta       zenith angle (degrees)
  !! @param[in]     windspeed   wind speed (m/s)
  !! @param[in]     tskin       skin temperature (K)
  !! @param[in]     salinity    salinity (practical salinity unit)
  !! @param[in]     e_Nv        vertical component of e Neutral
  !! @param[in]     e_Nh        horizontal component of e Neutral
  !! @param[out]    e_out       calculated emissivity components
  SUBROUTINE net_surfem(f_GHz, theta, windspeed, tskin, salinity, e_Nv, e_Nh, e_out)

    USE rttov_surfem_ocean_coef_mod, ONLY : &
        surfem_nh, surfem_ni, surfem_no_iso, iso_net, aniso_net

    REAL(jprb), INTENT(IN)  :: f_GHz, theta, windspeed, tskin, salinity, e_Nv, e_Nh
    REAL(jprb), INTENT(OUT) :: e_out(:)

    REAL(jprb) :: x(surfem_ni)
    REAL(jprb) :: tran_1(surfem_ni)
    REAL(jprb) :: prop_1(surfem_nh)
    REAL(jprb) :: prop_2(surfem_nh)

    REAL(jprb), POINTER :: ti_ymin(:), ti_xoffset(:), ti_gain(:)
    REAL(jprb), POINTER :: b1(:), w1(:,:), b2(:), w2(:,:)
    REAL(jprb), POINTER :: to_ymin(:), to_xoffset(:), to_gain(:)

    INTEGER(jpim) :: i, j, surfem_no

    x(:) = (/ f_GHz, theta, windspeed, tskin, salinity, e_Nv, e_Nh /)

    surfem_no = SIZE(e_out)

    IF (surfem_no == surfem_no_iso) THEN
      ti_ymin    => iso_net%ti_ymin
      ti_xoffset => iso_net%ti_xoffset
      ti_gain    => iso_net%ti_gain

      b1         => iso_net%b1
      w1         => iso_net%w1
      b2         => iso_net%b2
      w2         => iso_net%w2

      to_ymin    => iso_net%to_ymin
      to_xoffset => iso_net%to_xoffset
      to_gain    => iso_net%to_gain
    ELSE
      ti_ymin    => aniso_net%ti_ymin
      ti_xoffset => aniso_net%ti_xoffset
      ti_gain    => aniso_net%ti_gain

      b1         => aniso_net%b1
      w1         => aniso_net%w1
      b2         => aniso_net%b2
      w2         => aniso_net%w2

      to_ymin    => aniso_net%to_ymin
      to_xoffset => aniso_net%to_xoffset
      to_gain    => aniso_net%to_gain
    ENDIF

    ! Step 1: Transform input 
    DO i = 1, surfem_ni
      tran_1(i)= ti_ymin(1) + (x(i) - ti_xoffset(i)) * ti_gain(i)
    ENDDO

    ! Step 2: Propagate through the first layer => x1  = W1 * aux + b1
    DO i = 1, surfem_nh
      prop_1(i) = b1(i)
      DO j = 1, surfem_ni ! sum up all input parameters for each node 
        prop_1(i) = prop_1(i) + W1(j,i) * tran_1(j)  
      ENDDO 
      ! Step 3: Transform with tansig(tanh) - would TANH be better to use?
      prop_1(i) = 2._jprb / (1._jprb + EXP(-2._jprb * prop_1(i))) - 1._jprb
    ENDDO

    ! Step 4: Propagate through the second layer => aux  = W2 * aux + b2
    DO i = 1, surfem_no
      prop_2(i) = b2(i)
      DO j = 1, surfem_nh
        prop_2(i) = prop_2(i) + W2(j,i) * prop_1(j)
      ENDDO
    ENDDO

    ! Step 5. Output Y transformed back with mapminmax => Y = ((aux - oymi ) ./ ogai ) + ooff
    DO i = 1, surfem_no
      e_out(i) = ((prop_2(i) - to_ymin(1)) / to_gain(i)) + to_xoffset(i) 
    ENDDO

  END SUBROUTINE net_surfem

  SUBROUTINE e_total(e_Nv, e_Nh, e_v0, e_h0, e_aniso, phi, ev, eh, e3, e4)

    USE rttov_surfem_ocean_coef_mod, ONLY : surfem_no_aniso

    REAL(jprb), INTENT(IN)  :: e_Nv, e_Nh, e_v0, e_h0, e_aniso(surfem_no_aniso)
    REAL(jprb), INTENT(IN)  :: phi ! relative wind direction in degrees
    REAL(jprb), INTENT(OUT) :: ev, eh, e3, e4

    REAL(jprb) :: e_v1, e_h1, e_31, e_41, e_v2, e_h2, e_32, e_42
    REAL(jprb) :: phi_r ! phi in radians

    e_v1 = e_aniso(1)
    e_h1 = e_aniso(2)
    e_31 = e_aniso(3)
    e_41 = e_aniso(4)
    e_v2 = e_aniso(5)
    e_h2 = e_aniso(6)
    e_32 = e_aniso(7)
    e_42 = e_aniso(8)

    phi_r = phi * pi / 180._jprb 

    ev = e_Nv + e_v0 + e_v1 * COS(phi_r) + e_v2 * COS(2.0_jprb*phi_r)
    eh = e_Nh + e_h0 + e_h1 * COS(phi_r) + e_h2 * COS(2.0_jprb*phi_r)
    e3 = e_31 * SIN(phi_r) + e_32 * SIN(2.0_jprb*phi_r)
    e4 = e_41 * SIN(phi_r) + e_42 * SIN(2.0_jprb*phi_r)

  END SUBROUTINE e_total

END MODULE rttov_surfem_ocean_mod
