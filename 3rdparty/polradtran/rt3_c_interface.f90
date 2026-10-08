! ARTS3: ISO_C_BINDING entry points for the RT3 library.
!
! The Fortran 77 sources are compiled with default (4-byte) INTEGER and
! with REAL*8 / COMPLEX*16 arguments.  These wrappers take explicit C
! types (int64_t scalars and the one-character flags by value, int64_t
! index arrays, doubles, the complex ground index as two doubles),
! convert them to the default Fortran kinds and call the original
! routines.  This removes the need for -fdefault-integer-8 and for the
! compiler-specific hidden string length arguments of CHARACTER dummies.
!
! The C names (rt3_radtran, rt3_*_quadrature, and rt3_ and the name of
! every other routine that RADTRAN or a routine ported to C++ calls,
! without its RT3_ prefix) differ from
! the linker names of the Fortran routines (rt3_*_quadrature_, with the
! trailing underscore of gfortran's external names).  The Fortran names of
! the wrappers are c_rt3_*.
!
! The real arrays are passed through unchanged.  Their layouts are
! documented at the RADTRAN header in radtran3.f and in
! src/core/polradtran/rt3/rt3.h.  None of these routines is reentrant: RT3 uses
! COMMON blocks, SAVEd FFT tables and large static local arrays, so the
! caller must serialise all calls.
module rt3_c_interface
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_double, c_char, c_bool
  implicit none
  private

  interface
    subroutine radtran(nstokes, nummu, aziorder, max_delta_tau, &
                       src_code, quad_type, deltam, direct_flux, &
                       direct_mu, ground_temp, ground_type, &
                       ground_albedo, ground_index, sky_temp, wavelength, &
                       num_layers, height, temperatures, gas_extinct, &
                       nsl, scat_extinct, scat_scatter, scat_nlegen, &
                       ldcoef, scat_coef, scatlayers, noutlevels, &
                       outlevels, mu_values, up_flux, down_flux, up_rad, &
                       down_rad)
      integer :: nstokes, nummu, aziorder, src_code, num_layers
      integer :: nsl, ldcoef, noutlevels
      integer :: scat_nlegen(*), scatlayers(*), outlevels(*)
      real(8) :: max_delta_tau, direct_flux, direct_mu
      real(8) :: ground_temp, ground_albedo, sky_temp, wavelength
      complex(8) :: ground_index
      character(len=1) :: quad_type, deltam, ground_type
      real(8) :: height(*), temperatures(*), gas_extinct(*)
      real(8) :: scat_extinct(*), scat_scatter(*), scat_coef(*)
      real(8) :: mu_values(*), up_flux(*), down_flux(*)
      real(8) :: up_rad(*), down_rad(*)
    end subroutine radtran

    subroutine rt3_double_gauss_quadrature(num, abscissas, weights)
      integer :: num
      real(8) :: abscissas(*), weights(*)
    end subroutine rt3_double_gauss_quadrature

    subroutine rt3_gauss_legendre_quadrature(num, abscissas, weights)
      integer :: num
      real(8) :: abscissas(*), weights(*)
    end subroutine rt3_gauss_legendre_quadrature

    subroutine rt3_lobatto_quadrature(num, abscissas, weights)
      integer :: num
      real(8) :: abscissas(*), weights(*)
    end subroutine rt3_lobatto_quadrature

    ! The routines RADTRAN calls, for the C++ port of RADTRAN
    ! (src/core/polradtran/rt3/radtran3.cc)
    subroutine get_scat_set(deltam, nummu, nlegin, coefin, extin, scatin, &
        nlegen, coef, extinction, scatter)
      integer :: nummu, nlegin, nlegen
      real(8) :: extin, scatin, extinction, scatter
      character(len=1) :: deltam
      real(8) :: coefin(*), coef(*)
    end subroutine get_scat_set

    subroutine scattering(nummu, aziorder, nstokes, mu_values, quad_weights, &
        numlegendre, legendre_coef, scat_num, scatbuf)
      integer :: nummu, aziorder, nstokes, numlegendre, scat_num
      real(8) :: mu_values(*), quad_weights(*), legendre_coef(*), scatbuf(*)
    end subroutine scattering

    subroutine direct_scattering(nummu, aziorder, nstokes, mu_values, &
        numlegendre, legendre_coef, direct_mu, scat_num, directbuf)
      integer :: nummu, aziorder, nstokes, numlegendre, scat_num
      real(8) :: direct_mu
      real(8) :: mu_values(*), legendre_coef(*), directbuf(*)
    end subroutine direct_scattering

    subroutine get_scattering(nstokes, nummu, mode, aziorder, scat_num, &
        scatbuf, scatter_matrix)
      integer :: nstokes, nummu, mode, aziorder, scat_num
      real(8) :: scatbuf(*), scatter_matrix(*)
    end subroutine get_scattering

    subroutine rt3_check_norm(nstokes, nummu, quad_weights, scatter_matrix)
      integer :: nstokes, nummu
      real(8) :: quad_weights(*), scatter_matrix(*)
    end subroutine rt3_check_norm

    subroutine get_direct(nstokes, nummu, mode, aziorder, scat_num, directbuf, &
        direct_vector)
      integer :: nstokes, nummu, mode, aziorder, scat_num
      real(8) :: directbuf(*), direct_vector(*)
    end subroutine get_direct

    subroutine rt3_thermal_radiance(nstokes, nummu, mode, temperature, albedo, &
        wavelength, radiance)
      integer :: nstokes, nummu, mode
      real(8) :: temperature, albedo, wavelength
      real(8) :: radiance(*)
    end subroutine rt3_thermal_radiance

    subroutine rt3_nonscatter_layer(nstokes, nummu, mode, deltatau, mu_values, &
        planck0, planck1, reflect, trans, source)
      integer :: nstokes, nummu, mode
      real(8) :: deltatau, planck0, planck1
      real(8) :: mu_values(*), reflect(*), trans(*), source(*)
    end subroutine rt3_nonscatter_layer

    subroutine rt3_initial_source(nstokes, nummu, n, delta_z, mu_values, &
        extinction, source_vector, source)
      integer :: nstokes, nummu, n
      real(8) :: delta_z, extinction
      real(8) :: mu_values(*), source_vector(*), source(*)
    end subroutine rt3_initial_source

    subroutine rt3_initialize(nstokes, nummu, n, delta_z, mu_values, &
        extinction, albedo, phase_function, reflect, trans)
      integer :: nstokes, nummu, n
      real(8) :: delta_z, extinction, albedo
      real(8) :: mu_values(*), phase_function(*), reflect(*), trans(*)
    end subroutine rt3_initialize

    subroutine rt3_doubling_integration(n, num_doubles, src_code, symmetric, &
        reflect, trans, exp_source, expfactor, lin_source, linfactor, &
        t_reflect, t_trans, t_source)
      integer :: n, num_doubles, src_code
      real(8) :: expfactor, linfactor
      logical :: symmetric
      real(8) :: reflect(*), trans(*), exp_source(*), lin_source(*), &
          t_reflect(*), t_trans(*), t_source(*)
    end subroutine rt3_doubling_integration

    subroutine rt3_combine_layers(n, reflect1, trans1, source1, reflect2, &
        trans2, source2, out_reflect, out_trans, out_source)
      integer :: n
      real(8) :: reflect1(*), trans1(*), source1(*), reflect2(*), trans2(*), &
          source2(*), out_reflect(*), out_trans(*), out_source(*)
    end subroutine rt3_combine_layers

    subroutine rt3_internal_radiance(n, upreflect, uptrans, upsource, &
        downreflect, downtrans, downsource, intoprad, inbottomrad, uprad, &
        downrad)
      integer :: n
      real(8) :: upreflect(*), uptrans(*), upsource(*), downreflect(*), &
          downtrans(*), downsource(*), intoprad(*), inbottomrad(*), uprad(*), &
          downrad(*)
    end subroutine rt3_internal_radiance

    subroutine rt3_lambert_surface(nstokes, nummu, mode, mu_values, &
        quad_weights, ground_albedo, reflect, trans, source)
      integer :: nstokes, nummu, mode
      real(8) :: ground_albedo
      real(8) :: mu_values(*), quad_weights(*), reflect(*), trans(*), &
          source(*)
    end subroutine rt3_lambert_surface

    subroutine rt3_lambert_radiance(nstokes, nummu, mode, src_code, &
        ground_albedo, ground_temp, wavelength, direct_sfc_flux, radiance)
      integer :: nstokes, nummu, mode, src_code
      real(8) :: ground_albedo, ground_temp, wavelength, direct_sfc_flux
      real(8) :: radiance(*)
    end subroutine rt3_lambert_radiance

    subroutine rt3_fresnel_surface(nstokes, nummu, mu_values, index, reflect, &
        trans, source)
      integer :: nstokes, nummu
      complex(8) :: index
      real(8) :: mu_values(*), reflect(*), trans(*), source(*)
    end subroutine rt3_fresnel_surface

    subroutine rt3_fresnel_radiance(nstokes, nummu, mode, mu_values, index, &
        ground_temp, wavelength, radiance)
      integer :: nstokes, nummu, mode
      real(8) :: ground_temp, wavelength
      complex(8) :: index
      real(8) :: mu_values(*), radiance(*)
    end subroutine rt3_fresnel_radiance

    ! The routines SCATTERING calls, for its C++ port
    ! (src/core/polradtran/rt3/radscat3.cc)
    subroutine number_sums(nstokes, nlegen, coef, dosum)
      integer :: nstokes, nlegen, dosum(6)
      real(8) :: coef(*)
    end subroutine number_sums

    subroutine sum_legendre(nlegen, coef, x, dosum, phase_matrix)
      integer :: nlegen, dosum(6)
      real(8) :: x
      real(8) :: coef(*), phase_matrix(*)
    end subroutine sum_legendre

    subroutine rotate_phase_matrix(phase_matrix1, mu1, mu2, delphi, cos_scat, &
        phase_matrix2, nstokes)
      integer :: nstokes
      real(8) :: mu1, mu2, delphi, cos_scat
      real(8) :: phase_matrix1(*), phase_matrix2(*)
    end subroutine rotate_phase_matrix

    subroutine matrix_symmetry(nstokes, matrix1, matrix2)
      integer :: nstokes
      real(8) :: matrix1(*), matrix2(*)
    end subroutine matrix_symmetry

    subroutine fourier_matrix(aziorder, numpts, nstokes, real_matrix, &
        basis_matrix)
      integer :: aziorder, numpts, nstokes
      real(8) :: real_matrix(*), basis_matrix(*)
    end subroutine fourier_matrix

    subroutine combine_phase_modes(nstokes, aziorder, m, tmp, basis_matrix, &
        out_matrix)
      integer :: nstokes, aziorder, m
      real(8) :: tmp
      real(8) :: basis_matrix(*), out_matrix(*)
    end subroutine combine_phase_modes

    subroutine fourier_basis(numbasis, order, numpts, direction, basis_vector, &
        real_vector)
      integer :: numbasis, order, numpts, direction
      real(8) :: basis_vector(*), real_vector(*)
    end subroutine fourier_basis

    subroutine fft1dr(data, n, isign)
      integer :: n, isign
      real(8) :: data(*)
    end subroutine fft1dr

    subroutine fftc(data, n, phase)
      integer :: n
      real(8) :: data(*), phase(*)
    end subroutine fftc

    subroutine fixreal(data, nyquist, n, isign, phase)
      integer :: n, isign
      real(8) :: data(*), nyquist(*), phase(*)
    end subroutine fixreal

    subroutine makephase(phase, nmax)
      integer :: nmax
      real(8) :: phase(*)
    end subroutine makephase

    subroutine scatter_symmetry(nstokes, nummu, scat)
      integer :: nstokes, nummu
      real(8) :: scat(*)
    end subroutine scatter_symmetry
  end interface

  public :: c_rt3_radtran
  public :: c_rt3_double_gauss_quadrature, c_rt3_gauss_legendre_quadrature
  public :: c_rt3_lobatto_quadrature
  public :: c_rt3_get_scat_set, c_rt3_scattering, c_rt3_direct_scattering
  public :: c_rt3_get_scattering, c_rt3_check_norm, c_rt3_get_direct
  public :: c_rt3_thermal_radiance, c_rt3_nonscatter_layer, c_rt3_initial_source
  public :: c_rt3_initialize, c_rt3_doubling_integration, c_rt3_combine_layers
  public :: c_rt3_internal_radiance, c_rt3_lambert_surface, c_rt3_lambert_radiance
  public :: c_rt3_fresnel_surface, c_rt3_fresnel_radiance
  public :: c_rt3_number_sums, c_rt3_sum_legendre, c_rt3_rotate_phase_matrix
  public :: c_rt3_matrix_symmetry, c_rt3_fourier_matrix, c_rt3_combine_phase_modes
  public :: c_rt3_fourier_basis, c_rt3_fft1dr, c_rt3_fftc, c_rt3_fixreal
  public :: c_rt3_makephase, c_rt3_scatter_symmetry

contains

  ! RADTRAN.  scat_nlegen has nsl entries, scatlayers num_layers (1-based
  ! set, 0 for no scattering) and outlevels noutlevels (1 = top,
  ! num_layers + 1 = bottom).  scat_coef is (6, ldcoef, nsl).  mu_values
  ! must hold nummu values; for quad_type 'E' its first entries are 0 and
  ! the extra angles follow, otherwise it is output only.  gas_extinct is
  ! clipped at zero by RADTRAN (the values are not changed).  up_flux and
  ! down_flux hold nstokes*noutlevels values, up_rad and down_rad
  ! nstokes*nummu*(aziorder+1)*noutlevels.
  subroutine c_rt3_radtran(nstokes, nummu, aziorder, max_delta_tau, &
                           src_code, quad_type, deltam, direct_flux, &
                           direct_mu, ground_temp, ground_type, &
                           ground_albedo, ground_index_re, ground_index_im, &
                           sky_temp, wavelength, num_layers, height, &
                           temperatures, gas_extinct, nsl, scat_extinct, &
                           scat_scatter, scat_nlegen, ldcoef, scat_coef, &
                           scatlayers, noutlevels, outlevels, mu_values, &
                           up_flux, down_flux, up_rad, down_rad) &
      bind(C, name="rt3_radtran")
    integer(c_int64_t), value :: nstokes, nummu, aziorder, src_code
    integer(c_int64_t), value :: num_layers, nsl, ldcoef, noutlevels
    integer(c_int64_t), intent(in) :: scat_nlegen(*), scatlayers(*)
    integer(c_int64_t), intent(in) :: outlevels(*)
    real(c_double), value :: max_delta_tau, direct_flux, direct_mu
    real(c_double), value :: ground_temp, ground_albedo
    real(c_double), value :: ground_index_re, ground_index_im
    real(c_double), value :: sky_temp, wavelength
    character(kind=c_char), value :: quad_type, deltam, ground_type
    real(c_double) :: height(*), temperatures(*), gas_extinct(*)
    real(c_double) :: scat_extinct(*), scat_scatter(*), scat_coef(*)
    real(c_double) :: mu_values(*), up_flux(*), down_flux(*)
    real(c_double) :: up_rad(*), down_rad(*)

    integer :: ns, nm, na, sc, nl, ny, ld, no
    integer, allocatable :: nlegen(:), layers(:), levels(:)
    real(8) :: mdt, dflux, dmu, gtemp, galb, stemp, wl
    complex(8) :: gindex
    character(len=1) :: qtype, dm, gtype

    ns = int(nstokes)
    nm = int(nummu)
    na = int(aziorder)
    sc = int(src_code)
    nl = int(num_layers)
    ny = int(nsl)
    ld = int(ldcoef)
    no = int(noutlevels)
    ! At least one element each, so that no zero-size actual argument is
    ! passed for nsl = 0
    allocate(nlegen(max(ny, 1)), layers(max(nl, 1)), levels(max(no, 1)))
    nlegen = 0
    layers = 0
    levels = 0
    nlegen(1:ny) = int(scat_nlegen(1:ny))
    layers(1:nl) = int(scatlayers(1:nl))
    levels(1:no) = int(outlevels(1:no))
    mdt = max_delta_tau
    dflux = direct_flux
    dmu = direct_mu
    gtemp = ground_temp
    galb = ground_albedo
    stemp = sky_temp
    wl = wavelength
    gindex = cmplx(ground_index_re, ground_index_im, kind=8)
    qtype = quad_type
    dm = deltam
    gtype = ground_type

    call radtran(ns, nm, na, mdt, sc, qtype, dm, dflux, dmu, gtemp, &
                 gtype, galb, gindex, stemp, wl, nl, height, temperatures, &
                 gas_extinct, ny, scat_extinct, scat_scatter, nlegen, ld, &
                 scat_coef, layers, no, levels, mu_values, up_flux, &
                 down_flux, up_rad, down_rad)
  end subroutine c_rt3_radtran

  ! The quadratures write num ascending nodes in (0,1] and their weights
  ! for the integral over [0,1].
  subroutine c_rt3_double_gauss_quadrature(num, abscissas, weights) &
      bind(C, name="rt3_double_gauss_quadrature")
    integer(c_int64_t), value :: num
    real(c_double) :: abscissas(*), weights(*)

    integer :: n

    n = int(num)
    call rt3_double_gauss_quadrature(n, abscissas, weights)
  end subroutine c_rt3_double_gauss_quadrature

  subroutine c_rt3_gauss_legendre_quadrature(num, abscissas, weights) &
      bind(C, name="rt3_gauss_legendre_quadrature")
    integer(c_int64_t), value :: num
    real(c_double) :: abscissas(*), weights(*)

    integer :: n

    n = int(num)
    call rt3_gauss_legendre_quadrature(n, abscissas, weights)
  end subroutine c_rt3_gauss_legendre_quadrature

  subroutine c_rt3_lobatto_quadrature(num, abscissas, weights) &
      bind(C, name="rt3_lobatto_quadrature")
    integer(c_int64_t), value :: num
    real(c_double) :: abscissas(*), weights(*)

    integer :: n

    n = int(num)
    call rt3_lobatto_quadrature(n, abscissas, weights)
  end subroutine c_rt3_lobatto_quadrature

  ! GET_SCAT_SET (radscat3.f): nlegen, extinction and scatter are output
  subroutine c_rt3_get_scat_set(deltam, nummu, nlegin, coefin, extin, scatin, &
      nlegen, coef, extinction, scatter) &
      bind(C, name="rt3_get_scat_set")
    integer(c_int64_t), value :: nummu, nlegin
    integer(c_int64_t), intent(out) :: nlegen
    real(c_double), value :: extin, scatin
    real(c_double), intent(out) :: extinction, scatter
    character(kind=c_char), value :: deltam
    real(c_double) :: coefin(*), coef(*)

    integer :: f_nummu, f_nlegin, f_nlegen
    real(8) :: f_extin, f_scatin, f_extinction, f_scatter
    character(len=1) :: f_deltam

    f_deltam = deltam
    f_nummu = int(nummu)
    f_nlegin = int(nlegin)
    f_extin = extin
    f_scatin = scatin
    call get_scat_set(f_deltam, f_nummu, f_nlegin, coefin, f_extin, f_scatin, &
        f_nlegen, coef, f_extinction, f_scatter)
    nlegen = int(f_nlegen, c_int64_t)
    extinction = f_extinction
    scatter = f_scatter
  end subroutine c_rt3_get_scat_set

  ! SCATTERING (radscat3.f)
  subroutine c_rt3_scattering(nummu, aziorder, nstokes, mu_values, &
      quad_weights, numlegendre, legendre_coef, scat_num, scatbuf) &
      bind(C, name="rt3_scattering")
    integer(c_int64_t), value :: nummu, aziorder, nstokes, numlegendre, &
        scat_num
    real(c_double) :: mu_values(*), quad_weights(*), legendre_coef(*), &
        scatbuf(*)

    integer :: f_nummu, f_aziorder, f_nstokes, f_numlegendre, f_scat_num

    f_nummu = int(nummu)
    f_aziorder = int(aziorder)
    f_nstokes = int(nstokes)
    f_numlegendre = int(numlegendre)
    f_scat_num = int(scat_num)
    call scattering(f_nummu, f_aziorder, f_nstokes, mu_values, quad_weights, &
        f_numlegendre, legendre_coef, f_scat_num, scatbuf)
  end subroutine c_rt3_scattering

  ! DIRECT_SCATTERING (radscat3.f)
  subroutine c_rt3_direct_scattering(nummu, aziorder, nstokes, mu_values, &
      numlegendre, legendre_coef, direct_mu, scat_num, directbuf) &
      bind(C, name="rt3_direct_scattering")
    integer(c_int64_t), value :: nummu, aziorder, nstokes, numlegendre, &
        scat_num
    real(c_double), value :: direct_mu
    real(c_double) :: mu_values(*), legendre_coef(*), directbuf(*)

    integer :: f_nummu, f_aziorder, f_nstokes, f_numlegendre, f_scat_num
    real(8) :: f_direct_mu

    f_nummu = int(nummu)
    f_aziorder = int(aziorder)
    f_nstokes = int(nstokes)
    f_numlegendre = int(numlegendre)
    f_direct_mu = direct_mu
    f_scat_num = int(scat_num)
    call direct_scattering(f_nummu, f_aziorder, f_nstokes, mu_values, &
        f_numlegendre, legendre_coef, f_direct_mu, f_scat_num, directbuf)
  end subroutine c_rt3_direct_scattering

  ! GET_SCATTERING (radscat3.f)
  subroutine c_rt3_get_scattering(nstokes, nummu, mode, aziorder, scat_num, &
      scatbuf, scatter_matrix) &
      bind(C, name="rt3_get_scattering")
    integer(c_int64_t), value :: nstokes, nummu, mode, aziorder, scat_num
    real(c_double) :: scatbuf(*), scatter_matrix(*)

    integer :: f_nstokes, f_nummu, f_mode, f_aziorder, f_scat_num

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_mode = int(mode)
    f_aziorder = int(aziorder)
    f_scat_num = int(scat_num)
    call get_scattering(f_nstokes, f_nummu, f_mode, f_aziorder, f_scat_num, &
        scatbuf, scatter_matrix)
  end subroutine c_rt3_get_scattering

  ! RT3_CHECK_NORM (radscat3.f); STOPs if the phase function is not normalised
  subroutine c_rt3_check_norm(nstokes, nummu, quad_weights, scatter_matrix) &
      bind(C, name="rt3_check_norm")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double) :: quad_weights(*), scatter_matrix(*)

    integer :: f_nstokes, f_nummu

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    call rt3_check_norm(f_nstokes, f_nummu, quad_weights, scatter_matrix)
  end subroutine c_rt3_check_norm

  ! GET_DIRECT (radscat3.f)
  subroutine c_rt3_get_direct(nstokes, nummu, mode, aziorder, scat_num, &
      directbuf, direct_vector) &
      bind(C, name="rt3_get_direct")
    integer(c_int64_t), value :: nstokes, nummu, mode, aziorder, scat_num
    real(c_double) :: directbuf(*), direct_vector(*)

    integer :: f_nstokes, f_nummu, f_mode, f_aziorder, f_scat_num

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_mode = int(mode)
    f_aziorder = int(aziorder)
    f_scat_num = int(scat_num)
    call get_direct(f_nstokes, f_nummu, f_mode, f_aziorder, f_scat_num, &
        directbuf, direct_vector)
  end subroutine c_rt3_get_direct

  ! RT3_THERMAL_RADIANCE (radutil3.f)
  subroutine c_rt3_thermal_radiance(nstokes, nummu, mode, temperature, albedo, &
      wavelength, radiance) &
      bind(C, name="rt3_thermal_radiance")
    integer(c_int64_t), value :: nstokes, nummu, mode
    real(c_double), value :: temperature, albedo, wavelength
    real(c_double) :: radiance(*)

    integer :: f_nstokes, f_nummu, f_mode
    real(8) :: f_temperature, f_albedo, f_wavelength

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_mode = int(mode)
    f_temperature = temperature
    f_albedo = albedo
    f_wavelength = wavelength
    call rt3_thermal_radiance(f_nstokes, f_nummu, f_mode, f_temperature, &
        f_albedo, f_wavelength, radiance)
  end subroutine c_rt3_thermal_radiance

  ! RT3_NONSCATTER_LAYER (radintg3.f)
  subroutine c_rt3_nonscatter_layer(nstokes, nummu, mode, deltatau, mu_values, &
      planck0, planck1, reflect, trans, source) &
      bind(C, name="rt3_nonscatter_layer")
    integer(c_int64_t), value :: nstokes, nummu, mode
    real(c_double), value :: deltatau, planck0, planck1
    real(c_double) :: mu_values(*), reflect(*), trans(*), source(*)

    integer :: f_nstokes, f_nummu, f_mode
    real(8) :: f_deltatau, f_planck0, f_planck1

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_mode = int(mode)
    f_deltatau = deltatau
    f_planck0 = planck0
    f_planck1 = planck1
    call rt3_nonscatter_layer(f_nstokes, f_nummu, f_mode, f_deltatau, &
        mu_values, f_planck0, f_planck1, reflect, trans, source)
  end subroutine c_rt3_nonscatter_layer

  ! RT3_INITIAL_SOURCE (radintg3.f)
  subroutine c_rt3_initial_source(nstokes, nummu, n, delta_z, mu_values, &
      extinction, source_vector, source) &
      bind(C, name="rt3_initial_source")
    integer(c_int64_t), value :: nstokes, nummu, n
    real(c_double), value :: delta_z, extinction
    real(c_double) :: mu_values(*), source_vector(*), source(*)

    integer :: f_nstokes, f_nummu, f_n
    real(8) :: f_delta_z, f_extinction

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_n = int(n)
    f_delta_z = delta_z
    f_extinction = extinction
    call rt3_initial_source(f_nstokes, f_nummu, f_n, f_delta_z, mu_values, &
        f_extinction, source_vector, source)
  end subroutine c_rt3_initial_source

  ! RT3_INITIALIZE (radintg3.f)
  subroutine c_rt3_initialize(nstokes, nummu, n, delta_z, mu_values, &
      extinction, albedo, phase_function, reflect, trans) &
      bind(C, name="rt3_initialize")
    integer(c_int64_t), value :: nstokes, nummu, n
    real(c_double), value :: delta_z, extinction, albedo
    real(c_double) :: mu_values(*), phase_function(*), reflect(*), trans(*)

    integer :: f_nstokes, f_nummu, f_n
    real(8) :: f_delta_z, f_extinction, f_albedo

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_n = int(n)
    f_delta_z = delta_z
    f_extinction = extinction
    f_albedo = albedo
    call rt3_initialize(f_nstokes, f_nummu, f_n, f_delta_z, mu_values, &
        f_extinction, f_albedo, phase_function, reflect, trans)
  end subroutine c_rt3_initialize

  ! RT3_DOUBLING_INTEGRATION (radintg3.f); also overwrites reflect, trans,
  ! exp_source and lin_source
  subroutine c_rt3_doubling_integration(n, num_doubles, src_code, symmetric, &
      reflect, trans, exp_source, expfactor, lin_source, linfactor, t_reflect, &
      t_trans, t_source) &
      bind(C, name="rt3_doubling_integration")
    integer(c_int64_t), value :: n, num_doubles, src_code
    real(c_double), value :: expfactor, linfactor
    logical(c_bool), value :: symmetric
    real(c_double) :: reflect(*), trans(*), exp_source(*), lin_source(*), &
        t_reflect(*), t_trans(*), t_source(*)

    integer :: f_n, f_num_doubles, f_src_code
    real(8) :: f_expfactor, f_linfactor
    logical :: f_symmetric

    f_n = int(n)
    f_num_doubles = int(num_doubles)
    f_src_code = int(src_code)
    f_symmetric = symmetric
    f_expfactor = expfactor
    f_linfactor = linfactor
    call rt3_doubling_integration(f_n, f_num_doubles, f_src_code, f_symmetric, &
        reflect, trans, exp_source, f_expfactor, lin_source, f_linfactor, &
        t_reflect, t_trans, t_source)
  end subroutine c_rt3_doubling_integration

  ! RT3_COMBINE_LAYERS (radintg3.f)
  subroutine c_rt3_combine_layers(n, reflect1, trans1, source1, reflect2, &
      trans2, source2, out_reflect, out_trans, out_source) &
      bind(C, name="rt3_combine_layers")
    integer(c_int64_t), value :: n
    real(c_double) :: reflect1(*), trans1(*), source1(*), reflect2(*), &
        trans2(*), source2(*), out_reflect(*), out_trans(*), out_source(*)

    integer :: f_n

    f_n = int(n)
    call rt3_combine_layers(f_n, reflect1, trans1, source1, reflect2, trans2, &
        source2, out_reflect, out_trans, out_source)
  end subroutine c_rt3_combine_layers

  ! RT3_INTERNAL_RADIANCE (radintg3.f)
  subroutine c_rt3_internal_radiance(n, upreflect, uptrans, upsource, &
      downreflect, downtrans, downsource, intoprad, inbottomrad, uprad, &
      downrad) &
      bind(C, name="rt3_internal_radiance")
    integer(c_int64_t), value :: n
    real(c_double) :: upreflect(*), uptrans(*), upsource(*), downreflect(*), &
        downtrans(*), downsource(*), intoprad(*), inbottomrad(*), uprad(*), &
        downrad(*)

    integer :: f_n

    f_n = int(n)
    call rt3_internal_radiance(f_n, upreflect, uptrans, upsource, downreflect, &
        downtrans, downsource, intoprad, inbottomrad, uprad, downrad)
  end subroutine c_rt3_internal_radiance

  ! RT3_LAMBERT_SURFACE (radutil3.f)
  subroutine c_rt3_lambert_surface(nstokes, nummu, mode, mu_values, &
      quad_weights, ground_albedo, reflect, trans, source) &
      bind(C, name="rt3_lambert_surface")
    integer(c_int64_t), value :: nstokes, nummu, mode
    real(c_double), value :: ground_albedo
    real(c_double) :: mu_values(*), quad_weights(*), reflect(*), trans(*), &
        source(*)

    integer :: f_nstokes, f_nummu, f_mode
    real(8) :: f_ground_albedo

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_mode = int(mode)
    f_ground_albedo = ground_albedo
    call rt3_lambert_surface(f_nstokes, f_nummu, f_mode, mu_values, &
        quad_weights, f_ground_albedo, reflect, trans, source)
  end subroutine c_rt3_lambert_surface

  ! RT3_LAMBERT_RADIANCE (radutil3.f)
  subroutine c_rt3_lambert_radiance(nstokes, nummu, mode, src_code, &
      ground_albedo, ground_temp, wavelength, direct_sfc_flux, radiance) &
      bind(C, name="rt3_lambert_radiance")
    integer(c_int64_t), value :: nstokes, nummu, mode, src_code
    real(c_double), value :: ground_albedo, ground_temp, wavelength, &
        direct_sfc_flux
    real(c_double) :: radiance(*)

    integer :: f_nstokes, f_nummu, f_mode, f_src_code
    real(8) :: f_ground_albedo, f_ground_temp, f_wavelength, f_direct_sfc_flux

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_mode = int(mode)
    f_src_code = int(src_code)
    f_ground_albedo = ground_albedo
    f_ground_temp = ground_temp
    f_wavelength = wavelength
    f_direct_sfc_flux = direct_sfc_flux
    call rt3_lambert_radiance(f_nstokes, f_nummu, f_mode, f_src_code, &
        f_ground_albedo, f_ground_temp, f_wavelength, f_direct_sfc_flux, &
        radiance)
  end subroutine c_rt3_lambert_radiance

  ! RT3_FRESNEL_SURFACE (radutil3.f)
  subroutine c_rt3_fresnel_surface(nstokes, nummu, mu_values, index_re, &
      index_im, reflect, trans, source) &
      bind(C, name="rt3_fresnel_surface")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: index_re, index_im
    real(c_double) :: mu_values(*), reflect(*), trans(*), source(*)

    integer :: f_nstokes, f_nummu
    complex(8) :: f_index

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_index = cmplx(index_re, index_im, kind=8)
    call rt3_fresnel_surface(f_nstokes, f_nummu, mu_values, f_index, reflect, &
        trans, source)
  end subroutine c_rt3_fresnel_surface

  ! RT3_FRESNEL_RADIANCE (radutil3.f)
  subroutine c_rt3_fresnel_radiance(nstokes, nummu, mode, mu_values, index_re, &
      index_im, ground_temp, wavelength, radiance) &
      bind(C, name="rt3_fresnel_radiance")
    integer(c_int64_t), value :: nstokes, nummu, mode
    real(c_double), value :: ground_temp, wavelength, index_re, index_im
    real(c_double) :: mu_values(*), radiance(*)

    integer :: f_nstokes, f_nummu, f_mode
    real(8) :: f_ground_temp, f_wavelength
    complex(8) :: f_index

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    f_mode = int(mode)
    f_index = cmplx(index_re, index_im, kind=8)
    f_ground_temp = ground_temp
    f_wavelength = wavelength
    call rt3_fresnel_radiance(f_nstokes, f_nummu, f_mode, mu_values, f_index, &
        f_ground_temp, f_wavelength, radiance)
  end subroutine c_rt3_fresnel_radiance

  ! NUMBER_SUMS (radscat3.f).  dosum(6) is output.
  subroutine c_rt3_number_sums(nstokes, nlegen, coef, dosum) &
      bind(C, name="rt3_number_sums")
    integer(c_int64_t), value :: nstokes, nlegen
    real(c_double) :: coef(*)
    integer(c_int64_t) :: dosum(6)

    integer :: f_nstokes, f_nlegen, f_dosum(6)

    f_nstokes = int(nstokes)
    f_nlegen = int(nlegen)
    call number_sums(f_nstokes, f_nlegen, coef, f_dosum)
    dosum = int(f_dosum, c_int64_t)
  end subroutine c_rt3_number_sums

  ! SUM_LEGENDRE (radscat3.f).  dosum(6) is input.
  subroutine c_rt3_sum_legendre(nlegen, coef, x, dosum, phase_matrix) &
      bind(C, name="rt3_sum_legendre")
    integer(c_int64_t), value :: nlegen
    real(c_double), value :: x
    real(c_double) :: coef(*), phase_matrix(*)
    integer(c_int64_t) :: dosum(6)

    integer :: f_nlegen, f_dosum(6)
    real(8) :: f_x

    f_nlegen = int(nlegen)
    f_dosum = int(dosum)
    f_x = x
    call sum_legendre(f_nlegen, coef, f_x, f_dosum, phase_matrix)
  end subroutine c_rt3_sum_legendre

  ! ROTATE_PHASE_MATRIX (radscat3.f)
  subroutine c_rt3_rotate_phase_matrix(phase_matrix1, mu1, mu2, delphi, &
      cos_scat, phase_matrix2, nstokes) &
      bind(C, name="rt3_rotate_phase_matrix")
    real(c_double), value :: mu1, mu2, delphi, cos_scat
    integer(c_int64_t), value :: nstokes
    real(c_double) :: phase_matrix1(*), phase_matrix2(*)

    integer :: f_nstokes
    real(8) :: f_mu1, f_mu2, f_delphi, f_cos_scat

    f_nstokes = int(nstokes)
    f_mu1 = mu1
    f_mu2 = mu2
    f_delphi = delphi
    f_cos_scat = cos_scat
    call rotate_phase_matrix(phase_matrix1, f_mu1, f_mu2, f_delphi, &
        f_cos_scat, phase_matrix2, f_nstokes)
  end subroutine c_rt3_rotate_phase_matrix

  ! MATRIX_SYMMETRY (radscat3.f).  SCATTERING calls it once with matrix1
  ! and matrix2 the same array.
  subroutine c_rt3_matrix_symmetry(nstokes, matrix1, matrix2) &
      bind(C, name="rt3_matrix_symmetry")
    integer(c_int64_t), value :: nstokes
    real(c_double) :: matrix1(*), matrix2(*)

    integer :: f_nstokes

    f_nstokes = int(nstokes)
    call matrix_symmetry(f_nstokes, matrix1, matrix2)
  end subroutine c_rt3_matrix_symmetry

  ! FOURIER_MATRIX (radscat3.f)
  subroutine c_rt3_fourier_matrix(aziorder, numpts, nstokes, real_matrix, &
      basis_matrix) &
      bind(C, name="rt3_fourier_matrix")
    integer(c_int64_t), value :: aziorder, numpts, nstokes
    real(c_double) :: real_matrix(*), basis_matrix(*)

    integer :: f_aziorder, f_numpts, f_nstokes

    f_aziorder = int(aziorder)
    f_numpts = int(numpts)
    f_nstokes = int(nstokes)
    call fourier_matrix(f_aziorder, f_numpts, f_nstokes, real_matrix, &
        basis_matrix)
  end subroutine c_rt3_fourier_matrix

  ! COMBINE_PHASE_MODES (radscat3.f)
  subroutine c_rt3_combine_phase_modes(nstokes, aziorder, m, tmp, &
      basis_matrix, out_matrix) &
      bind(C, name="rt3_combine_phase_modes")
    integer(c_int64_t), value :: nstokes, aziorder, m
    real(c_double), value :: tmp
    real(c_double) :: basis_matrix(*), out_matrix(*)

    integer :: f_nstokes, f_aziorder, f_m
    real(8) :: f_tmp

    f_nstokes = int(nstokes)
    f_aziorder = int(aziorder)
    f_m = int(m)
    f_tmp = tmp
    call combine_phase_modes(f_nstokes, f_aziorder, f_m, f_tmp, &
        basis_matrix, out_matrix)
  end subroutine c_rt3_combine_phase_modes

  ! FOURIER_BASIS (radscat3.f).  real_vector is overwritten (FFT1DR works
  ! in place) when order > 0.
  subroutine c_rt3_fourier_basis(numbasis, order, numpts, direction, &
      basis_vector, real_vector) &
      bind(C, name="rt3_fourier_basis")
    integer(c_int64_t), value :: numbasis, order, numpts, direction
    real(c_double) :: basis_vector(*), real_vector(*)

    integer :: f_numbasis, f_order, f_numpts, f_direction

    f_numbasis = int(numbasis)
    f_order = int(order)
    f_numpts = int(numpts)
    f_direction = int(direction)
    call fourier_basis(f_numbasis, f_order, f_numpts, f_direction, &
        basis_vector, real_vector)
  end subroutine c_rt3_fourier_basis

  ! FFT1DR (radscat3.f), in place.  It keeps its phase table in SAVEd
  ! variables.
  subroutine c_rt3_fft1dr(data, n, isign) bind(C, name="rt3_fft1dr")
    integer(c_int64_t), value :: n, isign
    real(c_double) :: data(*)

    integer :: f_n, f_isign

    f_n = int(n)
    f_isign = int(isign)
    call fft1dr(data, f_n, f_isign)
  end subroutine c_rt3_fft1dr

  ! FFTC (radscat3.f): complex FFT of n points in place, with the phase
  ! table of MAKEPHASE (its + or - half)
  subroutine c_rt3_fftc(data, n, phase) bind(C, name="rt3_fftc")
    integer(c_int64_t), value :: n
    real(c_double) :: data(*), phase(*)

    integer :: f_n

    f_n = int(n)
    call fftc(data, f_n, phase)
  end subroutine c_rt3_fftc

  ! FIXREAL (radscat3.f).  nyquist(2) is output for isign > 0, nyquist(1)
  ! input otherwise.
  subroutine c_rt3_fixreal(data, nyquist, n, isign, phase) &
      bind(C, name="rt3_fixreal")
    integer(c_int64_t), value :: n, isign
    real(c_double) :: data(*), nyquist(*), phase(*)

    integer :: f_n, f_isign

    f_n = int(n)
    f_isign = int(isign)
    call fixreal(data, nyquist, f_n, f_isign, phase)
  end subroutine c_rt3_fixreal

  ! MAKEPHASE (radscat3.f): phase holds 4*nmax values
  subroutine c_rt3_makephase(phase, nmax) bind(C, name="rt3_makephase")
    integer(c_int64_t), value :: nmax
    real(c_double) :: phase(*)

    integer :: f_nmax

    f_nmax = int(nmax)
    call makephase(phase, f_nmax)
  end subroutine c_rt3_makephase

  ! SCATTER_SYMMETRY (radscat3.f): scat is (nstokes, nummu, nstokes, nummu,
  ! 4); parts 3 and 4 are made from 2 and 1
  subroutine c_rt3_scatter_symmetry(nstokes, nummu, scat) &
      bind(C, name="rt3_scatter_symmetry")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double) :: scat(*)

    integer :: f_nstokes, f_nummu

    f_nstokes = int(nstokes)
    f_nummu = int(nummu)
    call scatter_symmetry(f_nstokes, f_nummu, scat)
  end subroutine c_rt3_scatter_symmetry
end module rt3_c_interface
