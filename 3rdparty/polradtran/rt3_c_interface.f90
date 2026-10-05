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
! The C names (rt3_radtran, rt3_*_quadrature) differ from the linker
! names of the renamed Fortran routines (rt3_*_quadrature_, with the
! trailing underscore of gfortran's external names).
!
! The real arrays are passed through unchanged.  Their layouts are
! documented at the RADTRAN header in radtran3.f and in
! src/core/rt3/rt3.h.  None of these routines is reentrant: RT3 uses
! COMMON blocks, SAVEd FFT tables and large static local arrays, so the
! caller must serialise all calls.
module rt3_c_interface
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_double, c_char
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
  end interface

  public :: c_rt3_radtran
  public :: c_rt3_double_gauss_quadrature, c_rt3_gauss_legendre_quadrature
  public :: c_rt3_lobatto_quadrature

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
end module rt3_c_interface
