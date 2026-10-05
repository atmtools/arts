! ARTS3: ISO_C_BINDING entry points for the RT4 library.
!
! The Fortran 77 sources are compiled with default (4-byte) INTEGER and
! with REAL*8 / COMPLEX*16 arguments.  These wrappers take explicit C
! types (int64_t scalars and the one-character flags by value, doubles,
! the complex ground index as two doubles), convert them to the default
! Fortran kinds and call the original routines.  This removes the need
! for -fdefault-integer-8 and for the compiler-specific hidden string
! length arguments of CHARACTER dummies.
!
! The array arguments are passed through unchanged.  Their layouts are
! documented at the RADTRANO header in radtran4.f and in
! src/core/rt4/rt4.h.  None of these routines is reentrant: RT4 uses
! COMMON blocks and large static local arrays, so the caller must
! serialise all calls.
module rt4_c_interface
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_double, c_char
  implicit none
  private

  interface
    subroutine radtrano(nstokes, nummu, nuummu, max_delta_tau, &
                        quad_type, ground_temp, ground_type, &
                        ground_albedo, ground_index, ground_reflec, &
                        surf_reflect, gnd_radiance, sky_temp, wavelength, &
                        num_layers, height, temperatures, gas_extinct, &
                        nsl, scatlayers, extinct_matrix, emis_vector, &
                        scatter_matrix, mu_values, up_rad, down_rad)
      integer :: nstokes, nummu, nuummu, num_layers, nsl
      real(8) :: max_delta_tau, ground_temp, ground_albedo
      real(8) :: sky_temp, wavelength
      complex(8) :: ground_index
      character(len=1) :: quad_type, ground_type
      real(8) :: ground_reflec(*), surf_reflect(*), gnd_radiance(*)
      real(8) :: height(*), temperatures(*), gas_extinct(*)
      real(8) :: scatlayers(*), extinct_matrix(*), emis_vector(*)
      real(8) :: scatter_matrix(*), mu_values(*), up_rad(*), down_rad(*)
    end subroutine radtrano

    subroutine planck_function(temp, units, wavelength, planck)
      real(8) :: temp, wavelength, planck
      character(len=1) :: units
    end subroutine planck_function

    subroutine double_gauss_quadrature(num, abscissas, weights)
      integer :: num
      real(8) :: abscissas(*), weights(*)
    end subroutine double_gauss_quadrature

    subroutine gauss_legendre_quadrature(num, abscissas, weights)
      integer :: num
      real(8) :: abscissas(*), weights(*)
    end subroutine gauss_legendre_quadrature

    subroutine lobatto_quadrature(num, abscissas, weights)
      integer :: num
      real(8) :: abscissas(*), weights(*)
    end subroutine lobatto_quadrature
  end interface

  public :: rt4_radtrano, rt4_planck_function
  public :: rt4_double_gauss_quadrature, rt4_gauss_legendre_quadrature
  public :: rt4_lobatto_quadrature

contains

  ! RADTRANO.  gnd_radiance is input for ground_type 'A' and is overwritten
  ! for 'F', 'L' and 'S'; gas_extinct is clipped at zero in place;
  ! mu_values returns the quadrature nodes in its first nummu-nuummu
  ! entries and must hold the extra angles in the remaining ones.
  ! up_rad and down_rad must hold nstokes*nummu*(num_layers+1) values.
  subroutine rt4_radtrano(nstokes, nummu, nuummu, max_delta_tau, &
                          quad_type, ground_temp, ground_type, &
                          ground_albedo, ground_index_re, ground_index_im, &
                          ground_reflec, surf_reflect, gnd_radiance, &
                          sky_temp, wavelength, num_layers, height, &
                          temperatures, gas_extinct, nsl, scatlayers, &
                          extinct_matrix, emis_vector, scatter_matrix, &
                          mu_values, up_rad, down_rad) &
      bind(C, name="rt4_radtrano")
    integer(c_int64_t), value :: nstokes, nummu, nuummu, num_layers, nsl
    real(c_double), value :: max_delta_tau, ground_temp, ground_albedo
    real(c_double), value :: ground_index_re, ground_index_im
    real(c_double), value :: sky_temp, wavelength
    character(kind=c_char), value :: quad_type, ground_type
    real(c_double) :: ground_reflec(*), surf_reflect(*), gnd_radiance(*)
    real(c_double) :: height(*), temperatures(*), gas_extinct(*)
    real(c_double) :: scatlayers(*), extinct_matrix(*), emis_vector(*)
    real(c_double) :: scatter_matrix(*), mu_values(*), up_rad(*), down_rad(*)

    integer :: ns, nm, nu, nl, ny
    real(8) :: mdt, gtemp, galb, stemp, wl
    complex(8) :: gindex
    character(len=1) :: qtype, gtype

    ns = int(nstokes)
    nm = int(nummu)
    nu = int(nuummu)
    nl = int(num_layers)
    ny = int(nsl)
    mdt = max_delta_tau
    gtemp = ground_temp
    galb = ground_albedo
    stemp = sky_temp
    wl = wavelength
    gindex = cmplx(ground_index_re, ground_index_im, kind=8)
    qtype = quad_type
    gtype = ground_type

    call radtrano(ns, nm, nu, mdt, qtype, gtemp, gtype, galb, gindex, &
                  ground_reflec, surf_reflect, gnd_radiance, stemp, wl, &
                  nl, height, temperatures, gas_extinct, ny, scatlayers, &
                  extinct_matrix, emis_vector, scatter_matrix, mu_values, &
                  up_rad, down_rad)
  end subroutine rt4_radtrano

  ! PLANCK_FUNCTION with radiance units ('R'): W m-2 sr-1 um-1 for a
  ! temperature in K and a wavelength in um; 0 for temp <= 0.
  function rt4_planck_function(temp, wavelength) result(planck) &
      bind(C, name="rt4_planck_function")
    real(c_double), value :: temp, wavelength
    real(c_double) :: planck

    real(8) :: t, wl, b

    t = temp
    wl = wavelength
    call planck_function(t, 'R', wl, b)
    planck = b
  end function rt4_planck_function

  ! The quadratures write num ascending nodes in (0,1] and their weights
  ! for the integral over [0,1].
  subroutine rt4_double_gauss_quadrature(num, abscissas, weights) &
      bind(C, name="rt4_double_gauss_quadrature")
    integer(c_int64_t), value :: num
    real(c_double) :: abscissas(*), weights(*)

    integer :: n

    n = int(num)
    call double_gauss_quadrature(n, abscissas, weights)
  end subroutine rt4_double_gauss_quadrature

  subroutine rt4_gauss_legendre_quadrature(num, abscissas, weights) &
      bind(C, name="rt4_gauss_legendre_quadrature")
    integer(c_int64_t), value :: num
    real(c_double) :: abscissas(*), weights(*)

    integer :: n

    n = int(num)
    call gauss_legendre_quadrature(n, abscissas, weights)
  end subroutine rt4_gauss_legendre_quadrature

  subroutine rt4_lobatto_quadrature(num, abscissas, weights) &
      bind(C, name="rt4_lobatto_quadrature")
    integer(c_int64_t), value :: num
    real(c_double) :: abscissas(*), weights(*)

    integer :: n

    n = int(num)
    call lobatto_quadrature(n, abscissas, weights)
  end subroutine rt4_lobatto_quadrature
end module rt4_c_interface
