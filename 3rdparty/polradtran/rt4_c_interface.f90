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
! src/core/polradtran/rt4/rt4.h.  None of these routines is reentrant: RT4 uses
! COMMON blocks and large static local arrays, so the caller must
! serialise all calls.
!
! Besides RADTRANO itself, every subroutine that RADTRANO calls has an
! entry point here (rt4_ and the lower-case name), for the C++ port of
! RADTRANO in src/core/polradtran/rt4/radtran4.cc.
module rt4_c_interface
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_double, c_char, c_bool
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

    ! radintg4.f
    subroutine initialize(nstokes, nummu, delta_z, mu_values, &
                          quad_weights, gas_extinct, extinct_matrix, &
                          scatter_matrix, reflect, trans)
      integer :: nstokes, nummu
      real(8) :: delta_z, gas_extinct
      real(8) :: mu_values(*), quad_weights(*), extinct_matrix(*)
      real(8) :: scatter_matrix(*), reflect(*), trans(*)
    end subroutine initialize

    subroutine initial_source(nstokes, nummu, delta_z, mu_values, &
                              planck, emis_vector, gas_extinct, source)
      integer :: nstokes, nummu
      real(8) :: delta_z, planck, gas_extinct
      real(8) :: mu_values(*), emis_vector(*), source(*)
    end subroutine initial_source

    subroutine nonscatter_layer(nstokes, nummu, deltatau, mu_values, &
                                planck0, planck1, reflect, trans, source)
      integer :: nstokes, nummu
      real(8) :: deltatau, planck0, planck1
      real(8) :: mu_values(*), reflect(*), trans(*), source(*)
    end subroutine nonscatter_layer

    subroutine internal_radiance(n, upreflect, uptrans, upsource, &
                                 downreflect, downtrans, downsource, &
                                 intoprad, inbottomrad, uprad, downrad)
      integer :: n
      real(8) :: upreflect(*), uptrans(*), upsource(*)
      real(8) :: downreflect(*), downtrans(*), downsource(*)
      real(8) :: intoprad(*), inbottomrad(*), uprad(*), downrad(*)
    end subroutine internal_radiance

    subroutine doubling_integration(n, num_doubles, symmetric, reflect, &
                                    trans, lin_source, linfactor, &
                                    t_reflect, t_trans, t_source)
      integer :: n, num_doubles
      logical :: symmetric
      real(8) :: linfactor
      real(8) :: reflect(*), trans(*), lin_source(*)
      real(8) :: t_reflect(*), t_trans(*), t_source(*)
    end subroutine doubling_integration

    subroutine combine_layers(n, reflect1, trans1, source1, reflect2, &
                              trans2, source2, out_reflect, out_trans, &
                              out_source)
      integer :: n
      real(8) :: reflect1(*), trans1(*), source1(*)
      real(8) :: reflect2(*), trans2(*), source2(*)
      real(8) :: out_reflect(*), out_trans(*), out_source(*)
    end subroutine combine_layers

    ! radutil4.f
    subroutine lambert_surface(nstokes, nummu, mode, mu_values, &
                               quad_weights, ground_albedo, reflect, &
                               trans, source)
      integer :: nstokes, nummu, mode
      real(8) :: ground_albedo
      real(8) :: mu_values(*), quad_weights(*)
      real(8) :: reflect(*), trans(*), source(*)
    end subroutine lambert_surface

    subroutine lambert_radiance(nstokes, nummu, ground_albedo, &
                                ground_temp, wavelength, radiance)
      integer :: nstokes, nummu
      real(8) :: ground_albedo, ground_temp, wavelength
      real(8) :: radiance(*)
    end subroutine lambert_radiance

    subroutine fresnel_surface(nstokes, nummu, mu_values, index, &
                               reflect, trans, source)
      integer :: nstokes, nummu
      complex(8) :: index
      real(8) :: mu_values(*), reflect(*), trans(*), source(*)
    end subroutine fresnel_surface

    subroutine fresnel_radiance(nstokes, nummu, mu_values, index, &
                                ground_temp, wavelength, radiance)
      integer :: nstokes, nummu
      complex(8) :: index
      real(8) :: ground_temp, wavelength
      real(8) :: mu_values(*), radiance(*)
    end subroutine fresnel_radiance

    subroutine specular_surface(nstokes, nummu, ground_reflec, reflect, &
                                trans, source)
      integer :: nstokes, nummu
      real(8) :: ground_reflec(*), reflect(*), trans(*), source(*)
    end subroutine specular_surface

    subroutine specular_radiance(nstokes, nummu, ground_reflec, &
                                 ground_temp, wavelength, radiance)
      integer :: nstokes, nummu
      real(8) :: ground_temp, wavelength
      real(8) :: ground_reflec(*), radiance(*)
    end subroutine specular_radiance

    subroutine external_surface(nstokes, nummu, surf_refl, radiance, &
                                reflect, trans, source)
      integer :: nstokes, nummu
      real(8) :: surf_refl(*), radiance(*)
      real(8) :: reflect(*), trans(*), source(*)
    end subroutine external_surface

    subroutine thermal_radiance(nstokes, nummu, temperature, albedo, &
                                wavelength, radiance)
      integer :: nstokes, nummu
      real(8) :: temperature, albedo, wavelength
      real(8) :: radiance(*)
    end subroutine thermal_radiance

    ! radmat.f
    subroutine mcopy(n, m, matrix1, matrix2)
      integer :: n, m
      real(8) :: matrix1(*), matrix2(*)
    end subroutine mcopy

    subroutine mzero(n, m, matrix1)
      integer :: n, m
      real(8) :: matrix1(*)
    end subroutine mzero

    subroutine midentity(n, matrix)
      integer :: n
      real(8) :: matrix(*)
    end subroutine midentity
  end interface

  public :: rt4_radtrano, rt4_planck_function
  public :: rt4_double_gauss_quadrature, rt4_gauss_legendre_quadrature
  public :: rt4_lobatto_quadrature
  public :: rt4_initialize, rt4_initial_source, rt4_nonscatter_layer
  public :: rt4_internal_radiance, rt4_doubling_integration
  public :: rt4_combine_layers
  public :: rt4_lambert_surface, rt4_lambert_radiance
  public :: rt4_fresnel_surface, rt4_fresnel_radiance
  public :: rt4_specular_surface, rt4_specular_radiance
  public :: rt4_external_surface, rt4_thermal_radiance
  public :: rt4_mcopy, rt4_mzero, rt4_midentity

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

  ! INITIALIZE (radintg4.f)
  subroutine rt4_initialize(nstokes, nummu, delta_z, mu_values, &
                            quad_weights, gas_extinct, extinct_matrix, &
                            scatter_matrix, reflect, trans) &
      bind(C, name="rt4_initialize")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: delta_z, gas_extinct
    real(c_double) :: mu_values(*), quad_weights(*), extinct_matrix(*)
    real(c_double) :: scatter_matrix(*), reflect(*), trans(*)

    integer :: ns, nm
    real(8) :: dz, gext

    ns = int(nstokes)
    nm = int(nummu)
    dz = delta_z
    gext = gas_extinct
    call initialize(ns, nm, dz, mu_values, quad_weights, gext, &
                    extinct_matrix, scatter_matrix, reflect, trans)
  end subroutine rt4_initialize

  ! INITIAL_SOURCE (radintg4.f)
  subroutine rt4_initial_source(nstokes, nummu, delta_z, mu_values, &
                                planck, emis_vector, gas_extinct, source) &
      bind(C, name="rt4_initial_source")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: delta_z, planck, gas_extinct
    real(c_double) :: mu_values(*), emis_vector(*), source(*)

    integer :: ns, nm
    real(8) :: dz, b, gext

    ns = int(nstokes)
    nm = int(nummu)
    dz = delta_z
    b = planck
    gext = gas_extinct
    call initial_source(ns, nm, dz, mu_values, b, emis_vector, gext, &
                        source)
  end subroutine rt4_initial_source

  ! NONSCATTER_LAYER (radintg4.f)
  subroutine rt4_nonscatter_layer(nstokes, nummu, deltatau, mu_values, &
                                  planck0, planck1, reflect, trans, &
                                  source) &
      bind(C, name="rt4_nonscatter_layer")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: deltatau, planck0, planck1
    real(c_double) :: mu_values(*), reflect(*), trans(*), source(*)

    integer :: ns, nm
    real(8) :: dtau, b0, b1

    ns = int(nstokes)
    nm = int(nummu)
    dtau = deltatau
    b0 = planck0
    b1 = planck1
    call nonscatter_layer(ns, nm, dtau, mu_values, b0, b1, reflect, &
                          trans, source)
  end subroutine rt4_nonscatter_layer

  ! INTERNAL_RADIANCE (radintg4.f)
  subroutine rt4_internal_radiance(n, upreflect, uptrans, upsource, &
                                   downreflect, downtrans, downsource, &
                                   intoprad, inbottomrad, uprad, &
                                   downrad) &
      bind(C, name="rt4_internal_radiance")
    integer(c_int64_t), value :: n
    real(c_double) :: upreflect(*), uptrans(*), upsource(*)
    real(c_double) :: downreflect(*), downtrans(*), downsource(*)
    real(c_double) :: intoprad(*), inbottomrad(*), uprad(*), downrad(*)

    integer :: nn

    nn = int(n)
    call internal_radiance(nn, upreflect, uptrans, upsource, &
                           downreflect, downtrans, downsource, &
                           intoprad, inbottomrad, uprad, downrad)
  end subroutine rt4_internal_radiance

  ! DOUBLING_INTEGRATION (radintg4.f)
  subroutine rt4_doubling_integration(n, num_doubles, symmetric, &
                                      reflect, trans, lin_source, &
                                      linfactor, t_reflect, t_trans, &
                                      t_source) &
      bind(C, name="rt4_doubling_integration")
    integer(c_int64_t), value :: n, num_doubles
    logical(c_bool), value :: symmetric
    real(c_double), value :: linfactor
    real(c_double) :: reflect(*), trans(*), lin_source(*)
    real(c_double) :: t_reflect(*), t_trans(*), t_source(*)

    integer :: nn, nd
    logical :: sym
    real(8) :: linf

    nn = int(n)
    nd = int(num_doubles)
    sym = symmetric
    linf = linfactor
    call doubling_integration(nn, nd, sym, reflect, trans, lin_source, &
                              linf, t_reflect, t_trans, t_source)
  end subroutine rt4_doubling_integration

  ! COMBINE_LAYERS (radintg4.f)
  subroutine rt4_combine_layers(n, reflect1, trans1, source1, reflect2, &
                                trans2, source2, out_reflect, out_trans, &
                                out_source) &
      bind(C, name="rt4_combine_layers")
    integer(c_int64_t), value :: n
    real(c_double) :: reflect1(*), trans1(*), source1(*)
    real(c_double) :: reflect2(*), trans2(*), source2(*)
    real(c_double) :: out_reflect(*), out_trans(*), out_source(*)

    integer :: nn

    nn = int(n)
    call combine_layers(nn, reflect1, trans1, source1, reflect2, &
                        trans2, source2, out_reflect, out_trans, &
                        out_source)
  end subroutine rt4_combine_layers

  ! LAMBERT_SURFACE (radutil4.f)
  subroutine rt4_lambert_surface(nstokes, nummu, mode, mu_values, &
                                 quad_weights, ground_albedo, reflect, &
                                 trans, source) &
      bind(C, name="rt4_lambert_surface")
    integer(c_int64_t), value :: nstokes, nummu, mode
    real(c_double), value :: ground_albedo
    real(c_double) :: mu_values(*), quad_weights(*)
    real(c_double) :: reflect(*), trans(*), source(*)

    integer :: ns, nm, md
    real(8) :: alb

    ns = int(nstokes)
    nm = int(nummu)
    md = int(mode)
    alb = ground_albedo
    call lambert_surface(ns, nm, md, mu_values, quad_weights, alb, &
                         reflect, trans, source)
  end subroutine rt4_lambert_surface

  ! LAMBERT_RADIANCE (radutil4.f)
  subroutine rt4_lambert_radiance(nstokes, nummu, ground_albedo, &
                                  ground_temp, wavelength, radiance) &
      bind(C, name="rt4_lambert_radiance")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: ground_albedo, ground_temp, wavelength
    real(c_double) :: radiance(*)

    integer :: ns, nm
    real(8) :: alb, gtemp, wl

    ns = int(nstokes)
    nm = int(nummu)
    alb = ground_albedo
    gtemp = ground_temp
    wl = wavelength
    call lambert_radiance(ns, nm, alb, gtemp, wl, radiance)
  end subroutine rt4_lambert_radiance

  ! FRESNEL_SURFACE (radutil4.f)
  subroutine rt4_fresnel_surface(nstokes, nummu, mu_values, index_re, &
                                 index_im, reflect, trans, source) &
      bind(C, name="rt4_fresnel_surface")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: index_re, index_im
    real(c_double) :: mu_values(*), reflect(*), trans(*), source(*)

    integer :: ns, nm
    complex(8) :: gindex

    ns = int(nstokes)
    nm = int(nummu)
    gindex = cmplx(index_re, index_im, kind=8)
    call fresnel_surface(ns, nm, mu_values, gindex, reflect, trans, &
                         source)
  end subroutine rt4_fresnel_surface

  ! FRESNEL_RADIANCE (radutil4.f)
  subroutine rt4_fresnel_radiance(nstokes, nummu, mu_values, index_re, &
                                  index_im, ground_temp, wavelength, &
                                  radiance) &
      bind(C, name="rt4_fresnel_radiance")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: index_re, index_im, ground_temp, wavelength
    real(c_double) :: mu_values(*), radiance(*)

    integer :: ns, nm
    complex(8) :: gindex
    real(8) :: gtemp, wl

    ns = int(nstokes)
    nm = int(nummu)
    gindex = cmplx(index_re, index_im, kind=8)
    gtemp = ground_temp
    wl = wavelength
    call fresnel_radiance(ns, nm, mu_values, gindex, gtemp, wl, radiance)
  end subroutine rt4_fresnel_radiance

  ! SPECULAR_SURFACE (radutil4.f)
  subroutine rt4_specular_surface(nstokes, nummu, ground_reflec, &
                                  reflect, trans, source) &
      bind(C, name="rt4_specular_surface")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double) :: ground_reflec(*), reflect(*), trans(*), source(*)

    integer :: ns, nm

    ns = int(nstokes)
    nm = int(nummu)
    call specular_surface(ns, nm, ground_reflec, reflect, trans, source)
  end subroutine rt4_specular_surface

  ! SPECULAR_RADIANCE (radutil4.f)
  subroutine rt4_specular_radiance(nstokes, nummu, ground_reflec, &
                                   ground_temp, wavelength, radiance) &
      bind(C, name="rt4_specular_radiance")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: ground_temp, wavelength
    real(c_double) :: ground_reflec(*), radiance(*)

    integer :: ns, nm
    real(8) :: gtemp, wl

    ns = int(nstokes)
    nm = int(nummu)
    gtemp = ground_temp
    wl = wavelength
    call specular_radiance(ns, nm, ground_reflec, gtemp, wl, radiance)
  end subroutine rt4_specular_radiance

  ! EXTERNAL_SURFACE (radutil4.f)
  subroutine rt4_external_surface(nstokes, nummu, surf_refl, radiance, &
                                  reflect, trans, source) &
      bind(C, name="rt4_external_surface")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double) :: surf_refl(*), radiance(*)
    real(c_double) :: reflect(*), trans(*), source(*)

    integer :: ns, nm

    ns = int(nstokes)
    nm = int(nummu)
    call external_surface(ns, nm, surf_refl, radiance, reflect, trans, &
                          source)
  end subroutine rt4_external_surface

  ! THERMAL_RADIANCE (radutil4.f)
  subroutine rt4_thermal_radiance(nstokes, nummu, temperature, albedo, &
                                  wavelength, radiance) &
      bind(C, name="rt4_thermal_radiance")
    integer(c_int64_t), value :: nstokes, nummu
    real(c_double), value :: temperature, albedo, wavelength
    real(c_double) :: radiance(*)

    integer :: ns, nm
    real(8) :: t, alb, wl

    ns = int(nstokes)
    nm = int(nummu)
    t = temperature
    alb = albedo
    wl = wavelength
    call thermal_radiance(ns, nm, t, alb, wl, radiance)
  end subroutine rt4_thermal_radiance

  ! MCOPY (radmat.f)
  subroutine rt4_mcopy(n, m, matrix1, matrix2) bind(C, name="rt4_mcopy")
    integer(c_int64_t), value :: n, m
    real(c_double) :: matrix1(*), matrix2(*)

    integer :: nn, mm

    nn = int(n)
    mm = int(m)
    call mcopy(nn, mm, matrix1, matrix2)
  end subroutine rt4_mcopy

  ! MZERO (radmat.f)
  subroutine rt4_mzero(n, m, matrix1) bind(C, name="rt4_mzero")
    integer(c_int64_t), value :: n, m
    real(c_double) :: matrix1(*)

    integer :: nn, mm

    nn = int(n)
    mm = int(m)
    call mzero(nn, mm, matrix1)
  end subroutine rt4_mzero

  ! MIDENTITY (radmat.f)
  subroutine rt4_midentity(n, matrix) bind(C, name="rt4_midentity")
    integer(c_int64_t), value :: n
    real(c_double) :: matrix(*)

    integer :: nn

    nn = int(n)
    call midentity(nn, matrix)
  end subroutine rt4_midentity
end module rt4_c_interface
