! Fortran 2003 interface to the SHTns library (with contributions from T. Gastine).

  integer(C_INT), parameter :: SHT_NATIVE_LAYOUT=0
  integer(C_INT), parameter :: SHT_THETA_CONTIGUOUS=256
  integer(C_INT), parameter :: SHT_PHI_CONTIGUOUS=512
  integer(C_INT), parameter :: SHT_NO_CS_PHASE=1024
  integer(C_INT), parameter :: SHT_REAL_NORM=2048
  integer(C_INT), parameter :: SHT_SOUTH_POLE_FIRST=8192
  integer(C_INT), parameter :: SHT_SCALAR_ONLY=4096
  integer(C_INT), parameter :: SHT_LOAD_SAVE_CFG=16384
  integer(C_INT), parameter :: SHT_ALLOW_GPU=32768
  integer(C_INT), parameter :: SHT_ALLOW_PADDING=65536

  enum, bind(C)
    enumerator :: SHT_GAUSS=0, SHT_AUTO=1, SHT_REG_FAST=2
    enumerator :: SHT_REG_DCT=3, SHT_QUICK_INIT=4, SHT_REG_POLES=5
    enumerator :: SHT_GAUSS_FLY=6
  end enum

  enum, bind(C)
    enumerator :: SHT_ORTHONORMAL=0,SHT_FOURPI=1,SHT_SCHMIDT=2
  end enum

  type, bind(C) :: shtns_info
     integer(C_INT) :: nlm
     integer(C_SHORT) :: lmax, mmax, mres
     integer(C_SHORT) :: nlat_2
     integer(C_INT) :: nlat, nphi, nspat
     type(C_PTR) :: li, mi
     type(C_PTR) :: ct, st
     integer(C_INT) :: nlat_padded, nlm_cplx
  end type shtns_info

  interface

    subroutine shtns_verbose(i) bind(C, name='shtns_verbose')
      import
      integer(C_INT), value :: i
    end subroutine shtns_verbose
    
    subroutine shtns_print_version() bind(C, name='shtns_print_version')
      import
    end subroutine shtns_print_version
    
    subroutine shtns_print_cfg(shtns) bind(C, name='shtns_print_cfg')
      import
      type(C_PTR) :: shtns
    end subroutine shtns_print_cfg
    
    type(C_PTR) function shtns_init(flags,lmax,mmax,mres,nlat,nphi) bind(C, name='shtns_init')
      import
      integer(C_INT), value :: flags
      integer(C_INT), value :: lmax
      integer(C_INT), value :: mmax
      integer(C_INT), value :: mres
      integer(C_INT), value :: nlat
      integer(C_INT), value :: nphi
    end function shtns_init
    
    type(C_PTR) function shtns_create(lmax,mmax,mres,norm) bind(C, name='shtns_create')
      import
      integer(C_INT), value :: lmax
      integer(C_INT), value :: mmax
      integer(C_INT), value :: mres
      integer(C_INT), value :: norm
    end function shtns_create
    
    subroutine shtns_set_grid(shtns,flags,eps,nlat,nphi) bind(C, name='shtns_set_grid')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: flags
      real(C_DOUBLE), value :: eps
      integer(C_INT), value :: nlat
      integer(C_INT), value :: nphi
    end subroutine shtns_set_grid

    subroutine shtns_set_grid_auto(shtns,flags,eps,nl_order,nlat,nphi) bind(C, name='shtns_set_grid_auto')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: flags
      real(C_DOUBLE), value :: eps
      integer(C_INT), value :: nl_order
      integer(C_INT), value :: nlat
      integer(C_INT), value :: nphi
    end subroutine shtns_set_grid_auto
    
    type(C_PTR) function shtns_create_with_grid(shtns,mmax,nofft) bind(C, name='shtns_create_with_grid')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: mmax
      integer(C_INT), value :: nofft
    end function shtns_create_with_grid
    
    integer(C_INT) function shtns_use_threads(num_threads) bind(C, name='shtns_use_threads')
      import
      integer(C_INT), value :: num_threads
    end function shtns_use_threads
    
    subroutine shtns_use_gpu(device_id) bind(C, name='shtns_use_gpu')
      import
      integer(C_INT), value :: device_id
    end subroutine shtns_use_gpu
    
    subroutine shtns_reset() bind(C, name='shtns_reset')
      import
    end subroutine shtns_reset
    
    subroutine shtns_destroy(shtns) bind(C, name='shtns_destroy')
      import
      type(C_PTR), value :: shtns
    end subroutine shtns_destroy
    
    subroutine shtns_unset_grid(shtns) bind(C, name='shtns_unset_grid')
      import
      type(C_PTR), value :: shtns
    end subroutine shtns_unset_grid
    
    subroutine shtns_robert_form(shtns,robert) bind(C, name='shtns_robert_form')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: robert
    end subroutine shtns_robert_form
    
    type(C_PTR) function shtns_malloc(bytes) bind(C, name='shtns_malloc')
      import
      integer(C_SIZE_T), value :: bytes
    end function shtns_malloc
    
    subroutine shtns_free(p) bind(C, name='shtns_free')
      import
      type(C_PTR), value :: p
    end subroutine shtns_free
    
    real(C_DOUBLE) function sh00_1(shtns) bind(C, name='sh00_1')
      import
      type(C_PTR), value :: shtns
    end function sh00_1
    
    real(C_DOUBLE) function sh10_ct(shtns) bind(C, name='sh10_ct')
      import
      type(C_PTR), value :: shtns
    end function sh10_ct
    
    real(C_DOUBLE) function sh11_st(shtns) bind(C, name='sh11_st')
      import
      type(C_PTR), value :: shtns
    end function sh11_st
    
    real(C_DOUBLE) function shlm_e1(shtns,l,m) bind(C, name='shlm_e1')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: l
      integer(C_INT), value :: m
    end function shlm_e1
    
    subroutine shtns_gauss_wts(shtns,wts) bind(C, name='shtns_gauss_wts')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), intent(inout) :: wts(*)
    end subroutine
    
    subroutine SH_Yrotate(shtns,Qlm,alpha,Rlm) bind(C, name='SH_Yrotate')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      real(C_DOUBLE), value :: alpha
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Rlm(*)
    end subroutine SH_Yrotate
    
    subroutine SH_Yrotate90(shtns,Qlm,Rlm) bind(C, name='SH_Yrotate90')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Rlm(*)
    end subroutine SH_Yrotate90
    
    subroutine SH_Xrotate90(shtns,Qlm,Rlm) bind(C, name='SH_Xrotate90')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Rlm(*)
    end subroutine SH_Xrotate90
    
    subroutine mul_ct_matrix(shtns,mx) bind(C, name='mul_ct_matrix')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), dimension(*) :: mx
    end subroutine mul_ct_matrix
    
    subroutine st_dt_matrix(shtns,mx) bind(C, name='st_dt_matrix')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), dimension(*) :: mx
    end subroutine st_dt_matrix
    
    subroutine SH_mul_mx(shtns,mx,Qlm,Rlm) bind(C, name='SH_mul_mx')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), dimension(*) :: mx
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Rlm(*)
    end subroutine SH_mul_mx

    subroutine spat_to_SH(shtns, Vr, Qlm) bind(C, name='spat_to_SH')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), intent(inout) :: Vr(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Qlm(*)
    end subroutine spat_to_SH
    
    subroutine SH_to_spat(shtns,Qlm,Vr) bind(C, name='SH_to_spat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(in) :: Qlm(*)
      real(C_DOUBLE), intent(out) :: Vr(*)
    end subroutine SH_to_spat
    
    subroutine SH_to_spat_cplx(shtns,alm,z) bind(C, name='SH_to_spat_cplx')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: alm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: z(*)
    end subroutine SH_to_spat_cplx
    
    subroutine spat_cplx_to_SH(shtns,z,alm) bind(C, name='spat_cplx_to_SH')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: z(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: alm(*)
    end subroutine spat_cplx_to_SH
    
    subroutine spat_to_SHsphtor(shtns,Vt,Vp,Slm,Tlm) bind(C, name='spat_to_SHsphtor')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), intent(inout) :: Vt(*)
      real(C_DOUBLE), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tlm(*)
    end subroutine spat_to_SHsphtor
    
    subroutine SHsphtor_to_spat(shtns,Slm,Tlm,Vt,Vp) bind(C, name='SHsphtor_to_spat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(in) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(in) :: Tlm(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
    end subroutine SHsphtor_to_spat

    subroutine SHsph_to_spat(shtns,Slm,Vt,Vp) bind(C, name='SHsph_to_spat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(in) :: Slm(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
    end subroutine SHsph_to_spat
    
    subroutine SHtor_to_spat(shtns,Tlm,Vt,Vp) bind(C, name='SHtor_to_spat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(in) :: Tlm(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
    end subroutine SHtor_to_spat

    subroutine SH_to_grad_spat(shtns,Qlm,Gt,Gp) bind(C, name='SHsph_to_spat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(in) :: Qlm(*)
      real(C_DOUBLE), intent(out) :: Gt(*)
      real(C_DOUBLE), intent(out) :: Gp(*)
    end subroutine SH_to_grad_spat
    
    subroutine spat_cplx_to_SHsphtor(shtns,Vt,Vp,Slm,Tlm) bind(C, name='spat_cplx_to_SHsphtor')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tlm(*)
    end subroutine spat_cplx_to_SHsphtor
    
    subroutine SHsphtor_to_spat_cplx(shtns,Slm,Tlm,Vt,Vp) bind(C, name='SHsphtor_to_spat_cplx')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vp(*)
    end subroutine SHsphtor_to_spat_cplx
    
    subroutine SHqst_to_spat(shtns,Qlm,Slm,Tlm,Vr,Vt,Vp) bind(C, name='SHqst_to_spat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      real(C_DOUBLE), intent(out) :: Vr(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
    end subroutine SHqst_to_spat

    subroutine spat_to_SHqst(shtns,Vr,Vt,Vp,Qlm,Slm,Tlm) bind(C, name='spat_to_SHqst')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), intent(inout) :: Vr(*)
      real(C_DOUBLE), intent(inout) :: Vt(*)
      real(C_DOUBLE), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tlm(*)
    end subroutine spat_to_SHqst
    
    subroutine spat_cplx_to_SHqst(shtns,Vr,Vt,Vp,Qlm,Slm,Tlm) bind(C, name='spat_cplx_to_SHqst')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vr(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tlm(*)
    end subroutine spat_cplx_to_SHqst
    
    subroutine SHqst_to_spat_cplx(shtns,Qlm,Slm,Tlm,Vr,Vt,Vp) bind(C, name='SHqst_to_spat_cplx')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vr(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vp(*)
    end subroutine SHqst_to_spat_cplx
    
    subroutine spat_to_SH_l(shtns,Vr,Qlm,ltr) bind(C, name='spat_to_SH_l')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), intent(inout) :: Vr(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Qlm(*)
      integer(C_INT), value :: ltr
    end subroutine spat_to_SH_l
    
    subroutine SH_to_spat_l(shtns,Qlm,Vr,ltr) bind(C, name='SH_to_spat_l')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      real(C_DOUBLE), intent(out) :: Vr(*)
      integer(C_INT), value :: ltr
    end subroutine SH_to_spat_l
    
    subroutine SHsphtor_to_spat_l(shtns,Slm,Tlm,Vt,Vp,ltr) bind(C, name='SHsphtor_to_spat_l')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHsphtor_to_spat_l
    
    subroutine SHsph_to_spat_l(shtns,Slm,Vt,Vp,ltr) bind(C, name='SHsph_to_spat_l')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHsph_to_spat_l
    
    subroutine SHtor_to_spat_l(shtns,Tlm,Vt,Vp,ltr) bind(C, name='SHtor_to_spat_l')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHtor_to_spat_l
    
    subroutine spat_to_SHsphtor_l(shtns,Vt,Vp,Slm,Tlm,ltr) bind(C, name='spat_to_SHsphtor_l')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), intent(inout) :: Vt(*)
      real(C_DOUBLE), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tlm(*)
      integer(C_INT), value :: ltr
    end subroutine spat_to_SHsphtor_l
    
    subroutine spat_to_SHqst_l(shtns,Vr,Vt,Vp,Qlm,Slm,Tlm,ltr) bind(C, name='spat_to_SHqst_l')
      import
      type(C_PTR), value :: shtns
      real(C_DOUBLE), intent(inout) :: Vr(*)
      real(C_DOUBLE), intent(inout) :: Vt(*)
      real(C_DOUBLE), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tlm(*)
      integer(C_INT), value :: ltr
    end subroutine spat_to_SHqst_l
    
    subroutine SHqst_to_spat_l(shtns,Qlm,Slm,Tlm,Vr,Vt,Vp,ltr) bind(C, name='SHqst_to_spat_l')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      real(C_DOUBLE), intent(out) :: Vr(*)
      real(C_DOUBLE), intent(out) :: Vt(*)
      real(C_DOUBLE), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHqst_to_spat_l
    
    subroutine spat_to_SH_ml(shtns,im,Vr,Ql,ltr) bind(C, name='spat_to_SH_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vr(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Ql(*)
      integer(C_INT), value :: ltr
    end subroutine spat_to_SH_ml
    
    subroutine SH_to_spat_ml(shtns,im,Ql,Vr,ltr) bind(C, name='SH_to_spat_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Ql(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vr(*)
      integer(C_INT), value :: ltr
    end subroutine SH_to_spat_ml
    
    subroutine spat_to_SHsphtor_ml(shtns,im,Vt,Vp,Sl,Tl,ltr) bind(C, name='spat_to_SHsphtor_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Sl(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tl(*)
      integer(C_INT), value :: ltr
    end subroutine spat_to_SHsphtor_ml
    
    subroutine SHsphtor_to_spat_ml(shtns,im,Sl,Tl,Vt,Vp,ltr) bind(C, name='SHsphtor_to_spat_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Sl(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tl(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHsphtor_to_spat_ml
    
    subroutine SHsph_to_spat_ml(shtns,im,Sl,Vt,Vp,ltr) bind(C, name='SHsph_to_spat_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Sl(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHsph_to_spat_ml
    
    subroutine SHtor_to_spat_ml(shtns,im,Tl,Vt,Vp,ltr) bind(C, name='SHtor_to_spat_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tl(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHtor_to_spat_ml
    
    subroutine spat_to_SHqst_ml(shtns,im,Vr,Vt,Vp,Ql,Sl,Tl,ltr) bind(C, name='spat_to_SHqst_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vr(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Vp(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Ql(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Sl(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Tl(*)
      integer(C_INT), value :: ltr
    end subroutine spat_to_SHqst_ml
    
    subroutine SHqst_to_spat_ml(shtns,im,Ql,Sl,Tl,Vr,Vt,Vp,ltr) bind(C, name='SHqst_to_spat_ml')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), value :: im
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Ql(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Sl(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tl(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vr(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vt(*)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vp(*)
      integer(C_INT), value :: ltr
    end subroutine SHqst_to_spat_ml
    
    subroutine SH_to_grad_point(shtns,DrSlm,Slm,cost,phi,vr,vt,vp) bind(C, name='SH_to_grad_point')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: DrSlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      real(C_DOUBLE), value :: cost
      real(C_DOUBLE), value :: phi
      real(C_DOUBLE), intent(out) :: vr(*)
      real(C_DOUBLE), intent(out) :: vt(*)
      real(C_DOUBLE), intent(out) :: vp(*)
    end subroutine SH_to_grad_point
    
    subroutine SHqst_to_point(shtns,Qlm,Slm,Tlm,cost,phi,vr,vt,vp) bind(C, name='SHqst_to_point')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      real(C_DOUBLE), value :: cost
      real(C_DOUBLE), value :: phi
      real(C_DOUBLE), intent(out) :: vr(*)
      real(C_DOUBLE), intent(out) :: vt(*)
      real(C_DOUBLE), intent(out) :: vp(*)
    end subroutine SHqst_to_point
    
    subroutine SH_to_lat(shtns,Qlm,cost,vr,nphi,ltr,mtr) bind(C, name='SH_to_lat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      real(C_DOUBLE), value :: cost
      real(C_DOUBLE), intent(out) :: vr(*)
      integer(C_INT), value :: nphi
      integer(C_INT), value :: ltr
      integer(C_INT), value :: mtr
    end subroutine SH_to_lat
    
    subroutine SHqst_to_lat(shtns,Qlm,Slm,Tlm,cost,vr,vt,vp,nphi,ltr,mtr) bind(C, name='SHqst_to_lat')
      import
      type(C_PTR), value :: shtns
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Qlm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Slm(*)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: Tlm(*)
      real(C_DOUBLE), value :: cost
      real(C_DOUBLE), intent(out) :: vr(*)
      real(C_DOUBLE), intent(out) :: vt(*)
      real(C_DOUBLE), intent(out) :: vp(*)
      integer(C_INT), value :: nphi
      integer(C_INT), value :: ltr
      integer(C_INT), value :: mtr
    end subroutine SHqst_to_lat


    integer(C_INT) function shtns_lmidx(shtns,l,m) bind(C, name='shtns_lmidx_fortran')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), intent(in) :: l
      integer(C_INT), intent(in) :: m
    end function shtns_lmidx

    integer(C_INT) function shtns_lm2l(shtns,lm) bind(C, name='shtns_lm2l_fortran')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), intent(in) :: lm
    end function shtns_lm2l

    integer(C_INT) function shtns_lm2m(shtns,lm) bind(C, name='shtns_lm2m_fortran')
      import
      type(C_PTR), value :: shtns
      integer(C_INT), intent(in) :: lm
    end function shtns_lm2m

  end interface
