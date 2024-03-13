program SHT_example
   !
   ! This is a modern Fortran example that shows how to set-up the spherical
   ! harmonic transforms using SHTns.
   !

   use iso_c_binding
   use iso_fortran_env, only: real64

   implicit none

   include 'shtns.f03'

   integer, parameter :: dp=real64
   real(dp), parameter :: pi=acos(-1.0_dp)

   type(shtns_info), pointer :: shtns
   real(dp), pointer :: cosTheta(:), sinTheta(:)
   type(c_ptr) :: shtns_c

   integer :: lmax, mmax, mres, nthreads
   integer :: nlat, nphi, n_p, lm, l, m, lmStart, lmStop
   integer :: nlm, norm, nout, layout
   real(dp) :: eps_polar
   real(dp), allocatable :: Sh(:,:)
   complex(dp), allocatable :: Slm(:), ShT(:)


   !--  Set the size of the transforms
   lmax = 5
   mmax = 3
   mres = 1
   nphi = 10
   nlat = 64

   !-- Polar optimisation threshold
   eps_polar = 1.0e-10_dp

   !-- Norm and layout for SHTs
   norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
   layout = SHT_GAUSS + SHT_PHI_CONTIGUOUS

   call shtns_verbose(2)
   nthreads = shtns_use_threads(0)
   print*, 'nthreads=', nthreads

   shtns_c = shtns_create(lmax, mmax, mres, norm)
   call shtns_set_grid(shtns_c, layout, eps_polar, nlat, nphi)

   !-- C/Fortran pointer mapping
   call c_f_pointer(cptr=shtns_c, fptr=shtns)
   call c_f_pointer(cptr=shtns%ct, fptr=cosTheta, shape=[shtns%nlat])
   call c_f_pointer(cptr=shtns%st, fptr=sinTheta, shape=[shtns%nlat])

   print*, 'cosTheta', cosTheta
   print*, 'sinTheta', sinTheta

   print*, 'NLM=', shtns%nlm

   allocate( Sh(1:shtns%nphi,1:shtns%nlat), Slm(1:shtns%nlm), ShT(1:shtns%nlat) )
   Slm(:)=0.0_dp

   !-- Get lm index from the (l,m) pair
   l = 1 
   m = 0
   lm = shtns_lmidx(shtns_c,l,m)
   !-- Set the l=1, m=0 mode to 1.0
   Slm(lm)=(1.0_dp,0.0_dp)

   !-- Get (l,m) from lm
   lm = 17
   l = shtns_lm2l(shtns_c,lm)
   m = shtns_lm2m(shtns_c,lm)
   print*, '(l,m)=', l,m

   !-- Spec -> Spat
   call SH_to_spat(shtns_c, Slm, Sh)

   print*, 'Y(1,0)', Sh(1,:) ! It should be sqrt(3/4/pi) * cos(theta)

   !-- Spat -> Spec,  DOES NOT PRESERVE THE CONTENT OF INPUT 'Sh'
   call spat_to_SH(shtns_c, Sh, Slm)

   !-- print S(1,0) to check it 1.0 is recovered.
   print*, 'S(1,0)', Slm(lm)

   !-- Legendre only for m=0
   lmStart = shtns_lmidx(shtns_c,0,0)
   lmStop = shtns_lmidx(shtns_c,lmax,0)
   call SH_to_spat_ml(shtns_c, 0, Slm(lmStart:lmStop), ShT, lmax)
   print*, 'Legendre only', ShT

   deallocate(Sh, Slm, ShT)
   call shtns_unset_grid(shtns_c)
   call shtns_destroy(shtns_c)

end program
