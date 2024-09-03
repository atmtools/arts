!
!  Copyright (c) 2010-2018 Centre National de la Recherche Scientifique.
!  written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
!  
!  nathanael.schaeffer@univ-grenoble-alpes.fr
!  
!  This software is governed by the CeCILL license under French law and
!  abiding by the rules of distribution of free software. You can use,
!  modify and/or redistribute the software under the terms of the CeCILL
!  license as circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info".
!  
!  The fact that you are presently reading this means that you have had
!  knowledge of the CeCILL license and that you accept its terms.
!  


!! A Fortran example program that performs backward and forward Spherical Harmonic Transforms using SHTns
      PROGRAM SHT_example

      IMPLICIT NONE
! import useful parameters for shtns initialization
      include 'shtns.f'

      integer*4 lmax, mmax, mres
      integer*4 nlat, nphi
      integer*4 nlm
      integer*4 layout
      integer*4 norm
      real*8 eps_polar
      complex*16, allocatable :: Slm(:), Tlm(:)
      real*8, allocatable :: Sh(:,:), Th(:,:)

      integer i,lm, m

! set size of transform
      lmax = 5
      mmax = 2
      mres = 2
      nphi = 6
      nlat = 8

! compute sizes required for arrays.
      call shtns_calc_nlm(nlm, lmax, mmax, mres)
      print*,'NLM=',nlm

! displays information during initialization
      call shtns_verbose(1)

! enable multi-threaded transform (OpenMP) if supported.
      call shtns_use_threads(0)

! init SHT. SHT_PHI_CONTIGUOUS is defined in 'shtns.f'
      layout = SHT_PHI_CONTIGUOUS
      call shtns_init_sh_gauss(layout, lmax, mmax, mres, nlat, nphi)

! alternatively, you can use the two following calls, giving more control on the initialization process,
! namely you can choose a normalization with 'norm' and control the polar optimization
! with 'eps_polar' : from 0 (no optimization) to 1e-6 (agressive optimization).
!      norm = SHT_ORTHONORMAL + SHT_REAL_NORM
!      call shtns_set_size(lmax, mmax, mres, norm)
!      eps_polar = 1.e-10
!      call shtns_precompute(SHT_GAUSS, layout, eps_polar, nlat, nphi)

! display information:
       call shtns_print_cfg()

! allocate memory for spectral and spatial representation
      allocate ( Slm(1:nlm), Tlm(1:nlm) )
      allocate ( Sh(1:nphi,1:nlat), Th(1:nphi,1:nlat) )
      Slm(:) = 0.0

! get index for given l and m
      call shtns_lmidx(lm, 1, 0)
      Slm(lm) = 1.0
      print*
      print*,Slm(:)

      call shtns_sh_to_spat(Slm, Sh)

! print all spatial data
      print*
      do i=1,nphi
        print*,'phi index=',i
        print*,Sh(i,:)
      enddo

! spatial to spectral transform. DOES NOT PRESERVE THE CONTENT OF INPUT 'Sh'
      call shtns_spat_to_sh(Sh, Slm)

      print*
      print*,Slm(:)

! print SH data, grouping by m  (shows how to access SH coefficients)
      print*
      do m=0, mres*mmax, mres
        print*,'m=', m
        call shtns_lmidx(lm,m,m)
        print*,Slm(lm : lm+(lmax-m))
      enddo

! Legendre transform (fixed m, no fft):
! we need a larger spatial grid for these functions to work (nlat >= 16)
      nlat = 32 
      call shtns_init_sh_gauss(layout, lmax, mmax, mres, nlat, nphi)
      m = 1*mres
      call shtns_lmidx(lm, m, m)
      Slm(lm+1) = 1.0
      call shtns_sh_to_spat_ml(m/mres, Slm(lm), Tlm, lmax)
      print*
      print*,'spectral m=', m
      print*,Slm(lm : lm+lmax-m+1)
      print*,'spatial m=', m
      print*,Tlm(1 : nlat)

      stop
      END
