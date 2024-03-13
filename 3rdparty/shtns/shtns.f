!
!  Copyright (c) 2010-2015 Centre National de la Recherche Scientifique.
!  written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
!  
!  nathanael.schaeffer@ujf-grenoble.fr
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




      INTEGER SHT_NATIVE_LAYOUT
      PARAMETER (SHT_NATIVE_LAYOUT=0)
      INTEGER SHT_THETA_CONTIGUOUS
      PARAMETER (SHT_THETA_CONTIGUOUS=256)
      INTEGER SHT_PHI_CONTIGUOUS
      PARAMETER (SHT_PHI_CONTIGUOUS=512)

      INTEGER SHT_NO_CS_PHASE
      PARAMETER (SHT_NO_CS_PHASE=1024)
      INTEGER SHT_REAL_NORM
      PARAMETER (SHT_REAL_NORM=2048)

      INTEGER SHT_ORTHONORMAL
      PARAMETER (SHT_ORTHONORMAL=0)
      INTEGER SHT_FOURPI
      PARAMETER (SHT_FOURPI=1)
      INTEGER SHT_SCHMIDT
      PARAMETER (SHT_SCHMIDT=2)

      INTEGER SHT_GAUSS
      PARAMETER (SHT_GAUSS=0)
      INTEGER SHT_AUTO
      PARAMETER (SHT_AUTO=1)
      INTEGER SHT_REG_FAST
      PARAMETER (SHT_REG_FAST=2)
      INTEGER SHT_REG_DCT
      PARAMETER (SHT_REG_DCT=3)
      INTEGER SHT_QUICK_INIT
      PARAMETER (SHT_QUICK_INIT=4)
      INTEGER SHT_REG_POLES
      PARAMETER (SHT_REG_POLES=5)
      INTEGER SHT_GAUSS_FLY
      PARAMETER (SHT_GAUSS_FLY=6)

      INTEGER SHT_SOUTH_POLE_FIRST
      PARAMETER (SHT_SOUTH_POLE_FIRST=8192)
      INTEGER SHT_SCALAR_ONLY
      PARAMETER (SHT_SCALAR_ONLY=4096)
      INTEGER SHT_LOAD_SAVE_CFG
      PARAMETER (SHT_LOAD_SAVE_CFG=16384)
