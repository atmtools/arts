!--------------------------------------------------------------------------------------------------------------------
module module_phsub
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to 
!
    interface

        subroutine PopuCAL(dta1, nLines, molP,econ)
            use module_common_var
            implicit none
            integer*8              , intent(in   ) :: nLines
            type (dta_MOL)         , intent(in   ) :: molP
            type (dta_SDF), pointer, intent(inout) :: dta1
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine PopuCAL

        subroutine DipCAL(dta1,nLines,molP,econ)
            use module_common_var
            use module_molecSp
            use module_maths
            use module_error
            implicit none
            integer*8               , intent(in)    :: nLines
            type  (dta_SDF), pointer, intent(inout) :: dta1
            type  (dta_ERR)         , intent(inout) :: econ
            type  (dta_MOL)         , intent(in)    :: molP
        end subroutine DipCAL
        
        double precision function pureHund(caseHund, J_i, N_i, J_f, N_f)
            implicit none
            integer (kind=8), intent(in) :: N_i, N_f
            real*8          , intent(in) :: J_i, J_f
            character       , intent(in) :: caseHund
        end function pureHund

        double precision function Kpart1( Ji_p, Jf, Jf_p, li, lf, AF1)
            use module_common_var
            implicit none
            integer*8, intent(in)        :: li, lf
            double precision, intent(in) :: Ji_p, Jf, Jf_p, AF1
        end function Kpart1

        double precision function Kpart1_O2(Ji, Ji_p, Jf, Jf_p, Ni, Ni_p, Nf, Nf_p, &
                                            AF1, econ)
            use module_common_var
            implicit none
            integer*8, intent(in)        :: Ni, Ni_p, Nf, Nf_p
            double precision, intent(in) :: Ji, Ji_p, Jf, Jf_p, AF1 
            type (dta_ERR), intent(inout):: econ
        end function Kpart1_O2

        double precision function Kpart2(L, Ji, Ji_p, Jf, Jf_p, li, lf, AF2, econ)
            use module_common_var
            use module_error
            use module_maths
            implicit none
            integer*8, intent(in)        :: L, li, lf
            double precision, intent(in) :: Ji, Ji_p, Jf, Jf_p, AF2
            type (dta_ERR), intent(inout):: econ
        end function Kpart2

        double precision function Kpart2_O2(L, Ji, Ji_p, Jf, Jf_p, &
                                            Ni, Ni_p, Nf, Nf_p, Si, Sf, &
                                            AF2,econ)
            use module_common_var
            use module_error
            use module_maths
            implicit none
            integer*8, intent(in)        :: L, Ni, Ni_p, Nf, Nf_p
            double precision, intent(in) :: Ji, Ji_p, Jf, Jf_p, Si, Sf, AF2
            type (dta_ERR), intent(inout):: econ
        end function Kpart2_O2
        
        double precision function K_jkCalc(nt,j,k,h,nLines,molP,PerM,econ) 
            use module_common_var
            use module_error
            use module_molecSp
            use module_maths
            implicit none
            integer*8               , intent(in   ) :: nt,j,k,nLines
            type (dta_SDF), pointer , intent(in   ) :: h
            type (dta_MOL)          , intent(in   ) :: molP, PerM
            type (dta_ERR)          , intent(inout) :: econ
        end function K_jkCalc

        double precision function K_jkO2(j,k,h,nLines,molP,PerM,econ)
            use module_common_var
            use module_error
            use module_molecSp
            use module_maths
            implicit none
            integer*8              , intent(in   ) :: j,k,nLines
            type (dta_SDF), pointer, intent(in   ) :: h
            type (dta_MOL)         , intent(in   ) :: molP, PerM
            type (dta_ERR)         , intent(inout) :: econ

        end function K_jkO2
        
        double precision function Ql_mol_X(molP,L)
            use module_common_var
            use module_molecSp
            use module_maths
            implicit none
            integer*8    , intent(in) :: L
            type(dta_MOL), intent(in) :: molP
        end function Ql_mol_X

        double precision function Ql_mol_LLS(J1,J2,v1,v2,molP,L,econ)
            use module_common_var
            use module_molecSp
            use module_maths
            implicit none
            integer*8    , intent(in   ) :: L
            real*8       , intent(in   ) :: J1,J2,v1,v2
            type(dta_MOL), intent(in   ) :: molP
            type(dta_ERR), intent(inout) :: econ
        end function Ql_mol_LLS

        double precision function Ql_mol_AF_LLS(J1,J2,v1,v2,molP,PerM,L)
            use module_common_var
            use module_molecSp
            use module_maths
            implicit none
            integer*8    , intent(in) :: L
            type(dta_MOL), intent(in) :: molP,PerM
            real*8       , intent(in) :: J1,J2,v1,v2
        end function Ql_mol_AF_LLS
        
        double precision function AFmol_X(molP, PerM, L, step)
            use module_common_var
            implicit none
            integer*8     , intent(in) :: step
            real*8        , intent(in) :: L
            type (dta_MOL), intent(in) :: molP, PerM
        end function AFmol_X

        double precision function c_AF(molP, PerM, J, L, step)
            use module_common_var
            implicit none
            integer*8     , intent(in) :: step
            real*8        , intent(in) :: J,L
            type (dta_MOL), intent(in) :: molP, PerM
        end function c_AF

        subroutine sumRule(nLines,indexS,dipole,Wmat,dfact,econ)
            use module_common_var
            use module_error
            implicit none
            integer (kind=8), intent(in)    :: nLines
            integer (kind=8), intent(in)    :: indexS(nLines)
            real*8          , intent(in)    :: dfact
            Double Precision, intent(in)    :: dipole(nLines)
            Double Precision, intent(in)    :: Wmat(nLines,nLines)
            type (dta_ERR)  , intent(inout) :: econ
        end subroutine sumRule

    end interface

END module module_phsub
!--------------------------------------------------------------------------------------------------------------------
! RELAXATION MATRIX SUBROUTINES -------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE PopuCAL(dta1, nLines, molP,econ)
!--------------------------------------------------------------------------------------------------------------------
!"PopuCAL": POPUlation CALculation
!
! Detailed Description:
! ---------------------
! It calcules the Populations of the Lower state of the 
! Line/transition k at Temperature T0. Elements comes 
! from Subroutine "Readline".          
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
!   k       : marks the transition in study (Input).
!   E       : lower state energy (cm-1).
!   SWeig00 : State Statistical Weight (also called g00 in HITRAN). UNITLESS.
! molN    : = Na·Nmcon; molecular number density (molecules·cm-3)
! PFmol   : internal partition function. UNITLESS.
! ptype   : this parameter marks the result Population calculation 
!           in other words:
!           ('hit') Acording to HITRAN96 
!           ('tra') Acording to Hartmann et. al 2006 
!
! Other important Input Quantities (through "module_common_var")
! --------------------------------
!   c2      : second radiation constant = h·c/kb = 1.4388 cm·K
!   T0      : Reference Temperature in Kelvin, 296 (K, Input).
!
! Output Quantities (through Common Statements)
! -----------------
!   PopuL0  : Populations of the Lower state of the Lines
!             at Temperature T0                 (Output)
!           NOTE: 
!             1) with Ha Tran Code, the Population is UNITLESS.
!             2) HITRAN definition:
!                PopuT0(k) = molN*g00(k)*dexp(-c2*E(k)/T0)/(PFmol(iso,T0))
!                give the units of the molecular number density (molN)
!
! Accessed Files:  None
! --------------
!
! Called Routines: 'none'
! ---------------  
!
! Called By: 'DipCAL' (Dipole of the lower state CALculation)
! ---------
!
!   Double Precision Version
!
!   T.Mendaza, last change 25 Jan 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_molecSp
    Implicit None
    integer*8              , intent(in)    :: nLines
    type (dta_MOL)         , intent(in)    :: molP
    type (dta_SDF), pointer, intent(inout) :: dta1
    type (dta_ERR)         , intent(inout) :: econ
    integer*8                     :: iso
    integer*8                     :: i, j, k
    Double Precision              :: cte1,cte2, Erot_red, T
    Double Precision              :: molN, PFmol, PFr
!-----------------------------------------
    T = molP % Temp
!-----------------------------------------
!
! Partition function ratio:
!
    !old 
    !PFr=PFmol(molP,T0)/PFmol(molP,T)
    ! QT0 = PFmol(molP,T0)
    ! QT  = PFmol(molP,T )
    ! Now the partition function came as an input from "arts_interface"
    PFr=molP%QT0 / molP%QT
!
    DO k=1,nLines
        !Erot_red = dta1%J(k,1)*( dta1%J(k,1) + 1. )
        !cte1 = -c2*molP%B0*Erot_red/T0
        cte1 = -c2*dta1%E(k)/T0
        cte2 = -c2*dta1%E(k)*(1.d0/T-1.d0/T0)
        !
        dta1%PopuT0(k) = dta1%swei00(k)*dexp(cte1)/molP%QT0
        dta1%PopuT(k)  = dta1%PopuT0(k)*PFr*dexp(cte2)
    ENDDO
!
    IF (ptype .eq. 'tra') then
        ! the Population remains as it has been calculated.
    ELSEIF (ptype .eq. 'hit') then
        ! Here we follow the HITRAN 1996 definition of 
        ! Population of the lower level using the molecular number density:
        ! molN = Na·Nmcon
        ! CH4 example (Nmcon = Nch4):
        ! Na   = 6.022·10^23 molecules·mol-1; Avogadro constant
        ! Nch4 = 41.245 mol/m3 = 41.245 mol/(1E006 *cm3) 
        !      = 41.245*1E-06  mol/cm3 (value used in common_var_module)
        ! NIST website:
        ! http://webbook.nist.gov/
        !
        molN= Na*molP%Nmcon ! molecule/cm3
        !molN= Nmcon    ! mol/cm3
        !
        DO k=1,nLines
            dta1%PopuT0(k) = molN*dta1%PopuT0(k)
            dta1%PopuT(k)  = molN*dta1%PopuT(k)
        ENDDO
    ELSE
        call errorPType(ptype, econ)
    ENDIF
!
    Return
  END SUBROUTINE PopuCAL
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE DipCAL(dta1,nLines,molP,econ)
!--------------------------------------------------------------------------------------------------------------------
! "DipCAL": DIPole CALculation
! 
! Detailed Description:
! ---------------------
! This subroutine computes the dipole element d_k for a certain level 
! obtaining parameters from Subroutine "Readline". 
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
!	E         : lower state energy (cm-1)
!	Sig       : WaveNumber of the Computation (cm-1, Input)
!	PopuT0    : Populations of the Lower state of the Lines
!	            at Temperature T0 (input)
!
! Other important Input Constant (through "module_common_var")
! ------------------------------
!	T0	      : Reference Temperature in Kelvin, 296K (Input).
!	c2        : second radiation constant = hc/k = 1.4388 cm·K
!pure_b_case: --
!
! Output Quantities (through DTA structure)
! -----------------
!	DipoT 	  : Dipole transition Moments of the (cm/Molecule**0.5).
! DipoT0    : Dipole transition Moments of the Line at T0 [(J·cm3)**0.5].
! D0        : Reduce Matrix Element.
! Drigrotor : Rigid Rotor Dipole.
!
!
! Accessed Files:  None
! --------------
!
! Called Routines: None
! ---------------  
!
! Called By: Main interface
! ---------
!
!	Double Precision Version
!
!	T.Mendaza, last change 17 February 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_molecSp
    use module_maths
    use module_error
    Implicit None
    integer*8              , intent(in   ) :: nLines
    type (dta_SDF), pointer, intent(inout) :: dta1
    type (dta_ERR)         , intent(inout) :: econ
    type (dta_MOL)         , intent(in   ) :: molP
    integer*8                   :: N_i, N_f
    real*8                      :: J_i, J_f
    integer*8                   :: i, j, k
    real*8                      :: s_i, s_f
    double precision            :: cte,fA,fB, w3j, w6j
    double precision            :: purehund, d0
    double precision            :: signDIP, Ia_hit
!
!----------
! Compute the Dipole moment at every line 
! (Diatomic approximation... 
!
! INIT. VAR:
!-----------
    DO k = 1, nLines
      !
      ! DIPOLE MOMENT 
      !
       ! Hartmann's style (Using LINE-INTENSITY):
        fA  = -c2*dta1%Sig(k)/T0
        fB  = 1.D0 - dexp(fA)
        !dta1%DipoT(k)= dsqrt( dta1%Str(k)/(dta1%Sig(k)*dta1%PopuT0(k)*fB) )
        dta1%DipoT(k)= dsqrt( dta1%Str(k)/(dta1%Sig(k)*dta1%PopuT(k)*fB) )
            !cm/molecules**0.5
    !
    ! REDUCE DIPOLE MOMENT, D0
    !
        if ( molP%M .eq. 7 ) then
        ! O2 possible Dipole'scalculation:
        !
        ! a) Makarov et al. (2013):
        ! DATA:
            J_i  = dta1%J(k,1) !J_i  = J LOWER L.
            J_f  = dta1%J(k,2) !J_f  = J UPPER L.
            N_i  = dta1%N(k,1) !LOWER L.
            N_f  = dta1%N(k,2) !UPPER L.
            !
            !
            !if ( mode .eq. 'mak1') then
            ! a) [Makarov et al. 2013]
            !    signDIP = (-1)**(int(J_f)+N_i)
            !    w6j = wigner6j( 1.d0, 1.d0, 1.d0, J_i, J_f, real(N_i,dp))
            !    d0 = dsqrt(6.d0*(2*J_f+1.d0)*(2*J_i+1.d0))
            !    dta1%D0(k) = signDIP*w6j*d0
            !    !
            !elseif (mode .eq. 'mak2') then
            !
            !    signDIP = (-1)**(int(J_f)+N_f)
            !    w6j = wigner6j( 1.d0, 1.d0, 1.d0, J_i, J_f, real(N_f,dp))
            !    d0 = dsqrt(6.d0*(2*J_f+1.d0)*(2*J_i+1.d0))
            !    dta1%D0(k) = signDIP*w6j*d0 
            !    !
            !elseif (mode .eq. 'tran') then
            ! b) [Tran et al. 2006]    
                d0 = pureHund(caseHund, J_i, N_i, J_f, N_f)  
                signDIP = (d0/abs(d0))
                dta1%D0(k) = signDIP*dta1%DipoT(k)
            !endif
        else
            !Other molecules:
            if ( mod((int(J_f)+dta1%lv2(2)),2) .eq. 0) then
            ! j = (-1)**(J_f+lf)
                j = 1
            else
                j = -1
            endif
        ! K_t = marks the dimension of the tensor.
        ! isotropic Raman = 0 
        !   IR ----> K_t = 1
        !   Raman -> K_t = 2
        ! wigner simbols are implemented for -> IR case.
            w3j = wigner3j( dta1%J(k,1), real(K_t,dp), dta1%J(k,2), &
                        real(dta1%lv2(1), dp), real( dta1%lv2(2)-dta1%lv2(1),dp), &
                        real(-dta1%lv2(2), dp))
            d0 = j*dsqrt(2.d0*dta1%J(k,2) + 1.d0)*w3j
            dta1%D0(k) = abs(d0)
        endif
!      !
!      ! RIGID ROTOR DIPOLE
!      !
!        if ( mod(J_f,2) .eq. 0) then
!        ! j = (-1)**(J_f)
!          j = 1
!        else
!          j = -1
!        endif
!        w3j = wigner3j( dta1%J(k,1), real(K_t,dp), dta1%J(k,2), &
!                        0.0_dp, 0.0_dp, 0.0_dp)
!        dta1%Drigrotor(k) = j*dsqrt(2.d0*dta1%J(k,2) + 1.)*w3j
       !
    ENDDO
! 
  END SUBROUTINE DipCAL
!--------------------------------------------------------------------------------------------------------------------
  double precision function pureHund(caseHund, J_i, N_i, J_f, N_f)
!--------------------------------------------------------------------------------------------------------------------
! " pureHund " = Pure Hund's case
!
! Detailed description:
! ---------------------
! In rotational-vibrational and electronic spectroscopy of DIATOMIC MOLECULES, 
! Hund's coupling cases are idealized cases where specific terms appearing in 
! the molecular Hamiltonian and involving couplings between angular momenta 
! are assumed to dominate over all other terms. There are five cases, 
! traditionally notated with the letters (a) through (e). 
! Most diatomic molecules are somewhere between the idealized cases (a) and (b).
!
! Variables:
! J_i	  = Resultant total angular momentum quantum number L + S, 
!         excluding nuclear spins.(Subindex "i" = initial state) 
! J_f	  = Resultant total angular momentum quantum number L + S,
!         excluding nuclear spins. (Subindex "f" =final   state)
! N_i/f	= the total angular momentum minus the electron spin J - S (i/f)
!         Rotational angular momentum quantum number, 
!         excluding electron and nuclear spins, 
!         in the case where electron spin is present.
! auxiliary var:
! L	= the electronic orbital angular momentum
! S	= the electronic spin angular momentum
!
! Hund's case	Electrostatic	Spin-orbit	    Rotational
! -----------------------------------------------------------
! (a)		    strong	    	intermediate	   weak
! (b)		    strong		    weak		         intermediate
! (c)		    intermediate	strong		       weak
! (d)		    intermediate	weak		         strong
! (e)		    weak		      intermediate	   strong
!				    strong		    intermediate
!
      IMPLICIT NONE
      integer*8, parameter         :: dp = kind(1.0D0) !double precission
      integer (kind=8), intent(in) :: N_i, N_f
      double precision, intent(in) :: J_i, J_f
      character       , intent(in) :: caseHund
      double precision             :: cte1
! -----------------
    if (caseHund .eq. 'a') then
        !Case (a)
        pureHund = 1.0_dp
    else if (caseHund .eq. 'b') then
        !Case (b)
        if ( (N_f .eq. N_i+1 ) .and. ( J_f .eq. J_i+1.0_dp ) ) then
         cte1 = dsqrt(J_i)
      else if ( (N_f .eq. N_i+1 ) .and. ( J_f .eq. J_i ) ) then
         cte1 = (-1.0_dp)*dsqrt(J_i + 1.0_dp)
      else if ( (N_f .eq. N_i-1 ) .and. ( J_f .eq. J_i ) ) then
         cte1 = (-1.0_dp)*dsqrt(J_i)
      else if ( (N_f .eq. N_i-1 ) .and. ( J_f .eq. J_i-1.0_dp ) ) then
         cte1 = dsqrt(J_i + 1.0_dp)
      else ! N_f = N_i !Q band
         cte1 = 1.0_dp
      endif
      pureHund = dsqrt(1.0_dp/(2.0_dp*J_i +1.0_dp)) * cte1
    else if (caseHund .eq. 'c') then
        !Case (c)
      pureHund = 1.0_dp
    else if (caseHund .eq. 'd') then
        !Case (d)
      pureHund = 1.0_dp
    else if (caseHund .eq. 'e') then
        !Case (e)
      pureHund = 1.0_dp
    else
        pureHund = 1.0_dp
    endif      
    RETURN
  END function pureHund
!--------------------------------------------------------------------------------------------------------------------
  double precision function Kpart1(Ji_p, Jf, Jf_p, li, lf, AF1)
!--------------------------------------------------------------------------------------------------------------------
! "Kpart1": First part of the OFF-diagonal elements <<k'|W|k>>
        ! where:
        ! <<k'|W|k>> := transition k->k'
! NOTE:
! -----
! c1 = (2*Ji_p+1)*sqrt((2*Jf+1)*(2*Jf_p+1))*(-1)^(li+lf+K_t+1)
!
      use module_common_var
      IMPLICIT NONE
      integer*8, intent(in)         :: li, lf
      double precision, intent(in)  :: Ji_p, Jf, Jf_p, AF1
      double precision              :: cte1, cte2
      ! 
      !
      cte1 = (2.D0*Ji_p+1.d0) * dsqrt((2.D0*Jf+1.d0)*(2.D0*Jf_p+1.d0))
      !
      if ( mod( (li+lf+K_t+1),2 ) .eq. 0) then
        cte2 = 1.d0
      else 
        cte2 = -1.d0
      endif 
      Kpart1 = cte1*cte2*AF1
               
      RETURN
  END function Kpart1
!--------------------------------------------------------------------------------------------------------------------
  double precision function Kpart1_O2(Ji, Ji_p, Jf, Jf_p, Ni, Ni_p, Nf, Nf_p, AF1, econ)
!--------------------------------------------------------------------------------------------------------------------
! "Kpart1_O2": First part of the OFF-diagonal elements <<k'|W|k>> 
! where:
! <<k'|W|k>> := transition k->k'
! NOTE: mode mark the version of the calculation: tran/mak1/mak2
! -----
!
    use module_common_var
    IMPLICIT NONE
    integer*8, intent(in)         :: Ni, Ni_p, Nf, Nf_p
    double precision, intent(in)  :: Ji, Ji_p, Jf, Jf_p, AF1
    type (dta_ERR), intent(inout) :: econ
    double precision              :: cte1, cte2
    !
    Kpart1_O2 = 0.0_dp
    !
    if (mode .eq. 'tran') then
      !! [Tran et al. 2013]
      cte1 = (2.d0*Ni + 1.d0)*(2.d0*Ni_p + 1.d0)*(2.d0*Nf + 1.d0)*(2.d0*Nf_p + 1.d0)
      cte1 = cte1*(2.d0*Jf + 1.d0)*(2.d0*Jf_p + 1.d0)*(2.d0*Ji_p + 1.d0)**2
      if ( mod( int(Ji_p+Ji+K_t),2 ) .eq. 0) then !
        cte2 = 1.d0
      else 
        cte2 = -1.d0
      endif 
      !
    else if (mode .eq. 'mak1') then
      !! [Makarov et al. 2013] 
        cte1 = (2.d0*Ni + 1.d0)*(2.d0*Nf + 1.d0)
        cte2 = dsqrt((2.d0*Ji+1.d0)*(2.d0*Jf+1.d0)*(2.d0*Ji_p+1.d0)*(2.d0*Jf_p+1.d0))
      !
    else if (mode .eq. 'mak2') then
      !! [Makarov et al. 2013] CODE -> not equal to what is written in the paper.
        cte1 = (2.d0*Nf + 1.d0)*(2.d0*Nf_p + 1.d0)
        cte2 = dsqrt((2.d0*Ji+1.d0)*(2.d0*Jf+1.d0)*(2.d0*Ji_p+1.d0)*(2.d0*Jf_p+1.d0))
      !
      
    else
        call errorMode(mode,econ)
    endif
      
        Kpart1_O2 = cte1*cte2*AF1
                    
      RETURN
  END function Kpart1_O2
!--------------------------------------------------------------------------------------------------------------------
  double precision function Kpart2(L, Ji, Ji_p, Jf, Jf_p, li, lf, AF2,econ)
!--------------------------------------------------------------------------------------------------------------------
! "Kpart2": Calculates triangle coefficients for angular momenta.
! NOTE:
! -----
! This version returns 0 if the triangle inequalities are violated.  (RAH)
!
      use module_common_var
      use module_error
      use module_maths
      IMPLICIT NONE
      integer*8, intent(in)        :: L, li, lf
      double precision, intent(in) :: Ji, Ji_p, Jf, Jf_p, AF2
      type (dta_ERR), intent(inout):: econ
      double precision             :: w3j1,w3j2,w6j
      !
      !wigner 3j-Symbol 1:
      ! ( Ji  L  Ji' )
      ! ( li  0  -li )
      w3j1 = wigner3j( Ji_p, real(L, dp), Ji , real(li, dp), 0.0_dp, real(-li, dp) )
      if (abs(w3j1) .lt. TOL) w3j1 = 0.0D0
      ! CASE wig3j0(M,J1,J2,J)
      !                       m (j1 j2 j)  
      !                    (-)  (m  -m 0)  
      !
      !wigner 3j-Symbol 2:
      ! ( Jf  L  Jf' )
      ! ( lf  0  -lf )
      w3j2  = wigner3j( Jf_p, real(L, dp), Jf, real(-lf, dp), 0.0_dp, real(lf, dp) )
      if (abs(w3j2) .lt. TOL) w3j2 = 0.0D0
      !
      !wigner 6j-Symbol:
      ! { Ji   Jf   1 }
      ! { Jf'  Ji'  L }
      w6j   = wigner6j( Ji, Jf, real(K_t,dp), Jf_p, Ji_p, real(L, dp))
      if (abs(w6j) .lt. TOL) w6j = 0.0D0
      !
      !
      if ( .not.( isnan(w3j1  ) .or. isinf(w3j1  )) .and. &
           .not.( isnan(w3j2  ) .or. isinf(w3j2  )) .and. &
           .not.( isnan(w6j   ) .or. isinf(w6j   )) ) then
            Kpart2 = w3j1 * w3j2* w6j * (2.d0*L + 1.d0)/AF2  
      else
            call wignerS_ERROR(w3j1, w3j2, w6j, 0.0_dp, 0.0_dp,econ)     
      endif

      RETURN
  END function Kpart2
!--------------------------------------------------------------------------------------------------------------------
  double precision function Kpart2_O2(L, Ji, Ji_p, Jf, Jf_p, &
                                      Ni, Ni_p, Nf, Nf_p, Si, Sf, &
                                      AF2, econ)
!--------------------------------------------------------------------------------------------------------------------
! "Kpart2_O2": Calculates triangle coefficients for angular momenta.
! NOTE:
! -----
! This version returns 0 if the triangle inequalities are violated.  (RAH)
!
    use module_common_var
    use module_error
    use module_maths
    IMPLICIT NONE
    integer*8, intent(in)        :: L, Ni, Ni_p, Nf, Nf_p
    double precision, intent(in) :: Ji, Ji_p, Jf, Jf_p, Si, Sf, AF2
    type (dta_ERR), intent(inout):: econ
    double precision             :: sK
    double precision             :: w3j1,w3j2,w6j1,w6j2,w6j3 
      !
    Kpart2_O2 = 0.0_dp
      !
      !
    if (mode .eq. 'tran') then
      !! [Tran et al. 2006]
      !
      sK = 1.d0
      !
      !wigner 3j-Symbol 1:
      ! ( Ni'  Ni  L )
      ! ( 0    0   0 )
      w3j1 = wigner3j( real(Ni_p, dp), real(Ni, dp), real(L, dp), 0.0_dp, 0.0_dp, 0.0_dp )  
      !
      !wigner 3j-Symbol 2:
      ! ( Nf'  Nf  L )
      ! ( 0    0   0 )
      w3j2 = wigner3j( real(Nf_p, dp), real(Nf, dp), real(L, dp), 0.0_dp, 0.0_dp, 0.0_dp )
      !  
      !wigner 6j-Symbol 1: (and symmetries)
      ! Paper (Tran 2006)    -------------->     Addapted to subroutine
      ! { L   Ji    Ji'}  =  { Ji' Ji   L  }  =  { Ji' Ni'  Si}  
      ! { Si  Ni'   Ni }     { Ni  Ni'  Si }     { Ni  Ji   L }
      w6j1 = wigner6j(Ji_p, real(Ni_p, dp), Si, real(Ni, dp), Ji,real(L, dp))
      !
      !wigner 6j-Symbol 2:
      ! { L   Jf    Jf'}  =  { Jf' Jf   L  }  =  { Jf' Nf'  Sf}  
      ! { Sf  Nf'   Nf }     { Nf  Nf'  Sf }     { Nf  Jf   L }
      w6j2 = wigner6j(Jf_p, real(Nf_p, dp), Sf, real(Nf, dp), Jf, real(L, dp))
      !
      !wigner 6j-Symbol 3:
      ! { L   Ji   Ji' }  =  { Ji' Ji   L }  =  { Ji' Jf'  1 }  
      ! { 1   Jf'  Jf  }     { Jf  Jf'  1 }     { Jf  Ji   L }
      w6j3 = wigner6j(Ji_p, Jf_p, real(K_t,dp), Jf, Ji, real(L, dp))
      !
      !
    else if (mode .eq. 'mak1') then
      !! [Makarov et al. 2013] 
      !
      if ( mod( int(Ji+Jf+L+1),2 ) .eq. 0) then !corsign=(-1.)**(1+LL+ji+jip)
          sK = 1.d0
      else 
          sK = -1.d0
      endif 
      !
      !
      !wigner 3j-Symbol 1:
      ! ( Ni' Ni  L )
      ! ( 0   0   0 )
      w3j1 = wigner3j( real(Ni_p, dp), real(Ni, dp), real(L, dp), 0.0_dp, 0.0_dp, 0.0_dp ) 
      !w3j2 = 1.0_dp
      w3j2 = w3j1
      !
      !wigner 6j-Symbol 1: (and symmetries)
      ! Paper             -------------->     Addapted to subroutine
      ! { L   Ji  Ji'}  =  { Ji' Ji  L  }  =  { Ji' Ni' Si}  
      ! { Si  Ni' Ni }     { Ni  Ni' Si }     { Ni  Ji  L }
      w6j1 = wigner6j(Ji_p, real(Ni_p, dp), Si, real(Ni, dp), Ji,real(L, dp))
      !
      !wigner 6j-Symbol 2:
      ! { L   Jf   Jf'}  =  { Jf' Jf  L }  =  { Jf' Ni' Sf}  
      ! { Sf  Ni'  Ni }     { Ni  Ni' Sf}     { Ni  Jf  L }
      w6j2 = wigner6j(Jf_p, real(Ni_p, dp), Sf, real(Ni, dp), Jf, real(L, dp))
      !
      !wigner 6j-Symbol 3:
      ! { L   Ji   Ji'}  =  { Ji' Ji  L }  =  { Ji' Jf' 1 }  
      ! { 1   Jf'  Jf }     { Jf  Jf' 1 }     { Jf  Ji  L }
      w6j3 = wigner6j(Ji_p, Jf_p, real(K_t,dp), Jf, Ji, real(L, dp))
      !
    else if (mode .eq. 'mak2') then
      !! [Makarov et al. 2013] CODE -> not equal to what is written in the paper.
      !
      if ( mod( int(Ji+Ji_p+L+1),2 ) .eq. 0) then !corsign=(-1.)**(1+LL+ji+jip)
          sK = 1.d0
      else 
          sK = -1.d0
      endif 
      !
      !
      !wigner 3j-Symbol 1:
      ! ( Nf   Nf' L )
      ! ( 0    0   0 ) = GCM(0,Nf,Nfp,LL)
      w3j1 = wigner3j( real(Nf, dp), real(Nf_p, dp), real(L, dp), 0.0_dp, 0.0_dp, 0.0_dp ) 
      w3j2 = w3j1
      !
      !wigner 6j-Symbol 1: 
      ! { Ji' Nf' Si}  
      ! { Nf  Ji  L }  = sixj(jip,nfp,ji,nf,ll)
      w6j1 = wigner6j(Ji_p, real(Nf_p, dp), Si, real(Nf, dp), Ji,real(L, dp))
      !
      !wigner 6j-Symbol 2:
      ! { Jf' Nf' Sf}  
      ! { Nf  Jf  L }  = sixj(jfp,nfp,jf,nf,ll)
      w6j2 = wigner6j(Jf_p, real(Nf_p, dp), Sf, real(Nf, dp), Jf, real(L, dp))
      !
      !wigner 6j-Symbol 3:
      ! { Ji' Jf' 1 }  
      ! { Jf  Ji  L }  = sixj(jip,jfp,ji,jf,ll)
      w6j3 = wigner6j(Ji_p, Jf_p, real(K_t,dp), Jf, Ji, real(L, dp))
    else
        call errorMode(mode,econ)
    endif

    sK = sK*(2.d0*L + 1.d0)
                                                   
      if ( .not.( isnan(w3j1) .or. isinf(w3j1)) .and. &
           .not.( isnan(w3j2) .or. isinf(w3j2)) .and. &
           .not.( isnan(w6j1) .or. isinf(w6j1)) .and. &
           .not.( isnan(w6j2) .or. isinf(w6j2)) .and. &
           .not.( isnan(w6j3) .or. isinf(w6j3)) ) then
            Kpart2_O2 = w3j1 * w3j2* w6j1 * w6j2 * w6j3 * &
                        sK/AF2 
      else
            call wignerS_ERROR(w3j1, w3j2, w6j1, w6j2, w6j3,econ)
      endif
      
      RETURN
  END function Kpart2_O2
!--------------------------------------------------------------------------------------------------------------------
  double precision function K_jkCalc(nt,j,k,h,nLines,molP,PerM,econ)
!--------------------------------------------------------------------------------------------------------------------
! "K_jkCalc": OFF-diagonal matrix elements (Wk'k ∞ Q_k'k)
!
! Detailed Description:
! ---------------------
! this function gives the theoretic formulation of relaxation matrix elements:
! 			<< k' | W | k >> 
! They have the same functional form as the rotational state-to-state 
! cross sections within a single vibrational state.
! NOTE: this depends on Temperature but through the use of "Module_common_var"
!
! Accessed Files:  None
! --------------
!
! Called Routines: 
! ---------------  
!
! Called Functions: 'wigner3j' ; 'wigner6j' ; 'AFmol_X' ; 'Ql_mol_X'
! ---------------- 
!
! Called By: 'WelCAL' 
! ---------
!
! Double Precision Version
!
! T.Mendaza, last change 17 February 2017
!--------------------------------------------------------------------------------------------------------------------
!
      use module_common_var
      use module_molecSp
      use module_error
      use module_maths
      IMPLICIT none
      integer*8              , intent(in   ) :: nt,j,k,nLines
      type (dta_SDF), pointer, intent(in   ) :: h
      type (dta_MOL)         , intent(in   ) :: molP, PerM
      type (dta_ERR)         , intent(inout) :: econ
      !internal variables:
      integer*8       :: L, incr, step
      integer*8       :: i, iniL,endL
      double precision:: Ji, Ji_p, Jf, Jf_p 
      integer*8       :: pos_L_2, aux_j_L, pos_Ji_2
      integer*8       :: li,lf
      double precision:: suma, cte1, cte2
      double precision:: w3j1,w3j2,w6j,w3jaux
      double precision:: AF1 , AF2, Qaux
      double precision:: Afmol_X, Ql_mol_X, Ql_mol_LLS, Ql_mol_AF_LLS
      double precision:: K1, Kpart1, K2, Kpart2
!
!...co2  -- 
! This expresion is used for "downward" transitions (k"->k') only.
! That means that, given two transitions: 
! k := j_f , l_f , j_i , l_i
! k':= j'_f, l'_f, j'_i, l'_i
! so "downward" refears to that ones which (j'_i < j_i)
! the upwards ones can be obtained from detailed balance: 
!           rho_k'·<<k"|W|k'>> = rho_k"·<<k'|W|k">>
! NOTE: for further info about this formulation [Niro et al. 2004]
!
! STEP = 2 -> Because the quadrupole-quadrupole energy potential is dominant
!      in linear molecules.
  !
    !
      step = 2 
      ! j -> Ji(j) > Ji(k)=Ji'
      Ji   = h%J(j,1)
      Jf   = h%J(j,2)
      ! k-> 
      Ji_p = h%J(k,1)
      Jf_p = h%J(k,2)
      !
      ! li and lf are the vibrational angular momentum of the vibrational
      ! available for Linear molecules in HITRAN; 
      ! (due to its symmetry though)->
      li = h%lv2(1)
      lf = h%lv2(2) 
      ! 
      iniL=max(iabs(int(Ji)-int(Ji_p)),abs(int(Jf)-int(Jf_p)))
      endL=min(iabs(int(Ji)+int(Ji_p)),abs(int(Jf)+int(Jf_p)))
      !iniL=int(max(abs(Ji-Ji_p),abs(Jf-Jf_p)))
      If( mod(iniL,step) .ne. 0 )iniL=iniL+1
      !endL=int(min((Ji+Ji_p),(Jf+Jf_p)))
      ! 
      ! Adiabatic factor 1:
      if (molP%AF_ON) then
        AF1 = AFmol_X(molP, PerM, Ji, step)
      else
        AF1 = 1.0_dp
      endif
      K1 = Kpart1(Ji_p, Jf, Jf_p, li, lf, AF1)
      !
      suma = 0.d0
      !
      !do L = 1,130
      do L = iniL,endL,step
        ! Since the molecules to be considered are linear -> symmetrical 
        ! and the quadrupole - quadrupole interaction for the sistems mol-X
        ! are meant to be dominant the interaction potential is likely 
        ! dominated by EVEN rank contributions, then STEP = 2.
        ! NOTE: 
        ! (In this case) if STEP = odd; 
        ! the basis rates are expected to be small -> Q_mol_X=0.
        !
        ! Calculation of the basis rates "Q_mol-X":
          if (molP % QTy .eq. "TMC") then
            if (molP % LLSty .eq. "Li--AF") then
                Qaux = Ql_mol_AF_LLS(Ji,Ji_p,h%Sig(j),h%Sig(k),molP,PerM,L)
            else !if molP % LLSty = "Linear" or "Model1/2/3/4"
                Qaux = Ql_mol_LLS(Ji,Ji_p,h%Sig(j),h%Sig(k),molP,L,econ)
            endif
          else
            Qaux = Ql_mol_X(molP,L)
          endif
        ! ERASE NEXT LINE
        h%Qlt(nt,L) = Qaux

        if ( abs(Qaux) .gt. TOL ) then
          !
          ! Adiabatic factor 2:
          if (molP%AF_ON) then
            AF2 = AFmol_X(molP, PerM, real(L,dp), step)
          else
            AF2 = 1.0_dp
          endif
          K2 = Kpart2(L, Ji, Ji_p, Jf, Jf_p, li, lf, AF2,econ)
          !
          if ( .not.( isnan(K2) .or. isinf(K2)) .and. &
               .not.( isnan(suma) .or. isinf(suma)) .and. &
               .not.( isnan(AF2) .or. isinf(AF2)) ) then
              suma = suma + K2 * Qaux
          else
            call offdiagEL_IN_ERROR(AF2, Qaux, K2,econ)
          endif
        else 
        endif
      enddo
      !
      K_jkCalc = K1*suma 
      ! ERASE NEXT LINE
      ! h%Qlt(nt,L) = K_jkCalc

      if ( isnan( K_jkCalc ) .or. isinf( K_jkCalc  ) ) then
        call offdiagEL_ERROR(K1, AF1, suma,econ)
      endif
    RETURN
  END function K_jkCalc
!--------------------------------------------------------------------------------------------------------------------
  double precision function K_jkO2(j,k,h,nLines,molP,PerM,econ)
!--------------------------------------------------------------------------------------------------------------------
! "K_jkO2": OFF-diagonal matrix elements for O2 systems(Wk'k ∞ Q_k'k)
!
! Detailed Description:
! ---------------------
! Same as K_jkcalc
!
! Accessed Files:  None
! --------------
!
! Called Routines: 
! ---------------  
!
! Called Functions: 'wigner3j' ; 'wigner6j' ; 'AFmol_X' ; 'Ql_mol_X'
! ---------------- 
!
! Called By: 'WelCAL' 
! ---------
!
! Double Precision Version
!
! T.Mendaza, last change 26 January 2017
!--------------------------------------------------------------------------------------------------------------------
!
      use module_common_var
      use module_error
      use module_molecSp
      use module_maths
      IMPLICIT none
      integer*8              , intent(in   ) :: j,k,nLines
      type (dta_SDF), pointer, intent(in   ) :: h
      type (dta_MOL)         , intent(in   ) :: molP, PerM
      type (dta_ERR)         , intent(inout) :: econ
      !internal variables:
      integer*8             :: L, incr, step
      integer*8             :: i, iniL,endL
      double precision      :: Ji, Ji_p, Jf, Jf_p 
      integer*8             :: Ni, Ni_p, Nf, Nf_p 
      real*8                :: Si, Sf,n
      integer*8             :: pos_L_2, aux_j_L, pos_Ji_2
      integer*8             :: li,lf
      real*8                :: suma, cte1, cte2
      double precision      :: w3j1,w3j2
      double precision      :: w6j1, w6j2, w6j3
      double precision      :: AF1,AF2, Qaux
      double precision      :: Afmol_X, Ql_mol_X, Ql_mol_LLS, Ql_mol_AF_LLS
      double precision      :: K1, Kpart1_O2, K2, Kpart2_O2
!
!...O2  -- 
! This expresion is used for "downward" transitions (k"->k') only.
! That means that, given two transitions: 
! k := jf , Nf , ji , Ni, Si
! k':= jf', Nf', ji', Ni'
! so "downward" refears to that ones which (j'_i < j_i)
! the upwards ones can be obtained from detailed balance: 
!           rho_k'·<<k"|W|k'>> = rho_k"·<<k'|W|k">>
! NOTE: for further info about this formulation [Niro et al. 2004]
!
! STEP = 2 -> Because the quadrupole-quadrupole energy potential is dominant
!             in linear molecules.
  !
    !
      step = 2 
      ! j -> Ji(j) > Ji(k)=Ji'
      Ji   = h%J(j,1); Ni = h%N(j,1); Si = h%espin(j,1)
      Jf   = h%J(j,2); Nf = h%N(j,2); Sf = h%espin(j,2)
      ! k-> 
      Ji_p = h%J(k,1); Ni_p = h%N(k,1)
      Jf_p = h%J(k,2); Nf_p = h%N(k,2)
      !
      Si   = h%espin(j,1);
      Sf   = h%espin(j,2);
      !
      ! 
      iniL=int(max(abs(Ni-Ni_p),abs(Nf-Nf_p)))
      If( mod(iniL,step) .ne. 0 )iniL=iniL+1
      endL=int(min((Ni+Ni_p),(Nf+Nf_p)))
      ! 
      ! Adiabatic factor 1:
      if (molP%AF_ON) then
        AF1 = AFmol_X(molP, PerM, real(Ni,dp), step)
      else
        AF1 = 1.0_dp
      endif
      K1 = Kpart1_O2(Ji, Ji_p, Jf, Jf_p, Ni, Ni_p, Nf, Nf_p, AF1, econ)
      !
      suma = 0.d0
      !
      !do L = 1,130
      do L = iniL,endL,step
        ! Since the molecules to be considered are linear -> symmetrical 
        ! and the quadrupole - quadrupole interaction for the sistems mol-X
        ! are meant to be dominant the interaction potential is likely 
        ! dominated by EVEN rank contributions, then STEP = 2
        !
        ! Calculation of the basis rates "Q_mol-X":
        if (molP % QTy .eq. "TMC") then
            if (molP % LLSty .eq. "Li--AF") then
                !Qaux = Ql_mol_AF_LLS(Ji,Ji_p,h%Sig(j),h%Sig(k),molP,PerM,L)
                Qaux = Ql_mol_AF_LLS(real(Ni,dp),real(Ni_p,dp),h%Sig(j),h%Sig(k),molP,PerM,L)
            else !if molP % LLSty = "Linear" or "Model1/2/3/4"
                !Qaux = Ql_mol_LLS(Ji,Ji_p,h%Sig(j),h%Sig(k),molP,L)
                Qaux = Ql_mol_LLS(real(Ni,dp),real(Ni_p,dp),h%Sig(j),h%Sig(k),molP,L,econ)
            endif
        else
            Qaux = Ql_mol_X(molP,L)
        endif
        !
        if ( abs(Qaux) .gt. TOL ) then
          !
          ! Adiabatic factor 2:
          if (molP%AF_ON) then
            AF2 = AFmol_X(molP, PerM, real(L,dp), step)
          else
            AF2 = 1.0_dp
          endif
          K2 = Kpart2_O2(L, Ji, Ji_p, Jf, Jf_p, &
                         Ni, Ni_p, Nf, Nf_p, Si, Sf, &
                         AF2,econ)
          !
          if ( .not.( isnan(K2  ) .or. isinf(K2  )) .and. &
               .not.( isnan(suma) .or. isinf(suma)) .and. &
               .not.( isnan(AF2 ) .or. isinf(AF2 )) ) then
              suma = suma + K2 * Qaux
          else
            call offdiagEL_IN_ERROR(AF2, Qaux, K2,econ)
          endif
        !else 
        endif
      enddo
      !
      K_jkO2 = K1*suma
      if ( isnan( K_jkO2 ) .or. isinf( K_jkO2  ) ) then
        call offdiagEL_ERROR(K1, AF1, suma,econ)
      endif
    RETURN
  END function K_jkO2
!--------------------------------------------------------------------------------------------------------------------
  double precision function Ql_mol_X(molP,L)
!--------------------------------------------------------------------------------------------------------------------
! " Ql_mol_X": Basis rate.
! 
! Detailed description:
! ---------------------
! Basis rates: semi-empirical function.
!
      use module_common_var
      use module_molecSp
      use module_maths
      IMPLICIT none
      ! a1, a2, a3, dc were declared in module_common_var:
      integer*8        :: L
      real*8           :: E_l, T
      type(dta_MOL)    :: molP
!
! This expresion is used for "downward" transitions (k->k') only.
!
  !
    T = molP % Temp
    !
    if (L .lt. TOL) then
      Ql_mol_X=0.0_dp
    else
    !reduce rotational energy from the level L * 1/B0
        E_l = real((L*L + L), dp) 
    !
        ! O2: [Tran et al., 2006]
        !      DOES NOT WORK!
        !Ql_mol_X = molP%a1* &
        !       ( (  E_l )**(-molP%a2) )* &
        !       ( ( T/T0 )**(-molP%a3) )
        ! 
        !
        !CO2: [Niro et al., 2004]
        ! O2: [Makarov et al., 2013]      
        !Ql_mol_X = a1*( ( E_l )**(-a2) )*( dexp(-a3*c*hplank*molP%B0*E_l/(T*kb))  )
        Ql_mol_X = molP%a1* &
               ( ( E_l )**(-molP%a2) )* &
               ( dexp(-molP%a3*c2*molP%B0*E_l/T)  )
        !
    endif
    RETURN
  END function Ql_mol_X
!--------------------------------------------------------------------------------------------------------------------
  double precision function Ql_mol_LLS(J1,J2,v1,v2,molP,L,econ)
!--------------------------------------------------------------------------------------------------------------------
! " Ql_mol_LLS": Basis rate
! 
! Detailed description:
! ---------------------
! Basis rates: semi-empirical function.
!
      use module_common_var
      use module_molecSp
      use module_maths
      IMPLICIT none
      !Inputs
      integer*8        :: L
      type(dta_MOL)    :: molP
      type(dta_ERR)    :: econ
      real*8           :: J1,J2,v1,v2
      !internal
      real*8           :: E_l, T,Jaux, delta, coeffC
!
! This expresion is used for "downward" transitions (k->k') only.
!
  !
    T = molP % Temp
    Jaux = J1
    if (Jaux .lt. TOL) Jaux = 0.5_dp
    !
    if (L .lt. TOL) then
      Ql_mol_LLS=0.0_dp
    else
    !reduce rotational energy from the level L * 1/B0
      E_l = real((L*L + L), dp) 
    !
    ! 
    ![Mendaza et al., 2017]
    !
      if (molP%M .eq. 7) then! mol == O2
        !
        !Based in [Tran et al., 2006]
        !Ql_mol_LLS = A1 - A2*Ln( E_l ) - A3*( Ln(T) - Ln(T0) )
        !Ql_mol_LLS = molP%a1 - &
        !        ( (molP%a2)*log(E_l) ) - &
        !        ( molP%a3* ( log(T) - log(T0) ) )
        !
        !Based in [Makarov et al., 2013]:
        !Ql_mol_LLS = A1*Ln(2*L+1) - A2*Ln( E_l ) - A3*( c2*molP%B0*E_l/T )
        !
        Ql_mol_LLS = (molP%a1) - &
                ( (molP%a2)*log(E_l) ) - &
                ( molP%a3*c2*molP%B0*E_l/T ) + 1.0_dp
      else
        ! Bsed in [Niro et al. 2004]
        !Ql_mol_LLS = A1 - A2*Ln( E_l ) - A3*( c*hplank*B0*E_l/(T*kb) ) 
        !
        ![Mendaza et al., 2017]
        if (molP % LLSty .eq. "Linear") then
        ! 1) Linear
            Ql_mol_LLS = (molP%a1) - &
                ( (molP%a2)*log(E_l) ) - &
                ( molP%a3*c2*molP%B0*E_l/T ) + 1.0_dp
        else if (molP % LLSty .eq. "Model1") then 
        ! 2) correction 1
            Ql_mol_LLS = (molP%a1)/((dabs(v1-v2)/molP%B0)**(J2/Jaux)) - &
                ( (molP%a2)*log(E_l) ) - &
                ( molP%a3*c2*molP%B0*E_l/T ) + 1.0_dp
        else if (molP % LLSty .eq. "Model2") then
        ! 2) correction 2
            Ql_mol_LLS = (molP%a1)*T/(c2*(dabs(v1-v2)**(J2/Jaux))) - &
                ( (molP%a2)*log(E_l) ) - &
                ( molP%a3*c2*molP%B0*E_l/T ) + 1.0_dp
        else if (molP % LLSty .eq. "Model3") then
        ! 2) correction 3
            Ql_mol_LLS = (exp(-(c2/Jaux)*dabs(molP%v0-v1)/T))*(molP%a1) - &
                ( (molP%a2)*log(E_l) ) - &
                ( molP%a3*c2*molP%B0*E_l/T ) + 1.0_dp  
        else if (molP % LLSty .eq. "Model4") then
        ! 2) Correction 4
            delta = molP%v0-v1
            if (delta .lt. TOL) delta = 0.5
            Ql_mol_LLS = (molP%a1)/(exp(-molP%B0/(dabs(delta)**(1/Jaux)))) - &
                ( (molP%a2)*log(E_l) ) - &
                ( molP%a3*c2*molP%B0*E_l/T ) + 1.0_dp 
        else
            call errorLLStyType(molP % LLSty,econ)
        endif
        !
      endif
      !
    endif

    RETURN
  END function Ql_mol_LLS
!--------------------------------------------------------------------------------------------------------------------
  double precision function Ql_mol_AF_LLS(J1,J2,v1,v2,molP,PerM,L)
!--------------------------------------------------------------------------------------------------------------------
! " Ql_mol_AF_LLS": Basis rate
! 
! Detailed description:
! ---------------------
! Basis rates: semi-empirical function.
!
      use module_common_var
      use module_molecSp
      use module_maths
      IMPLICIT none
      !Inputs
      integer*8        :: L
      type(dta_MOL)    :: molP, PerM
      real*8           :: J1,J2,v1,v2
      !internal
      real*8           :: E_l, T,c_AF, coeffC
!
! This expresion is used for "downward" transitions (k->k') only.
!
  !
    T = molP % Temp
    !
    if (L .lt. TOL) then
      Ql_mol_AF_LLS=0.0_dp
    else
    !reduce rotational energy from the level L * 1/B0
      E_l = real((L*L + L), dp) 
    !
    ! 
    ![Mendaza et al., 2017]
    ! Correction 5: ---> 5
        coeffC = c_AF(molP, PerM, J1, real(L,dp), 2)
        Ql_mol_AF_LLS = (molP%a1) + &
                ( (molP%a2)*log(E_l) ) + &
                ( molP%a3*c2*molP%B0*E_l/T ) + &
                ( (molP%a4)*2*coeffC ) + &
                ( molP%a5*2*coeffC*log(E_l) ) + &
                ( molP%a6*2*coeffC*c2*molP%B0*E_l/T ) + &
                ( molP%a8*coeffC**2 ) + &
                ( (molP%a7)*2*coeffC ) + &
                ( molP%a10*(c2*molP%B0*E_l/T)*coeffC**2 ) + &
                ( (molP%a9)*log(E_l)*coeffC**2 ) + &
                ( molP%a11*coeffC**2 ) + &
                1.0_dp
    endif

    RETURN
  END function Ql_mol_AF_LLS
!--------------------------------------------------------------------------------------------------------------------
  double precision function AFmol_X(molP, PerM, L, step)
!--------------------------------------------------------------------------------------------------------------------
! "AFmol_X": Adiabaticity-Factor for the chemical reaction CH4-X 
!
! Description:
! The original expresion for this correction factor (a.k.a. Adiabaticity-Factor) 
! is given by Rodriguez et al. 1999.
! AF(J,T) = { 1 + 1/24 * [ (w_j_j2 * dc)/(v_mol-X) ]^2  }^-2
! UNITLESS
!
! where:
! w_j_jstep := frequency spacing 2*pi*c* |E(J) - E(J-step)|; E = Ev + Erot(J)
!           where the Vibrational energy (Ev) is the same for both J, and J-step
!           but not for the rotational part (Erot).
!           This leads to the reduced rotational energy Erot = B0*( J(J+1) )
!           |Erot(J) - Erot(J-step)| = step*B0[2J+1 - step]
! dc : scaling length (fitted parameter)
! v_mol-X := relative mean speed of the system. This value is obtained from the
!             Maxwell–Boltzmann distribution "most expected value":
!             <v> = ( (8*kb*T*mu)/pi)^1/2
!             where: kb is the Boltzmann cte. and mu is the reduced mass
!                    given by: mu = (1/m1 + 1/m2)
!
      use module_common_var
      IMPLICIT none
      integer*8    , intent(in)  :: step
      real*8       , intent(in)  :: L
      type(dta_MOL), intent(in)  :: molP, PerM
      double precision,Parameter :: cAF = 0.0006983_dp ! Adiabaticy Factor cte
                                 ! General constant non-molecule dependent included in the 
                                 ! Adiabaticy Factor. It is built as follows:
                                 ! * Velocity term 
                                 !   cte1 = 1./(Na·8·kb / (1E-03·pi))
                                 ! * Freq. spacing term 
                                 !   cte2 = (2·pi·c)^2 
                                 ! * AF correction form self-cte: 
                                 !   1/24
                                 ! * Unit conversion for dc (Å) to m 
                                 !                       dc = 2.2 Å = 2.2E-10 m
                                 !   (to cancel velocity units m/s):
                                 !   1 Å = 1E-10 m 
                                 !   NOTE: dc is squared in the formula
                                 !   cAF = (1./24.)*cte1*cte2*1E-20
                                 ! 
      double precision           :: v_mol_X, w_j_jstep
      double precision           :: mu, T
      !-----------------------------------------
            T = molP % Temp
      !-----------------------------------------
      !
      !
            if (L .lt. TOL) then
              AFmol_X = 1.0_dp
            else
              mu = ( 1.0_dp/molP%mms + 1.0_dp/PerM%mms ) ! 1/(g/mol)
            ! Velocity Term:
            ! Since it is a velocity its units are m/s. 
            ! Example: N2 at T=300K 
            ! N2_M = 28.006148 g/mol -> 1E-03 · 28.006148 / Na Kg
            ! mu = 1/(1E-03·N2_M/Na) Kg^-1
            ! kb is in J/K where J = Kg·m2/s2
            ! vp = sqrt(2·kb·T·mu) ~ 420 m/s = 420·1E02 cm/s
            ! vm = sqrt(8·kb·T·mu/pi) ~ 473 m/s = 473 · 1E02 cm/s
            ! cte1 = dsqrt(Na*8*kb / (1E-03*pi)) (J/(K·mol))^1/2 = (Kg·m2/(K·s2·mol))^1/2
            ! vm = cte1·sqrt(T/N2_M) ~ 473 m/s
            ! 
              v_mol_X  = 1.0_dp/T*mu ! ----> squared and inverted velocity term 
            !write(*,*) "v=", v_mol_X
            !
            ! Frequency spacing:
              w_j_jstep  = (molP%B0*( L + L + 1 - step )*step)**2 ! ---> squared freq. spacing
            !write(*,*) "w_j_j-2=", w_j_jstep
            ! Since we have to give this quantity in (s-1) then
            ! we need to multiply the: 
            ! |Erot(J) - Erot(J-step)| = step*B0[2J+1 - step]
            ! by cte2 = 2*pi*c <- for the units shake
            ! Regarding that B0 is in cm-1 then c must be in cm/s (which is the cte by default in this program)
            ! NOTE: careful if you change the units of the constants given in "module_common_var"
            !       lots of calculations in this program are fitted to the units given.
            ! w_j_jstep  = cte2*(molP%B0*( L + L + 1 - step )*step)**2 
            ! UNITS: s-1
            !
              AFmol_X = 1.0_dp/( 1.0_dp + cAF*( v_mol_X )*( w_j_jstep ) * (molP%dc * molP%dc) )**2
            !
            endif
      RETURN
  END function AFmol_X
!--------------------------------------------------------------------------------------------------------------------
  double precision function c_AF(molP, PerM, J, L, step)
!--------------------------------------------------------------------------------------------------------------------
! "c_AF": LLS factor
!
!
      use module_common_var
      IMPLICIT none
      integer*8    , intent(in)  :: step
      real*8       , intent(in)  :: J,L
      type(dta_MOL), intent(in)  :: molP, PerM
      double precision,Parameter :: cAF = 0.0006983_dp ! Adiabaticy Factor cte
                                 ! General constant non-molecule dependent included in the 
                                 ! Adiabaticy Factor. It is built as follows:
                                 ! * Velocity term 
                                 !   cte1 = 1./(Na·8·kb / (1E-03·pi))
                                 ! * Freq. spacing term 
                                 !   cte2 = (2·pi·c)^2 
                                 ! * AF correction form self-cte: 
                                 !   1/24
                                 ! * Unit conversion for dc (Å) to m 
                                 !                       dc = 2.2 Å = 2.2E-10 m
                                 !   (to cancel velocity units m/s):
                                 !   1 Å = 1E-10 m 
                                 !   NOTE: dc is squared in the formula
                                 !   cAF = (1./24.)*cte1*cte2*1E-20
                                 ! 
      double precision           :: v_mol_X, w_j_jstep, w_l_lstep
      double precision           :: mu, T
      !-----------------------------------------
            T = molP % Temp
      !-----------------------------------------
      !
      !
            if (L .lt. TOL) then
              c_AF = 1.0_dp
            else
              mu = ( 1.0_dp/molP%mms + 1.0_dp/PerM%mms ) ! 1/(g/mol)
            ! Velocity Term:
              v_mol_X  = 1.0_dp/T*mu ! ----> squared and inverted velocity term 
            !
            ! Frequency spacing:
              w_j_jstep  = (molP%B0*( J + J + 1 - step )*step)**2 ! ---> squared freq. spacing
              w_l_lstep  = (molP%B0*( L + L + 1 - step )*step)**2 ! ---> squared freq. spacing
            !
              c_AF = cAF*(( ( v_mol_X )*( w_l_lstep ) ) - ( ( v_mol_X )*( w_j_jstep ) ))
            !
            endif
      RETURN
  END function c_AF
!--------------------------------------------------------------------------------------------------------------------
  subroutine sumRule(nLines,indexS,dipole,Wmat,dfact,econ)
!--------------------------------------------------------------------------------------------------------------------
        use module_common_var
        use module_error
        implicit none
        integer (kind=8), intent(in)    :: nLines
        integer (kind=8), intent(in)    :: indexS(nLines)
        Double Precision, intent(in)    :: dfact
        Double Precision, intent(in)    :: dipole(nLines)
        Double Precision, intent(in)    :: Wmat(nLines,nLines)
        type (dta_ERR),   intent(inout) :: econ
        !local variables
        integer (kind=8)                :: i,j
        Double Precision                :: Saux, DipAux, Wij, Wii
        Double Precision, Parameter     :: TOLe2 = 1.0E-02!
        logical                         :: isinf, isnan
        logical                         :: testOK

!--------------------------------------------------------------------
        testOK = .true.
        DO i = 1, nLines
            Saux = 0.0D0
            DO j = 1, nLines
                if (j .ne. i) then
                    DipAux = dipole( indexS(j) )/dipole( indexS(i) )
                    Wij = Wmat(i,j)
                    if ( isnan(   Wij) .or. isinf(   Wij) ) stop
                    if ( isnan(DipAux) .or. isinf(DipAux) ) stop
                    Saux = Saux + DipAux*Wij
                else
                    Wii = Wmat(i,i)*dfact ! dfact = diagonal factor
                endif
            ENDDO
            if ( (abs(Wii+Saux) .gt. TOLe2) .and. (i .ne. nLines) ) then
                testOK = .false.
                if (econ % e(1) .ge. 1) then
                    write(*,'(a4,i3,a30)') "row#",i,"NOT meet SUM-RULE"
                    write(*,*) "half-Width", Wii , "?= SUMij{Wmat}", -Saux
                endif
            else
            ! Uncomment if you would also like to see 
            ! the lines that meet the requirement
            !
            !   if (econ % e(1) .ge. 1) then
            !       write(*,'(a4,i3,a30)') "row#",i,"meet SUM-RULE"
            !       write(*,*) "half-Width" Wii , "?= SUMij{Wmat}", -Saux
            !   endif
            endif
        ENDDO
      
        !write(*,*) "SUM-RULE TEST FINISHED"
        if ( .not.(testOK) ) then
            call SumRuleError(econ)
            econ % solu = 0
        else
            if (econ % e(1) .ge. 1) then
                write(*,*) "sumRule: The calculation correctly verifies the sum rule!"
            endif
            econ % solu = 1
        endif

  end subroutine sumRule
!--------------------------------------------------------------------------------------------------------------------