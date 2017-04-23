MODULE module_common_var

    IMPLICIT NONE
    save
    !
    !**************
    ! Kind helpers
    !**************
    integer, parameter  :: Int_8  = selected_int_kind(8)
    integer, parameter  :: Int_10 = selected_int_kind(10)
    integer, parameter  :: dp = kind(1.0D0) !double precission
    integer, parameter  :: sp = kind(1.0)   !simple precission
    integer, parameter  :: r8 = kind(1.0D0) !real*8
    integer, parameter  :: r4 = kind(1.0)   !real*4
    !
    !
    !*************************
    ! Population moment Transition
    !*************************
    ! ptype   : this parameter marks the result Population calculation 
    !           in other words:
    !           ('hit') Acording to HITRAN96 
    !           ('tra') Acording to Hartmann et. al 2006 
    character*3, parameter :: ptype = 'tra'
    !
    !*************************
    ! Dipole moment Transition
    !*************************
    !integer, parameter :: K_t = 0 !Isotropic-spectroscopy
    integer, parameter  :: K_t = 1 !IR-spectroscopy
    !integer, parameter :: K_t = 2 !Raman-spectroscopy
    ! Diatomic-molecule (or linear):
    ! HUND's CASE
    character, parameter   :: caseHund = 'b'
    ! Type of Dipole Calculation
    character, parameter   :: tdcal = 'S' 
    ! Transition Moment Type
    character*3, parameter :: tmt = 'edt' 
    !     1) Electric dipole transition (edt)
    !     2) magnetic dipole transition (mdt)
    !     3) electric-quadrupole transitions (eqt)
    character*4, parameter :: mode =  'mak1'
    !     1) [Tran et al. 2006] ECS matrix element           = 'tran'
    !     2) [Makarov et al. 2013] ECS matrix element        = 'mak1'
    !     3) [Makarov et al. 2013]'s CODE ECS matrix element = 'mak2'
    !
    double precision, parameter :: TOL= 1.0000000E-90! tolerance level
    !***********
    ! Constants
    !***********
    !
    ! Base Temperature:
    double precision, parameter :: T0= 296d0!K
    !
    ! Ph.ctes : 
    real(dp), parameter :: Pi      = 3.141592654d0 ! Pi number
    !real(dp) , parameter :: pi      = 3.1415926535897932384626433832795029_r8
    !real(dp) , parameter :: invpi   = 1./Pi
    !real(dp) , parameter :: twopi   = 2.*Pi
    !real(dp) , parameter :: halfpi  = Pi/2.
    real(dp), parameter :: c2      = 1.4387770 ! second radiation constant (cm·K)
                                  != hplank*c/kb
    real(dp), parameter :: c2_mK   = c2*1E-02! second radiation constant (m·K)
                                  != hplank*c/kb
    real(dp), parameter :: kb      = 1.380658E-23 !Bolzman constant ( J/K) 
                               !or = 8.617385E-05 !                 (eV/K)
    real(dp), parameter :: Rg      = 8.3144598 !Gas constant ( J·mol-1·K-1)                            
    real(dp), parameter :: hplank  = 6.6260755E-34!Plank constant ( J·s)
                               !or = 4.1356692E-15!               (eV·s)
    real(dp), parameter :: hp_cgs  = 6.6260755E-27!Plank constant (erg·s)
    real(dp), parameter :: c_ms    = 2.99792458E+08! speed of light ( m/s)
    real(dp), parameter :: c       = 2.99792458E+10! speed of light (cm/s)
    real(dp), parameter :: Na      = 6.022E23  ! Avogadro constant (molecules·mol-1)
    real(dp), parameter :: v_permit= 1.0/(4.0*Pi)!vacuum permittivity (cgs units)
                                     ! ε0 = 8.854187817 * 1E−12 F·m-1 (farads per metre)
                                     ! is the measure of the resistance that is encountered
                                     ! when forming an electric field in a medium.
    real(dp), parameter :: v_permea= 4*Pi/c**2 !vacuum permeability (N·A-2)
                                     ! 4*Pi*1.0E-07 !vacuum permeability (N·A-2)
                                     ! µ0 = 4π×10−7 N / A2
                                     ! the value of magnetic permeability in a classical vacuum
    !
    !**************
    ! IN/OUT Files
    !**************
    !
    ! Paths
    ! -----
    character (64) , parameter :: out_file_path   = "Output_files/"
    !
    ! Units of IN/OUT files
    ! ---------------------
    integer*8, parameter  :: u1 = 10 ! Unit input file#1  
    integer*8, parameter  :: u2 = 11 ! Unit input file#2
    integer*8, parameter  :: u3 = 20 ! Unit output file
    integer*8, parameter  :: u4 = 12 ! Unit input  file
    integer*8, parameter  :: u5 = 13 ! Unit input  file
    integer*8, parameter  :: u6 = 14 ! Unit input  file
    integer*8, parameter  :: u7 = 15
    integer*8, parameter  :: u8 = 16
    integer*8, parameter  :: u9 = 17
    !
    !********************************************
    ! Max Number of Lines and of matrix elements
    !********************************************
    integer*8, parameter  :: nLmx  = 500 
    integer*8, parameter  :: nMmx  = nLmx*nLmx
    !
    !**************************
    ! Program Type: STRUCTURES 
    !**************************
    type dta_SDF
    ! M        = HITRAN Molecule number
    ! SIG      = Vacuum wavenumber (cm-1)
    ! ISO      = Isotopologue number
    ! STR      = Intensity cm-1/(molecules·cm-2) at 296K
    ! PopuT0   = Populations of the Lower Levels of the Lines
    !	       at 296 K. UNITLESS.
    ! DipoT    = Dipole transition Moments of the Lines. 
    !          UNITS:
    !          Debye^(1/2) = (1E-36 ergs·cm3)^(1/2)
    !                      = (1E-36 (1E-07 J)·cm3)^(1/2)
    ! E        = Lower-state energy (cm-1)
    ! A21      = Einstein A-coefficient (s-1), 
    !            i.e. spontaneous emission (from upper to lower level)
    !           NOTATION (see NOTES bellow).
    ! HWT0     = Air-broadened Half-Widths (at 296 K) of the
    !	       Lines (cm-1/Atm)
    ! BHW      = Temperature Dependence Coefficients of HWT0
    ! SHIFT    = Air pressure-induced line shift
    !
    ! Qupp     = Upper-state ‘local’ quanta".
    !            Rotational information of the upper level.
    !  
    ! Qlow     = Lower-state ‘local’ quanta". 
    !            Rotational information of the lower level.    
    ! 
    ! swei0    = Statistical weight of the upper state
    ! 
    ! swei00   = Statistical weight of the lower state
    !
    ! hitBAND  = the band (Upper quanta) <- (Lower quanta)
    !
    !  NOTES: 
    !  (1) Prime and double primes refer, respectively, 
    !  to upper and lower states, respectively, i.e.
    !       * upper := 2; 0 ; or ' 
    !       * lower := 1;00 ; or " 
    !  (2) J is the quantum number associated with the total
    !  angular momentum excluding nuclear spin; 
    !  (3) F is the quantum number associated with the total
    !  angular momentum including nuclear spin. F is shown in 
    !  A5 FORTRAN format in order to accommodate integer (I5)
    !  or half-integer values (F5.1). 
    !  For group 3, the notations C and a are described in  
    !  [Brown et al. Methane line parameters in HITRAN. JQSRT 2003;82:219–38 ]. 
    !
    ! Swei00   = Statistical weight of the lower state
    !
    integer*8       :: M,iso
    integer*8       :: class
    double Precision,allocatable :: Sig(:), Str(:)
    double Precision,allocatable :: PopuT0(:), PopuT(:)
    double Precision,allocatable :: DipoT(:)!,DipoT0(:)
    double Precision,allocatable :: D0(:)!, Drigrotor(:)
    !double Precision,allocatable:: Y_RosT(:)
    double Precision,allocatable :: E(:), A21(:)
    double Precision,allocatable :: HWT0(:),BHW(:), &
                                    SHIFT(:)
    double Precision,allocatable :: swei0(:),swei00(:)
    character*38                 :: hitBAND
    ! GLOBAL quanta // Vibrational
    ! -------------
    !Class 1: Diatomic molecules
    ! Variables in use:
    ! integer    :: nu1 
    !               nu_j {j=1} is the quantum number associated 
    !               associated with the normal mode j. 
    !
    !Class 2: Diatomic molecules with different electronic levels
    ! Variables in use:
    ! integer    :: nu1
    !               nu_j {j=1} is the quantum number associated 
    !               associated with the normal mode j.
    ! character*1:: X
    !               electronic state of the molecule:
    !               * X = "X" -> ground state
    !               * X = "i" -> where i={1,2,...} any number 
    !                            that shows the exited state#
    !
    !Class 3: Diatomic molecules with doublet-PI electronic state
    ! Variables in use:
    ! integer    :: nu1
    !               nu_j {j=1} is the quantum number associated 
    !               associated with the normal mode j.
    ! character*1:: X 
    !               electronic state of the molecule
    ! character*3:: i_2 
    !               = 1/2 or 3/2 <- N = J +- 1/2 depending on the molecules.
    !
    !Class 4: Linear triatomic
    ! Variables in use:
    ! integer    :: nu1, nu2, nu3
    !               nu_j {j=1,2,3} is the quantum number associated 
    !               associated with the normal mode j.
    ! integer    :: l2
    !               vibrational angular momentum (vam)
    !               quantum number of the vibrational mode {j=2}
    !
    !Class 5: Linear triatomic with large Fermi resonance
    ! Variables in use:
    ! integer    :: nu1, nu2, nu3
    !               nu_j {j=1,2,3} is the quantum number associated 
    !               associated with the normal mode j.
    ! integer    :: l2, r
    !               l2: vibrational angular momentum (vam)
    !               quantum number of the vibrational mode {j=2}
    !               r : Ref[5] of HITRAN04 !!!!!! 
    !
    !Class 6: NON-Linear triatomic
    ! Variables in use:
    ! integer    :: nu1, nu2, nu3
    !               nu_j {j=1,2,3} is the quantum number associated 
    !               associated with the normal mode j.
    !
    !Class 7: Linear tretratomic 
    ! Variables in use:
    ! integer    :: nu1, nu2, nu3, nu4, nu5
    !               nu_j {j=1,...,5} is the quantum number associated 
    !               associated with the normal mode j.
    ! integer    :: lv, sign, r, Schar
    !               lv: absolute value of the sum of the vibrational angular momentum (vam)
    !               quantum number l_j {j=1,...,5}
    !               sign, r, Schar : Ref[7] of HITRAN04 !!!!!! 
    !
    !Class 8: Pyramidal tretratomic 
    ! Variables in use:
    ! integer    :: nu1, nu2, nu3, nu4
    !               nu_j {j=1,...,4} is the quantum number associated 
    !               associated with the normal mode j.
    ! integer    :: Sint
    !               Symmetry of the level 
    !               (only for NH3, for PH3 is blank)
    !
    !Class 9: NON-Linear tretratomic 
    ! Variables in use:
    ! integer    :: nu1, nu2, nu3, nu4, nu5, nu6
    !               nu_j {j=1,...,6} is the quantum number associated 
    !               associated with the normal mode j.
    !
    !Class 10: Pentatonic or greater polyatomic 
    ! Variables in use:
    ! integer    :: nu1, nu2, nu3, nu4
    !               nu_j {j=1,...,4} is the quantum number associated 
    !               associated with the normal mode j.
    ! character  :: Sv, mi
    !               Sv: symmetry
    !               mi: multiplicity index.    
    !
    ! NOTE:
    !       commented lines in front of some variables 
    !       are there to avoid its duplicity.
    !------------------------------------------------------------
    integer*8       :: lv2(2) ! lower = 1 or 'low'
                                           ! upper = 2 or 'upp'
    !------------------------------------------------------------
    ! LOCAL QUANTA:// Rotational
    ! -------------
    ! 
    !Group 1: Asymmetric rotors
    ! Variables in use:
    ! integer    :: J, Ju, Ka(2), Kc(2), F(2)
    ! character  :: Sym(2)
    !            J: quantum number associated with the 
    !               total angular momentum EXCLUDING nuclear SPIN.
    !            F: quantum number associated with the 
    !               total angular momentum INCLUDING nuclear SPIN. 
    !               (Integer, I5, or Half-Integer, F5.1). 
    !          Sym: is either symmetry (e,f), or +/- for requiered 
    !               symmetry symbols, or (d,q) for magnetic-dipole or
    !               electric-quadrupole transitions (ONLY FOR: O2, N2).
    ! NOTE: 
    ! for NO2 and HO2: J = N Sym 1/2 where Sym = + or -. 
    ! So N is in the place of J for those molecules.
    !
    !Group 2: Diatomic and Linear molecules.
    ! Variables in use:
    ! integer    :: J, F(2)
    ! character  :: Sym, Br ! Ju = J+Br
    !            J: Resultant total angular momentum quantum number, 
    !               excluding nuclear spins.(initial/lower level)
    !           Ju: rotational quantum number J (final/upper level)
    !               (Angular momentum without nuclear spins)
    !            F: quantum number associated with the 
    !               total angular momentum INCLUDING nuclear SPIN. 
    !               (Integer, I5, or Half-Integer, F5.1). 
    !           br: character that marks the branch of the band, 
    !               i.e. Br is the O-, P-, Q-, R-, or S branch
    !               symbol;
    !               The rotational selection rule gives rise to:
    !               *r = R-branch (when ∆J = +1) 
    !               *q = Q-branch (when ∆J =  0)
    !                i.e. J" = 0 and J' = 0, but  ν0≠0 is forbidden 
    !                so pure vibrational transition is not observed 
    !                in most cases
    !               *p = P-branch (when ∆J = -1). 
    !          Sym: takes (d,q) value for magnetic-dipole or
    !               electric-quadrupole transitions (ONLY FOR: O2, N2).
    !
    !
    !Group 3: Spherical rotors
    ! Variables in use:
    ! integer    :: J, Ju, alph(2), F(2)
    ! character*2:: Sr(2)
    !            J: Resultant total angular momentum quantum number, 
    !               excluding nuclear spins.(initial/lower level)
    !           Ju: rotational quantum number J (final/upper level)
    !               (Angular momentum without nuclear spins)
    !            F: quantum number associated with the 
    !               total angular momentum INCLUDING nuclear SPIN. 
    !               (Integer, I5, or Half-Integer, F5.1). 
    !           Sr: (or C in literature) Total rotational symmetry in T_d_(M), 
    !               Gamma = A_1 ,A_2 ,E,F_1 ,F_2 (num# = 1,2,3,4,5)
    !         alph: "α" (or "n" above 3400cm-1) are counting integers.
    !               for levels of the same J and C (== Sr); 
    ! 
    !      Q0:                   Q00:                
    !            J′ C′ α′ F´           J" C" α" F"    
    !         2x i3 a2 i3 a5        2x i3 a2 i3 a5
    !
    !Group 4: Symmetric rotors
    ! Variables in use:
    ! integer    :: J, Ju, K(2), l(2), F(2)
    ! character  :: Sym(2)
    ! character*2:: Sr(2)
    !            J: Resultant total angular momentum quantum number, 
    !               excluding nuclear spins but including electronic spin.
    !               (initial/lower level)
    !           Ju: rotational quantum number J (final/upper level)
    !               (Angular momentum without nuclear spins 
    !                but including electronic spin).
    !            l: Quantum number for vibrational angular momentum (vam).
    !            F: Quantum number associated with the 
    !               total angular momentum INCLUDING both nuclear SPIN 
    !               and electron SPIN. 
    !               (Integer, I5, or Half-Integer, F5.1). 
    !           Sr: (or C in literature) Total rotational symmetry in T_d_(M), 
    !               Gamma = A_+ ,A_- ,E (num# = 1,2,3)
    !          Sym: is either symmetry (e,f), or +/- for requiered 
    !               symmetry symbols.
    !
    !Group 5: Triplet-Sigma ground electronic rotors
    ! Variables in use:
    ! integer    :: J, Ka(2), Kc(2), F(2)
    ! character  :: Sym, Br 
    ! character  :: Br_N    ! Nu = N + Br_N
    !            J: rotational quantum number J (initial/lower level)
    !               (Angular momentum without nuclear spins)
    !           Ju= J + Br 
    !            N: Rotational angular momentum quantum number
    !               (excluding electron and nuclear spins, 
    !               in the case where electron spin is present).
    !               J = N +- 1/2.
    !           Nu= N + Br 
    !           br: character that marks the branch of the band, 
    !               i.e. Br is the O-, P-, Q-, R-, or S branch
    !               symbol;
    !               The rotational selection rule gives rise to:
    !               *r = R-branch (when ∆J = +1) 
    !               *q = Q-branch (when ∆J =  0)
    !                i.e. J" = 0 and J' = 0, but  ν0≠0 is forbidden 
    !                so pure vibrational transition is not observed 
    !                in most cases
    !               *p = P-branch (when ∆J = -1). 
    !          Sym: takes values (d,q) for magnetic-dipole or
    !               electric-quadrupole transitions (ONLY FOR: O2, N2).
    !
    !Group 6: Doublet-PI ground electronic states
    ! Variables in use:
    ! integer    :: F(2), J ! THIS IS F5.1
    ! character  :: Sym, Br ! Ju = J + Br
    !            J: rotational quantum number J (initial/lower level)
    !               (Angular momentum without nuclear spins 
    !                but including electronic spin). 
    !                NOTE: this molecules HAVE electronic spin!
    !           Ju= J + Br 
    !            F: Quantum number associated with the 
    !               total angular momentum INCLUDING both nuclear SPIN 
    !               and electron SPIN. 
    !               (Integer, I5, or Half-Integer, F5.1). 
    !           br: character that marks the branch of the band, 
    !               i.e. Br is the O-, P-, Q-, R-, or S branch
    !               symbol;
    !               The rotational selection rule gives rise to:
    !               *r = R-branch (when ∆J = +1) 
    !               *q = Q-branch (when ∆J =  0)
    !                i.e. J" = 0 and J' = 0, but  ν0≠0 is forbidden 
    !                so pure vibrational transition is not observed 
    !                in most cases
    !               *p = P-branch (when ∆J = -1). 
    !          Sym: takes values (d,q) for magnetic-dipole or
    !               electric-quadrupole transitions (ONLY FOR: O2, N2).
    !
    integer*8,allocatable        :: N(:,:)
    double precision,allocatable :: J(:,:)
    double precision,allocatable :: F(:,:)
    double precision,allocatable :: nspin(:,:), espin(:,:)
    character,allocatable        :: br(:), br_N(:)


    end type dta_SDF
    ! --------------
    type dta_RMF
    ! WT0      = Relaxation Matrix elements (at 296 K) 

        double Precision,allocatable :: WT0(:)

    end type dta_RMF
    ! --------------
    type HITRAN
    ! HITRAN STRUCTURE:-------------
    ! PARAMETER   MEANING			    FORTRAN Type   	Comments or units 
    ! M         = Molecule number    	      	I2		HITRAN chronological assignment
    ! ISO       = Isotopologue number   	    I1		Ordering within a molecule by terrestrial abundance
    ! wno       = Vacuum wavenumber   		    F12.6	cm-1 
    ! S		    = Intensity   			        E10.3	cm-1/(molecules·cm-2) at 296K
    ! A   	    = Einstein A-coefficient   	    E10.3	s-1 
    ! gair      = Air-broadened half-width   	F5.4   	HWHM at 296 K (in cm-1 atm-1)
    ! gself   	= Self-broadened half-width   	F5.4   	HWHM at 296 K (in cm-1 atm-1)
    ! E00   	= Lower-state energy   		    F10.4	cm-1 
    ! nair   	= Temp-depend exp for gair   	F4.2   	unitless, 
    ! Shift_air = Air press-induced line shift  F8.6   	cm-1 atm-1 at 296 K
    ! V0   	    = Upper-state "global" quanta 	A15		Format dependent on the molecule:
    !                                                   *CH4: Hollerith   see Table 3
    ! V00  	    = Lower-state "global" quanta 	A15		( same as V0 )
    ! Q0   	    = Upper-state "local"  quanta   A15		Format dependent on the molecule
    !                                                   *CH4: Brown et al.see Table 8
    ! Q00  	    = Lower-state "local"  quanta   A15		( same as Q0 )
    ! Ierr  	= Uncertainty indices   	    I6		Accuracy for 3/6 critical parameters
    ! 	     							                ðv; S; gair=v; S; gair; gself ; nair; 
    !								                    dairÞ, see Rothman et al. (2013)
    ! Iref   	= Reference indices   		    I12   	References for 3/6 critical parameters
    ! 	     							                ðv; S; gair=v; S; gair; gself ; nair; dairÞ
    ! *   	    = Flag   			            A1		Availability of program and data for the 
    !								                    case of line mixing
    ! g0   	    = Statistical weight of the    	F7.1	UNITLESS
    !		      upper state
    ! g00   	= Statistical weight of the     F7.1	UNITLESS
    !		      lower state

       integer*8              :: M, ISO, Ierr, Iref
       Double Precision       :: wno,S, A, gair, gself
       Double Precision       :: E00, nair, shift_air
       Double Precision       :: g0, g00
       character*15           :: V0, V00, Q0, Q00
       character              :: flag

    end type HITRAN
    ! --------------
    type dta_MOL
    ! MOLECULAR STRUCTURE:-------------
    ! PARAMETER   MEANING               FORTRAN Type    Comments or units 
    ! M         = Molecule number               I2      HITRAN chronological assignment
    ! ISO       = Isotopologue number           I1      Ordering within a molecule by terrestrial abundance
    ! Aco       = AFGL code                     I3      The old Air Force Geophysics Laboratory 
    !                                                   (AFGL) shorthand notation for isotopologues.
    ! IAb       = Iso Abundance                 --      HITRAN Isotopologue abundance.
    ! QT        = Partition function at T or T0 dp      T selected by the user, T0 = 296K
    ! mms       = Molar mass                    dp      Molar Mass [g·mol-1]
    ! Nmcon     = molar concentration at 296K   dp      mol·cm3
    ! B0        = Rotational constant B0        dp      cm-1 
    ! Temp      =  temperature of the Gas       dp      K
    ! Ptot      =  Pressure of the Gas          dp      atm
    !
    ! Further information (QT in different Temperatures) in: 
    ! http://hitran.org/docs/iso-meta/
    ! 
    ! PROGRAM EXTRA INFORMATION:-------------
    ! ai        =  fitting coeff of the         dp      cm-1/No-unit/No-unit
    !              base function      
    ! dc        =  Adiabatic Factor parameter   dp      Å
    ! exi       =  Temperature dependent        dp      No-unit
    !              exponents of ai param.
    ! availableParam = logical variable         lo      0/1 (Just used at "RM_LM_LLS_tmc_arts" )
    !              that checks whether the code 
    !              can access pre-registered 
    !              values of ai
    ! AF_ON      = allows to use the            lo      0/1 (Depends on availability and the freedom that the user needs)
    !              adiabatic factor 
    ! QTy        = "REG" or "TMC"               a3      depends on the called method:
    !                                                   "REG" = RM_LM_tmc_arts
    !                                                   "TMC" = RM_LM_LLS_tmc_arts
    !
    !
    integer*8 :: M
    integer*8 :: iso_m ! == column number.
    integer*8 :: Aco
    character*3            :: QTy
    character*6            :: chmol
    double Precision       :: Temp, Ptot
    Double Precision       :: mms!, IAb
    Double Precision       :: Nmcon, B0
    Double Precision       :: QT, QT0
    Double Precision       :: a1, a2, a3, dc, ex1, ex2
    logical                :: availableParam, AF_ON
    logical                :: error_flag

    end type dta_MOL
    ! --------------
    ! --------------
    type dta_ERR
    ! MOLECULAR STRUCTURE:-------------
    ! PARAMETER   MEANING               FORTRAN Type    Comments or units 
    ! e(1)      = debugging flag                I1      it allows this code to be verbose on screen 
    !                                                   for debugging purposes.
    ! e(2)      = error counter                 I1      It counts the number of errors, if e(2)>=1, 
    !                                                   then the error flag will be send back as 1. 
    !
    integer*8, dimension(2) :: e

    end type dta_ERR
    ! --------------
END MODULE module_common_var
