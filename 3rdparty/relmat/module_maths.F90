!--------------------------------------------------------------------------------------------------------------------
module module_maths
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to 
!
    interface        

        logical function isJb(dta1,i1,i2)
            use module_common_var
            implicit none
            integer (kind=8)      , intent(in) :: i1,i2
            type(dta_SDF), pointer, intent(in) :: dta1  
        end function isJb

        double precision function wigner3j( dJ1,dJ2,dJ3, dM1,dM2,dM3 )
            use module_common_var
#if USE_EXTERNAL_WIGNER
            use ffastwigxj
#endif
            implicit none
            double precision, intent(in)    :: dJ1,dJ2,dJ3
            double precision, intent(in)    :: dM1,dM2,dM3
        end function wigner3j

        double precision function wig3j0(M,J1,J2,J)
            use module_common_var  
            implicit none
            integer*8, intent(in) :: M, J1, J2, J
        end function wig3j0

        double precision function CombiJM(J,M)
            use module_common_var
            implicit none
            integer*8, intent(in) :: J, M
        END function CombiJM
        
        double precision function wigner6j( A,B,uJ3, D,C6,F )
            use module_common_var
#if USE_EXTERNAL_WIGNER
            use ffastwigxj
#endif
            implicit none
            double precision, intent(in)    :: A,B ,uJ3
            double precision, intent(in)    :: D,C6,F
        end function wigner6j
        
        subroutine bubble_index(N, array, indxo, ad, econ)
            use module_common_var
            use module_error
            implicit none
            integer (kind=8)  , intent(in)           :: N
            integer (kind=8)  , intent(inout)        :: indxo(N)
            real*8            , intent(in)           :: array(N)
            character(1)      , intent(in)           :: ad
            type (dta_ERR)    , intent(inout)        :: econ
        end subroutine bubble_index
        
        subroutine ibubble_index(N, array, indxo, ad, econ)
            use module_common_var
            use module_error
            implicit none
            integer (kind=8)  , intent(in)           :: N
            integer (kind=8)  , intent(inout)        :: indxo(N)
            integer (kind=8)  , intent(in)           :: array(N)
            character(1)      , intent(in)           :: ad
            type (dta_ERR)    , intent(inout)        :: econ
        end subroutine ibubble_index
        
        integer*8 function imax(N, iarray)
            implicit none
            integer (kind=8), intent(IN)           :: N
            integer (kind=8), intent(IN)           :: iarray(N)
        end function imax
        
        integer*8 function imin(N, iarray)
            implicit none
            integer (kind=8), intent(IN)           :: N
            integer (kind=8), intent(IN)           :: iarray(N)
        end function imin
        
        double precision function maxf(N, array)
            implicit none
            integer (kind=8), intent(IN)           :: N
            real*8 , intent(IN)                    :: array(N)
        end function maxf
        
        double precision function minf(N, array)
            implicit none
            integer (kind=8), intent(IN)           :: N
            real*8 , intent(IN)                    :: array(N)
        end function minf

        logical function isnan(a)
            implicit none
            double precision, intent(IN) :: a 
        end function isnan

        logical function isinf(a)
            implicit none
            double precision, intent(IN) :: a 
        end function isinf

    end interface

END module module_maths
!--------------------------------------------------------------------------------------------------------------------
! MATHEMATIC FUNCTIONS AND SUBROUTINES ------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
  logical function isJb(dta1,i1,i2)
!--------------------------------------------------------------------------------------------------------------------
!"IsJb": Is "J" bigger 
!
! Detailed description:
! ---------------------
! It says whether the energy level is bigger or not based on its J value
! (GROUND LEVEL! -> subindex == 1)
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! Q00i    : HITRAN's L State Local Quanta for CH4 (Input).
! Ji      : (Input).
! isJb    : it tells the result of the comparision:
!           1 -> greater
!           0 -> smaller
!          -1 -> same level
!
! Accessed Files:  None
! --------------
!
! Called Routines: 'none' 
! ---------------  
!
! Called By: 'WelCAL' (W Elements CALculation)
! ---------
!
! Double Precision Version
!
! T.Mendaza, last change 27 August 2015
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    implicit none
    integer*8             , intent(in) :: i1,i2
    type(dta_SDF), pointer, intent(in) :: dta1
    integer*8 :: J1, J2
!-----------------------------------------
! Q001:
    J1 = dta1%J(i1,1)
! Q002:
    J2 = dta1%J(i2,1)
!
    if (J1 .ge. J2) then
      isJb = .true.
    else
      isJb = .false.
    endif

    Return
  END function isJb
!--------------------------------------------------------------------------------------------------------------------
  double precision function wigner3j( dJ1,dJ2,dJ3, dM1,dM2,dM3 )
!--------------------------------------------------------------------------------------------------------------------
! "wigner3j": wigner symbols 3j.
!
! Detailed explanation: 
! For the calculation of coefficient 3j the 19-Messiah's recurrence  
! is used.
!                           (j1 j2 j3)  
!                           (m1 m2 m3) 
!
! Implementation by: Niro et al. 2004
!
      use module_common_var  
#if USE_EXTERNAL_WIGNER
      use ffastwigxj
#endif
      implicit none
      double precision, intent(in) :: dJ1,dJ2,dJ3,dM1,dM2,dM3
#if USE_EXTERNAL_WIGNER
      integer*4                    :: J1, J2, J3, M1, M2, M3
#else
      !function var.
      integer*8, parameter         :: iMaxi=100
      integer*8                    :: J0, J1, J2, J3, JS, JM, JJ(4)
      integer*8                    :: M0, M1, M2, M3, MM(4), M10
      integer*8                    :: MS,M1S,M20,M2S,MP0,MPS
      integer*8                    :: MM1, MM2, MM3
      integer*8                    :: i, i1, k, m, mp
      integer*8                    :: IM0, IM(4)
      double precision             :: AJ0, AJ1, AJ2, AJ3
      double precision             :: C1, C3, Cs2
      double precision             :: CC(0:iMaxi,0:iMaxi)
      double precision             :: UN, UNM, wig3j0
#endif
!
! Init.
!
#if USE_EXTERNAL_WIGNER
      J1 = 2*int(dJ1,4) ; J2 = 2*int(dJ2,4) ; J3 = 2*int(dJ3,4) 
      M1 = 2*int(dM1,4) ; M2 = 2*int(dM2,4) ; M3 = 2*int(dM3,4) 
      wigner3j = ffw3jja6(J1, J2, J3, M1, M2, M3)
    return
#else
      J1 = int(dJ1,8) ; J2 = int(dJ2,8) ; J3 = int(dJ3,8) 
      M1 = int(dM1,8) ; M2 = int(dM2,8) ; M3 = int(dM3,8) 

      JJ(1) = J1 ; JJ(2) = J2 ; JJ(3) = J3 
      MM(1) = M1 ; MM(2) = M2 ; MM(3) = M3 

      wigner3j=0.0_dp
      MS=0
      DO 1 i=1,3
        MS=MS+MM(i)
        IM(i)=iabs(MM(i))
        IF( .NOT.((JJ(i) .GE. 0) .AND. (IM(i) .LE. JJ(i)) ) )RETURN
1     CONTINUE     
!
! Other selection rules
!
      IF( MS .NE. 0 )RETURN
      J0=iabs(J1-J2)
      JM=J1+J2
      IF( (J3.GT.JM) .OR. (J3.LT.J0) )RETURN
!
! Parameter test 
!
      IM0=MIN0(IM(1),IM(2),IM(3))
      IF(IM0.GT.iMaxi)THEN
        write(8,100)IM0,iMaxi
100    FORMAT(' min(abs(m))=',I2,' > iMaxi=',I2)
        RETURN
      ENDIF
!
! Circular Permutations bringing the smallest abs (m) to the right
!
      IF(IM(3).EQ.IM0)GO TO 2
      DO 3 k=1,2
       JJ(4)=JJ(1)
       MM(4)=MM(1)
       IM(4)=IM(1)
       DO 4 i=1,3
        i1=i+1
        JJ(i)=JJ(i1)
        MM(i)=MM(i1)
4       IM(i)=IM(i1)
       IF(IM(3).EQ.IM0)GO TO 2
3      CONTINUE
!
! Change sign of m if m3<0
!
2     JS=JJ(1)+JJ(2)+JJ(3)
      UNM=1.0d0
      IF(IM(3).NE.MM(3))THEN
       UNM=1.0d0-2.0d0*MOD(JS,2)
       DO 5 i=1,3
5       MM(i)=-MM(i)
      ENDIF
!
! Initialisation of CG from m3=0
!
      M10=JJ(1)-MM(1)
      M1S=JJ(1)+MM(1)
      M20=JJ(2)-MM(1)
      M2S=JJ(2)+MM(1)
      DO 15 m=0,iMaxi
       DO 15 mp=0,iMaxi
15      CC(m,mp)=0.0d0
       M0=MAX0(-M1S,-M2S,0)
       MS=MIN0(M10,M20,MM(3))
       UN=1.0d0-2.0d0*MOD(iabs(MM(1)+M0),2)
       DO 10 m=M0,MS
        CC(MM(3),m)=UN*wig3j0(MM(1)+M,JJ(1),JJ(2),JJ(3))
10      UN=-UN
!
! Reccurence
!
      AJ1=real(JJ(1)*JJ(1)+JJ(1),dp)
      AJ2=real(JJ(2)*JJ(2)+JJ(2),dp)
      AJ3=real(JJ(3)*JJ(3)+JJ(3),dp)
      DO 20 m = (MM(3)-1),0,-1
       MM3=MM(3)-m
       C3=dSQRT(AJ3-MM3*MM3+MM3)
       MP0=MAX0(-M1S,-M2S-MM3,0)
       MPS=MIN0(M10,M20-MM3,m)
       DO 20 mp=MP0,MPS
        MM1=MM(1)+mp
        C1=dSQRT(AJ1-MM1*MM1-MM1)
        MM2=-MM1-MM3
        Cs2=dSQRT(AJ2-MM2*MM2-MM2)
20      CC(M,mp)=-(C1*CC(M+1,mp+1)+Cs2*CC(M+1,mp))/C3
      wigner3j=CC(0,0)*UNM
    return        
#endif
  END function wigner3j
!--------------------------------------------------------------------------------------------------------------------
  double precision function wig3j0(M,J1,J2,J)
!--------------------------------------------------------------------------------------------------------------------
!
!     calculation of the 3j coefficient:    
!                       m (j1 j2 j)  
!                    (-)  (m  -m 0)  
!
      use module_common_var  
      implicit none
      integer*8, intent(in) :: M, J1, J2, J
      !function var.
      integer*8, parameter         :: izero=0
      integer*8                    :: J0, JM, JS, JI
      integer*8                    :: MA
      integer*8                    :: i
      double precision             :: ACG, ACG0
      double precision             :: Ai, Ai1, Ai2
      double precision             :: AJ0, AJ02, AJM, AJM2
      double precision             :: CG0, CG1
      double precision             :: QM, UN, CombiJM
      !
      ! Init.
      !      
      wig3j0=0.0d0
      !
      ! 1st Zero rule:
      IF( (J1.lt.0) .OR. (J2.lt.0) .OR. (J.lt.0) ) RETURN
      !
      ! 2nd Zero rule:
      JM=J1+J2
      J0=abs(J1-J2)
      IF( (J.gt.JM) .OR. (J.lt.J0) ) RETURN
      !
      ! 3rd Zero rule:
      JS=MAX0(J1,J2)
      JI=MIN0(J1,J2)
      MA=iabs(M)
      IF( MA.gt.JI ) RETURN
      !
      ! w3j0 calcualtion:
      UN=1.0d0-2.0d0*MOD(JS,2)
      QM=real(M+M,dp)
      CG0=0.0d0
      wig3j0 = UN * &
               dSQRT(CombiJM(JI,MA)/CombiJM(JS,MA)*CombiJM(J0,izero)/(JS+JS+1.0d0))
      AJ0=real(J0,dp)
      AJM=real(JM+1,dp)
      AJ02=AJ0*AJ0
      AJM2=AJM*AJM
      ACG0=0.0d0
      DO 1 i=J0+1,J
        Ai=real(i,dp)
        Ai2=Ai*Ai
        ACG=dSQRT((AJM2-Ai2)*(Ai2-AJ02))
        CG1=(QM*(i+i-1.0d0)*wig3j0-ACG0*CG0)/ACG
        CG0=wig3j0
        wig3j0=CG1
1       ACG0=ACG
      RETURN
  END function wig3j0
!--------------------------------------------------------------------------------------------------------------------
  double precision function CombiJM(J,M)
!--------------------------------------------------------------------------------------------------------------------
!
!       calculation of: (2J)!/(J-M)!(J+M)!2**2J
!      ***********************************
      use module_common_var
      implicit none
      integer*8, intent(in) :: J, M
      ! In-function var.
      integer*8 :: i,k
      !
      CombiJM=1.0d0
      DO i=1,J
        CombiJM=CombiJM*(1.0d0-0.5d0/real(i,dp))
      ENDDO
      DO k=1,M
        CombiJM=(J+1.0d0-k)*CombiJM/(real(J+k,dp))
      ENDDO
      RETURN
  END function CombiJM
!--------------------------------------------------------------------------------------------------------------------                               
  double precision function wigner6j(A,B,uJ3,D,C6,F) 
!--------------------------------------------------------------------------------------------------------------------                                                                         
!wigner6j(uJ1,uJ2,uJ3,lJ1,lJ2,lJ3)                                                                 
!CALCUL OF 6J COEFFICIENTS FROM THE BOOKS OF MESSIAH, AND ROSE                    
!FOR THE CASES WHERE                                              
!           |  A  B  1  |       A+B+C+D                                
!           |  D  C  F  | = (-1)         W(ABCD;1F)                    
!           |_         _|                                              
!                                                                      
!WHERE THE RACAH COEFFICIENTS ARE TAKEN FROM THE BOOK OF ROSE (P.227)                    
!  
      use module_common_var
#if USE_EXTERNAL_WIGNER
      use ffastwigxj
#endif
      implicit none
      double precision :: A,B,uJ3,D,C6,F
#if USE_EXTERNAL_WIGNER
      integer (kind=4) :: J1, J2, J3, M1, M2, M3
#else
      integer (kind=8) :: E                                                 
      double precision :: term
#endif

#if USE_EXTERNAL_WIGNER
      J1 = 2*int(A,4) ; J2 = 2*int(B,4) ; J3 = 2*int(uJ3,4) 
      M1 = 2*int(D,4) ; M2 = 2*int(C6,4) ; M3 = 2*int(F,4) 
      wigner6j = ffw6jja(J1, J2, J3, M1, M2, M3)
      return
#else
!
!
      E = int(uJ3,8)
        
      
      IF( (abs(A-C6).GT.F) .OR. &
          (abs(B-D ).GT.F) .OR. &
          (A+C6     .LT.F) .OR. &      
          (B+D      .LT.F) .OR. &
          (E        .NE.1)) GOTO 1000  

      GOTO(1,2,3,1000) int(C6-D+2)                                             
!--------- CASE C=D-1 --------------                                     
   1  CONTINUE                                                          
      GOTO(10,11,12) int(A-B+2)                                               
!CASE A=B-1                                                             
10    CONTINUE                                                          
      TERM=((F+B+D+1.0_dp)*(F+B+D)*(B+D-F)*(B+D-(1+F)))                    
      TERM=TERM/(4.0_dp*(2*B+1)*B*(2*B-1)*D*(2*D-1)*(2*D+1))        
      wigner6j=((-1.0_dp)**(A+C6+F))*DSQRT(TERM)                                   
      RETURN                                                            
!CASE A=B                                                               
11    CONTINUE                                                          
      TERM=((F+B+D+1.0_dp)*(F+B-D+1)*(F+D-B)*(B+D-F))                      
      TERM=TERM/(4.0_dp*B*(2*B+1)*(B+1)*D*(2*D+1)*(2*D-1))           
      wigner6j=((-1.0_dp)**(A+C6-F-1))*DSQRT(TERM)                                 
      RETURN                                                            
!CASE A=B+1                                                             
12    CONTINUE                                                          
      TERM=((F+D-B)*(F+D-B-1.0_dp)*(F+B-D+2)*(F+B-D+1))                   
      TERM=TERM/(4.0_dp*(2*B+1)*(B+1)*(2*B+3)*(2*D-1)*D*(2*D+1))   
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!--------- CASE C=D   --------------                                     
   2  CONTINUE                                                          
      GOTO(20,21,22) int(A-B+2)                                               
!CASE A=B-1                                                             
20    CONTINUE                                                          
      TERM=((F+B+D+1)*(D+B-F)*(F+B-D)*(F+D-B+1))                      
      TERM=TERM/(4.0_dp*(2.0_dp*B+1)*B*(2*B-1)*D*(2*D+1)*(D+1))           
      wigner6j=((-1.0_dp)**(A+C6-F-1))*DSQRT(TERM)                                 
      RETURN                                                            
!CASE A=B                                                               
21    CONTINUE                                                          
      TERM=(B*(B+1.0_dp)+D*(D+1)-F*(F+1))                                 
      TERM=TERM/DSQRT(4.0_dp*B*(B+1)*(2*B+1)*D*(2*D+1)*(D+1))          
      wigner6j=((-1.0_dp)**(A+C6-F-1))*TERM                                       
      RETURN                                                            
!CASE A=B+1                                                             
22    CONTINUE                                                          
      TERM=((F+D+B+2.0_dp)*(F+B-D+1)*(B+D-F+1)*(F+D-B))                   
      TERM=TERM/(4.0_dp*(2*B+1)*(B+1)*(2*B+3)*D*(D+1)*(2*D+1))      
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!--------- CASE C=D+1 --------------                                     
   3  CONTINUE                                                          
      GOTO(30,31,32) int(A-B+2)                                               
!CASE A=B-1                                                             
30    CONTINUE                                                          
      TERM=((F+B-D)*(F+B-D-1.0_dp)*(F+D-B+2)*(F+D-B+1))                   
      TERM=TERM/(4.0_dp*(2*B+1)*B*(2*B-1)*(D+1)*(2*D+1)*(2*D+3))    
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!CASE A=B                                                               
31    CONTINUE                                                          
      TERM=((F+D+B+2.0_dp)*(B+D-F+1)*(F+D-B+1)*(F+B-D))                   
      TERM=TERM/(4.0_dp*B*(2*B+1)*(B+1)*(2*D+1)*(D+1)*(2*D+3))      
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!CASE A=B+1                                                             
32    CONTINUE                                                          
      TERM=((F+D+B+3.0_dp)*(F+B+D+2)*(B+D-F+2)*(B+D-F+1))                
      TERM=TERM/(4.0_dp*(2*B+3)*(B+1)*(2*B+1)*   &                      
                (2*D+3)*(D+1)*(2*D+1))                                       
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!                                                                       
! CASE DES 6J NULS                                                       
1000  wigner6j=0.0_dp                                                           
      RETURN        
#endif                                                    
  END function wigner6j
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE bubble_index(N, array, indxo, ad, econ)
!--------------------------------------------------------------------------------------------------------------------
! Bubble sorting algorithm                    * 
!*        (about n*n comparisons used).       *
!* ------------------------------------------ *
!* Reference: "A book on C By Al Kelley and   *
!* Ira Pohl, The Benjamin/Cummings Publishing *
!* Company, Inc, 1984" [BIBLI 09].            *
!*                                            *
!*                F90 version by J-P Moreau.  *
!*                    (www.jpmoreau.fr)       *
!* ------------------------------------------ *
  use module_common_var
  use module_error
	IMPLICIT NONE
	integer (kind=8), intent(IN)     :: N
	integer (kind=8), intent(INOUT)  :: indxo(N)
	double precision, intent(IN)     :: array(N)
  character(1)    , intent(IN)     :: ad ! ad = 'a' ; ad = 'd'
  type (dta_ERR)  , intent (inout) :: econ
	!Local var
	integer (kind=8)        	                 :: i, j, k
	double precision                 :: sorta(N), saux
  !-------------------------------------------------
	!Source code
	indxo=(/ (i,i=1,N)/)
	sorta=array
  if (ad .eq. 'a') then
  ! -> from lowest to highes| in ascending order|de menor a mayor
      do j=1, N-1
        do i=j+1,N
          if ( sorta(i) < sorta(j) ) then 
            k       = indxo(i)
            indxo(i)= indxo(j)
            indxo(j)= k
           
            saux     = sorta(i)
            sorta(i) = sorta(j)
            sorta(j) = saux 
          endif
        enddo
      enddo
  elseif (ad .eq. 'd') then
  ! -> in order of decreasing size|in descending order|de mayor a menor
      do j=1, N-1
        do i=j+1,N
        if ( sorta(i) > sorta(j) ) then 
            k       = indxo(i)
            indxo(i)= indxo(j)
            indxo(j)= k
           
            saux     = sorta(i)
            sorta(i) = sorta(j)
            sorta(j) = saux 
          endif
        enddo
      enddo
  else
    call errorBubble(econ)
  endif

  END SUBROUTINE bubble_index
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE ibubble_index(N, array, indxo, ad, econ)
!--------------------------------------------------------------------------------------------------------------------
! Bubble sorting algorithm                    * 
!*        (about n*n comparisons used).       *
!* ------------------------------------------ *
!* Reference: "A book on C By Al Kelley and   *
!* Ira Pohl, The Benjamin/Cummings Publishing *
!* Company, Inc, 1984" [BIBLI 09].            *
!*                                            *
!*                F90 version by J-P Moreau.  *
!*                    (www.jpmoreau.fr)       *
!* ------------------------------------------ *
  use module_common_var
  use module_error
	IMPLICIT NONE
	integer (kind=8), intent(IN)    :: N
	integer (kind=8), intent(INOUT) :: indxo(N)
	integer (kind=8), intent(IN)    :: array(N)
  character(1)    , intent(IN)    :: ad ! ad = 'a' ; ad = 'd'
  type (dta_ERR)  , intent(inout) :: econ
	!Local var
	integer (kind=8)        	       :: i, j, k
	integer (kind=8)                 :: sorta(N), saux
  !---------------------------------------
	!Source code
	indxo=(/ (i,i=1,N)/)
	sorta=array
  if (ad .eq. 'a') then
  ! in ascending order
      do j=1, N-1
        do i=j+1,N
          if ( sorta(i) < sorta(j) ) then 
            k       = indxo(i)
            indxo(i)= indxo(j)
            indxo(j)= k
           
            saux     = sorta(i)
            sorta(i) = sorta(j)
            sorta(j) = saux 
          endif
        enddo
      enddo
  elseif (ad .eq. 'd') then
  ! in descending order
      do j=1, N-1
        do i=j+1,N
        if ( sorta(i) > sorta(j) ) then 
            k       = indxo(i)
            indxo(i)= indxo(j)
            indxo(j)= k
           
            saux     = sorta(i)
            sorta(i) = sorta(j)
            sorta(j) = saux 
          endif
        enddo
      enddo
  else
    call errorBubble(econ)
  endif

  END SUBROUTINE ibubble_index
!--------------------------------------------------------------------------------------------------------------------
  integer*8 function imax(N, iarray)
!--------------------------------------------------------------------------------------------------------------------
! "imax": This will find a maximum within an array of integers.
!
	IMPLICIT NONE
	integer (kind=8), intent(IN)    :: N
	integer (kind=8), intent(IN)    :: iarray(N)
	!Local var
	integer (kind=8)        	       :: i, j, k
  !----------------------------------
	imax = iarray(1)
	do i= 2,N
		if (iarray(i) > imax) then
			imax = iarray(i)
		endif
	enddo
	
  END function imax
!--------------------------------------------------------------------------------------------------------------------
  integer*8 function imin(N, iarray)
!--------------------------------------------------------------------------------------------------------------------
! "imin": This will find a minimum within an array of integers.
!
	IMPLICIT NONE
	integer (kind=8), intent(IN)    :: N
	integer (kind=8), intent(IN)    :: iarray(N)
	!Local var
	integer (kind=8)        	       :: i, j, k
  !--------------------------------
	imin = iarray(1)
	do i= 2,N
		if (iarray(i) < imin) then
			imin = iarray(i)
		endif
	enddo
	
  END function imin
!--------------------------------------------------------------------------------------------------------------------
  double precision function maxf(N, array)
!--------------------------------------------------------------------------------------------------------------------
! "maxf": This will find a maximum within an array of reals. 
!
	IMPLICIT NONE
	integer (kind=8)         , intent(IN) :: N
	double precision, intent(IN) :: array(N)
	!Local var
	integer (kind=8)        	             :: i, j, k
  !---------------------------------------
	maxf = array(1)
	do i= 2,N
		if (array(i) > maxf) then
			maxf = array(i)
		endif
	enddo
	
  END function maxf
!--------------------------------------------------------------------------------------------------------------------
  double precision function minf(N, array)
!--------------------------------------------------------------------------------------------------------------------
! "minf": This will find a maximum within an array of reals. 
!
	IMPLICIT NONE
	integer (kind=8)         , intent(IN) :: N
  double precision, intent(IN) :: array(N)
	!Local var
	integer (kind=8)        	       :: i, j, k
  !---------------------------------------
	minf = array(1)
	do i= 2,N
		if (array(i) < minf) then
			minf = array(i)
		endif
	enddo
	
  END function minf
!--------------------------------------------------------------------------------------------------------------------
logical function isnan(a) 
!--------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
  double precision, intent(IN)  :: a 
!---------------------------------------
  if (a.ne.a) then 
    isnan = .true. 
  else 
    isnan = .false. 
  end if 
  return 
end function isnan
!--------------------------------------------------------------------------------------------------------------------
logical function isinf(a) 
!--------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
  double precision, intent(IN) :: a 
!---------------------------------------
  if ((a*0).ne.0) then 
    isinf = .true. 
  else 
    isinf = .false. 
  end if 
  return 
end function isinf
!--------------------------------------------------------------------------------------------------------------------
