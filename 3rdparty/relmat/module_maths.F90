!--------------------------------------------------------------------------------------------------------------------
module module_maths
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to 
!
    interface        

        logical function isJb(dta1,i1,i2)
            use module_common_var
            implicit none
            integer (kind=8), intent(in)             :: i1,i2
            type(dta_SDF),intent(in)        :: dta1  
        end function isJb
        
        integer*8 function compareC(C1, C12)
            implicit none
            character*2   , intent(in)      :: C1, C12
        end function compareC

        double precision function wigner3j( j1,j2,j3, m1,m2,m3 )
            use module_common_var
            implicit none
            double precision, intent(in)    :: j1, j2, j3
            double precision, intent(in)    :: m1, m2, m3
        end function wigner3j

        double precision function wig3j0(M,J1,J2,J)
            implicit double precision(a-h,o-z)
        end function wig3j0

        double precision function CombiJM(J,M)
            implicit double precision(a-h,o-z)
        END function CombiJM
        
        double precision function wigner6j(uJ1,uJ2,uJ3,lJ1,lJ2,lJ3)
            use module_common_var
            implicit none
            double precision, intent(in)    :: uJ1,uJ2,uJ3
            double precision, intent(in)    :: lJ1,lJ2,lJ3
        end function wigner6j

        double precision recursive function factr(n) RESULT(res)
            implicit none
            integer (kind=8), intent(in)             :: n
        end function factr
        
        double precision function fung(s, j1,j2,j3,Jd1,Jd2,Jd3)
            implicit none
            integer (kind=8), intent(in)             :: s
            integer (kind=8), intent(in)             :: j1 ,j2 ,j3
            integer (kind=8), intent(in)             :: Jd1,Jd2,Jd3
        end function fung
        
        double precision function triangle_coeff(a,b,c)
            implicit none
            integer (kind=8), intent(in)             :: a,b,c
        end function triangle_coeff
        
        subroutine bubble_index(N, array, indxo, ad)
            implicit none
            integer (kind=8)  , intent(in)           :: N
            integer (kind=8)  , intent(inout)        :: indxo(N)
            real*8   , intent(in)           :: array(N)
            character, intent(in)           :: ad
        end subroutine bubble_index
        
        subroutine ibubble_index(N, array, indxo, ad)
            implicit none
            integer (kind=8)  , intent(in)           :: N
            integer (kind=8)  , intent(inout)        :: indxo(N)
            integer (kind=8)  , intent(in)           :: array(N)
            character,intent(in)            :: ad
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
            real*8 , intent(IN)           :: array(N)
        end function maxf
        
        double precision function minf(N, array)
            implicit none
            integer (kind=8), intent(IN)           :: N
            real*8 , intent(IN)           :: array(N)
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
! Ci      : C part of the U/L State Local Quanta for CH4 (Input).
! ai      : "a" (or "n" above 3400cm-1) are counting integers.
!              for levels of the same J and C (Input).
! isJb    : it tells the result of the comparision:
!           1 -> greater
!           0 -> smaller
!          -1 -> same level
!
! Accessed Files:  None
! --------------
!
! Called Routines: 'compareC'   (Partition Function of CH4)
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
    integer*8, intent(in)     :: i1,i2
    type(dta_SDF),intent(in):: dta1
    !integer*8                  :: compareC
    integer*8                 :: J1, J2
    !integer*8                  :: alph1, alph2
    !character*2              :: C1, C12
!-----------------------------------------
! Q001:
    J1 = dta1%J(i1,1)!; C1 = dta1%Sr(i1,1); alph1= dta1%alph(i1,1)
! Q002:
    J2 = dta1%J(i2,1)!; C12= dta1%Sr(i2,1); alph2= dta1%alph(i2,1)
    if (J1 .ge. J2) then
      isJb = .true.
    else
      isJb = .false.
    endif

!    if (J1 .gt. J2) then
!      isJb = 1
!    else if (J1 .eq. J2) then
!      if (compareC(C1, C12) .eq. 1) then
!        isJb = 1
!      else if (compareC(C1, C12) .eq. -1) then
!        if (alph1 .gt. alph2) then
!          isJb = 1
!        else if (alph1 .lt. alph2) then
!          isJb = 0
!        else
!          !same initial level
!          isJb = -1
!        endif
!      else
!        isJb = 0
!      endif
!    else
!      isJb = 0
!    endif
    Return
  END function isJb
!--------------------------------------------------------------------------------------------------------------------
  integer*8 function compareC(C1, C12)
!--------------------------------------------------------------------------------------------------------------------
!"compareC":  
!
! Detailed description:
! ---------------------
! it compares the C part of the U/L State Local Quanta.
! the rotational quanta is given by J, C and ?a".
!
! 1) C can have the following values:
!               A1 > A2 > E > F1 > F2         (for CH4)
! 2) "a" (or "n" above 3400cm-1) are counting integer (kind=8)levels of the same J and C; 
! NOTE: the values are incremented in order of increasing energy.
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! Ci      : C part of the U/L State Local Quanta for CH4 (Input).
!     compareC: it tells the result of the comparision:
!     1 -> greater
!     0 -> smaller
!    -1 -> same level
!
! Accessed Files:  None
! --------------
!
! Called Routines: 'compareC'   (Partition Function of CH4)
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
    implicit none
    character*2, intent(in)     :: C1, C12
    integer (kind=8)            :: i, N, auxC1, auxC2
    character*2                 :: a(11)
!----------
    DATA a/'A ','A1','A2','A+','A-','E ','E1', 'E2','F ','F1','F2'/
!----------
!
! Check that Arrays for Results are Large Enough, Initialize
!  
    N = len(a) ! "a" length
    do i=1,N
      if ( adjustL(trim(C1)) .eq. adjustL(trim(a(i))) ) then
        auxC1 = i
      endif
      if ( adjustL(trim(C12)) .eq. adjustL(trim(a(i))) ) then
        auxC2 = i
      endif
    enddo
    
    if (auxC1 .gt. auxC2) then
      compareC = 1
      !print*, 'C1', '>', 'C2'
    else if (auxC1 .lt. auxC2) then
      compareC = 0
      !print*, 'C1', '<', 'C2'
    else !(auxC1 == auxC2)
      compareC = -1
      !print*, 'C1', '=', 'C2'
    endif
!
!
    Return
  END function compareC
!--------------------------------------------------------------------------------------------------------------------
  double precision function wigner3j(dJ1,dJ2,dJ3,dM1,dM2,dM3)
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

      
      implicit double precision(a-h,o-z)
      implicit integer*8 (i-n)
      PARAMETER (IMAX=100_Int_8)
      DIMENSION JJ(4),MM(4),IM(4),CC(0:IMAX,0:IMAX)

      J1 = dJ1 ; J2 = dJ2 ; J3 = dJ3 
      M1 = dM1 ; M2 = dM2 ; M3 = dM3 

      wigner3j=0.0_dp
      JJ(1)=J1 ; JJ(2)=J2 ; JJ(3)=J3
      MM(1)=M1 ; MM(2)=M2 ; MM(3)=M3
      MS=0
      DO 1 I=1,3
        MS=MS+MM(I)
        IM(I)=abs(MM(I))
        IF( .NOT.(JJ(I) .GE. 0 .AND. IM(I) .LE. JJ(I)) )RETURN
1     CONTINUE     
!
! Other selection rules
!
      IF( MS .NE. 0 )RETURN
      J0=abs(J1-J2)
      JM=J1+J2
      IF( J3.GT.JM .OR. J3.LT.J0 )RETURN
!
! Parameter test 
!
      IM0=MIN0(IM(1),IM(2),IM(3))
      IF(IM0.GT.IMAX)THEN
        write(8,100)IM0,IMAX
100    FORMAT(' min(abs(m))=',I2,' > IMAX=',I2)
        RETURN
      ENDIF
!
! Circular Permutations amenant le plus petit abs(m) a droite
!
      IF(IM(3).EQ.IM0)GO TO 2
      DO 3 K=1,2
       JJ(4)=JJ(1)
       MM(4)=MM(1)
       IM(4)=IM(1)
       DO 4 I=1,3
        I1=I+1
        JJ(I)=JJ(I1)
        MM(I)=MM(I1)
4       IM(I)=IM(I1)
       IF(IM(3).EQ.IM0)GO TO 2
3     CONTINUE
!
! Changement de signe des m si m3<0
!
2     JS=JJ(1)+JJ(2)+JJ(3)
      UNM=1.
      IF(IM(3).NE.MM(3))THEN
       UNM=1-2*MOD(JS,2)
       DO 5 I=1,3
5       MM(I)=-MM(I)
      ENDIF
!
! Initialisation des CG pour m3=0
!
      M10=JJ(1)-MM(1)
      M1S=JJ(1)+MM(1)
      M20=JJ(2)-MM(1)
      M2S=JJ(2)+MM(1)
      DO 15 M=0,IMAX
       DO 15 MP=0,IMAX
15      CC(M,MP)=0.
       M0=MAX0(-M1S,-M2S,0)
       MS=MIN0(M10,M20,MM(3))
       UN=1-2*MOD(abs(MM(1)+M0),2)
       DO 10 M=M0,MS
        CC(MM(3),M)=UN*wig3j0(MM(1)+M,JJ(1),JJ(2),JJ(3))
10      UN=-UN
!
! Reccurence
!
      AJ1=JJ(1)*JJ(1)+JJ(1)
      AJ2=JJ(2)*JJ(2)+JJ(2)
      AJ3=JJ(3)*JJ(3)+JJ(3)
      DO 20 M=MM(3)-1,0,-1
       MM3=MM(3)-M
       C3=dSQRT(AJ3-MM3*MM3+MM3)
       MP0=MAX0(-M1S,-M2S-MM3,0)
       MPS=MIN0(M10,M20-MM3,M)
       DO 20 MP=MP0,MPS
        MM1=MM(1)+MP
        C1=dSQRT(AJ1-MM1*MM1-MM1)
        MM2=-MM1-MM3
        Cs2=dSQRT(AJ2-MM2*MM2-MM2)
20      CC(M,MP)=-(C1*CC(M+1,MP+1)+Cs2*CC(M+1,MP))/C3
      wigner3j=CC(0,0)*UNM
    return        
  END function wigner3j
!--------------------------------------------------------------------------------------------------------------------
  double precision function wig3j0(M,J1,J2,J)
!--------------------------------------------------------------------------------------------------------------------
!
!     calculation of the 3j coefficient:    
!                       m (j1 j2 j)  
!                    (-)  (m  -m 0)  
!
      implicit double precision(a-h,o-z)
      implicit integer*8 (i-n)
      integer, parameter  :: Int_8  = selected_int_kind(8)
      
      izero=0
      wig3j0=0.
      IF(J1.LT.0 .OR. J2.LT.0 .OR. J.LT.0)RETURN
      JM=J1+J2
      J0=abs(J1-J2)
      IF(J.GT.JM .OR. J.LT.J0)RETURN
      JS=MAX0(J1,J2)
      JI=MIN0(J1,J2)
      MA=abs(M)
      IF(MA.GT.JI)RETURN
      UN=1-2*MOD(JS,2)
      QM=M+M
      CG0=0.
      wig3j0=UN*dSQRT(CombiJM(JI,MA)/CombiJM(JS,MA)* &
           CombiJM(J0,izero)/(JS+JS+1))
      AJ0=J0
      AJM=JM+1
      AJ02=AJ0*AJ0
      AJM2=AJM*AJM
      ACG0=0.
      DO 1 I=J0+1,J
        AI=I
        AI2=AI*AI
        ACG=dSQRT((AJM2-AI2)*(AI2-AJ02))
        CG1=(QM*(I+I-1)*wig3j0-ACG0*CG0)/ACG
        CG0=wig3j0
        wig3j0=CG1
1     ACG0=ACG
  END function wig3j0
!--------------------------------------------------------------------------------------------------------------------
  double precision function CombiJM(J,M)
!--------------------------------------------------------------------------------------------------------------------
!
!       calculation of: (2J)!/(J-M)!(J+M)!2**2J
!      ***********************************
      implicit double precision(a-h,o-z)
      implicit integer*8 (i-n)
      
      CombiJM=1.
      DO 1 I=1,J
1     CombiJM=CombiJM*(1.-.5/I)
      DO 2 K=1,M
2     CombiJM=(J+1-K)*CombiJM/(J+K)
  END function CombiJM
!--------------------------------------------------------------------------------------------------------------------                               
  double precision function wigner6j(uJ1,uJ2,uJ3,lJ1,lJ2,lJ3)    
!--------------------------------------------------------------------------------------------------------------------                                                                         
!                                                                      
!CALCUL OF 6J COEFFICIENTS FROM THE BOOKS OF MESSIAH, AND ROSE                    
!FOR THE CASES WHERE                                              
!           |  A  B  1  |       A+B+C+D                                
!           |  D  C  F  | = (-1)         W(ABCD;1F)                    
!           |_         _|                                              
!                                                                      
!WHERE THE RACAH COEFFICIENTS ARE TAKEN FROM THE BOOK OF ROSE (P.227)                    
!  
      use module_common_var
      implicit none
      double precision :: uJ1,uJ2,uJ3,lJ1,lJ2,lJ3
      integer (kind=8) :: A,B,C6,D,E,F                                                 
      double precision :: term
!
!
      A = int(uJ1) ; B = int(uJ2) ; E = int(uJ3)
      D = int(lJ1) ; C6 = int(lJ2) ; F = int(lJ3)  
      IF((abs(A-C6).GT.F).OR. &
        (abs(B-D).GT.F) .OR. &
        (A+C6.LT.F)       .OR. &      
        (B+D.LT.F)) GOTO 1000  

      GOTO(1,2,3,1000) C6-D+2                                             
!--------- CASE C=D-1 --------------                                     
   1  CONTINUE                                                          
      GOTO(10,11,12)A-B+2                                               
!CASE A=B-1                                                             
10    CONTINUE                                                          
      TERM=((F+B+D+1.0_dp)*(F+B+D)*(B+D-F)*(B+D-(1.+F)))                    
      TERM=TERM/(4.0_dp*(2.*B+1.)*B*(2.*B-1.)*D*(2.*D-1.)*(2.*D+1.))        
      wigner6j=((-1.0_dp)**(A+C6+F))*DSQRT(TERM)                                   
      RETURN                                                            
!CASE A=B                                                               
11    CONTINUE                                                          
      TERM=((F+B+D+1.0_dp)*(F+B-D+1.)*(F+D-B)*(B+D-F))                      
      TERM=TERM/(4.0_dp*B*(2.*B+1.)*(B+1.)*D*(2.*D+1.)*(2.*D-1.))           
      wigner6j=((-1.0_dp)**(A+C6-F-1))*DSQRT(TERM)                                 
      RETURN                                                            
!CASE A=B+1                                                             
12    CONTINUE                                                          
      TERM=((F+D-B)*(F+D-B-1.0_dp)*(F+B-D+2.)*(F+B-D+1.))                   
      TERM=TERM/(4.0_dp*(2.*B+1.)*(B+1.)*(2.*B+3.)*(2.*D-1.)*D*(2.*D+1.))   
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!--------- CASE C=D   --------------                                     
   2  CONTINUE                                                          
      GOTO(20,21,22)A-B+2                                               
!CASE A=B-1                                                             
20    CONTINUE                                                          
      TERM=((F+B+D+1.)*(D+B-F)*(F+B-D)*(F+D-B+1.))                      
      TERM=TERM/(4.0_dp*(2.0_dp*B+1.)*B*(2.*B-1.)*D*(2.*D+1.)*(D+1.))           
      wigner6j=((-1.0_dp)**(A+C6-F-1))*DSQRT(TERM)                                 
      RETURN                                                            
!CASE A=B                                                               
21    CONTINUE                                                          
      TERM=(B*(B+1.0_dp)+D*(D+1.)-F*(F+1.))                                 
      TERM=TERM/DSQRT(4.0_dp*B*(B+1.)*(2.*B+1.)*D*(2.*D+1.)*(D+1.))          
      wigner6j=((-1.0_dp)**(A+C6-F-1))*TERM                                       
      RETURN                                                            
!CASE A=B+1                                                             
22    CONTINUE                                                          
      TERM=((F+D+B+2.0_dp)*(F+B-D+1.)*(B+D-F+1.)*(F+D-B))                   
      TERM=TERM/(4.0_dp*(2.*B+1.)*(B+1.)*(2.*B+3.)*D*(D+1.)*(2.*D+1.))      
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!--------- CASE C=D+1 --------------                                     
   3  CONTINUE                                                          
      GOTO(30,31,32)A-B+2                                               
!CASE A=B-1                                                             
30    CONTINUE                                                          
      TERM=((F+B-D)*(F+B-D-1.0_dp)*(F+D-B+2.)*(F+D-B+1.))                   
      TERM=TERM/(4.0_dp*(2.*B+1.)*B*(2.*B-1.)*(D+1.)*(2.*D+1.)*(2*D+3.))    
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!CASE A=B                                                               
31    CONTINUE                                                          
      TERM=((F+D+B+2.0_dp)*(B+D-F+1.)*(F+D-B+1.)*(F+B-D))                   
      TERM=TERM/(4.0_dp*B*(2.*B+1.)*(B+1.)*(2.*D+1.)*(D+1.)*(2.*D+3.))      
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!CASE A=B+1                                                             
32    CONTINUE                                                          
      TERM=((F+D+B+3.0_dp)*(F+B+D+2.)*(B+D-F+2.)*(B+D-F+1.))                
      TERM=TERM/(4.0_dp*(2.*B+3.)*(B+1.)*(2.*B+1.)*   &                      
                (2.*D+3.)*(D+1.)*(2.*D+1.))                                       
      wigner6j=((-1.0_dp)**(A+C6-F))*DSQRT(TERM)                                   
      RETURN                                                            
!                                                                       
! CASE DES 6J NULS                                                       
1000  wigner6j=0.0_dp                                                           
      RETURN                                                            
  END function wigner6j
!--------------------------------------------------------------------------------------------------------------------
  double precision recursive function factr(n) RESULT(res)
!--------------------------------------------------------------------------------------------------------------------
! "factr": This is the Factorial function
! program in a recursive way. Double precission version.
!
      IMPLICIT NONE
      integer (kind=8) :: n
      !real*8  :: res
      IF (n .EQ. 0) THEN
        res = 1.0d0
      ELSE
        res = n * factr(n - 1)
      END IF
  END function factr
!--------------------------------------------------------------------------------------------------------------------
  double precision function fung(s, j1,j2,j3,Jd1,Jd2,Jd3)
!--------------------------------------------------------------------------------------------------------------------
! "fung": it calculates the denominator in Racah Formula.
!
      IMPLICIT NONE
      integer (kind=8),intent(in) :: s
      integer (kind=8),intent(in) :: j1 ,j2 ,j3
      integer (kind=8),intent(in) :: Jd1,Jd2,Jd3
      double precision   :: factr
! --------------------------------------
      fung = factr(s-j1-j2-j3)                              &
           * factr(s-j1-Jd2-J3)    * factr(s-Jd1-j2-Jd3)    &
           * factr(s-Jd1-Jd2-j3)   * factr(j1+j2+Jd1+Jd2-s) &
           * factr(j2+j3+Jd2+Jd3-s)* factr(j3+j1+Jd3+Jd1-s)
      RETURN
  END function fung
!--------------------------------------------------------------------------------------------------------------------
  double precision function triangle_coeff(a,b,h)
!--------------------------------------------------------------------------------------------------------------------
! "triangle_coeff": Calculates triangle coefficients for angular momenta.
! NOTE:
! -----
! This version returns 0 if the triangle inequalities are violated.  (RAH)
!
      use module_common_var
      IMPLICIT NONE
      integer (kind=8), intent(IN)  :: a,b,h
      integer (kind=8)              :: xa
      double precision             :: if1, if2, if3, if4
      double precision             :: factr
! -----------------------------------------------
      if ( (a .le. 0) .or. (b .le. 0) .or. (h .le. 0) ) then
       !print*, "a, b or h = 0"
	     triangle_coeff=0.0
      else
        do xa = abs(a-b),(a+b),1
          !print*, xa
	       if (h .eq. xa) then
            !print*, "h/xa:", h, xa
            if1 = factr(a+b-h)
            !print*, a+b-h, "=a+b-h; Factorial=", if1
            if2 = factr(a-b+h)
            !print*, a-b+h,"=a-b+h; Factorial=", if2
            if3 = factr(b+h-a)
            !print*, b+h-a,"=-a+b+h; Factorial=", if3
            if4 = factr(a+b+h+1)
            !print*, a+b+h+1,"=a+b+h+1; Factorial=", if4
	          triangle_coeff = dsqrt(if1*if2*if3/if4)
            !print*, "triangleCoeff=", triangle_coeff
	        endif
        enddo
      endif
      RETURN
  END function triangle_coeff
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE bubble_index(N, array, indxo, ad)
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
	IMPLICIT NONE
	integer (kind=8),  intent(IN)             :: N
	integer (kind=8),  intent(INOUT)          :: indxo(N)
	double precision, intent(IN)     :: array(N)
  character,intent(IN)             :: ad ! ad = 'a' ; ad = 'd'
	!Local var
	integer (kind=8)        	                 :: i, j, k
	double precision                 :: sorta(N), saux
  !-------------------------------------------------
	!Source code
	indxo=(/ (i,i=1,N)/)
	!print*, array(1), array(N)
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
    print*, 'Subroutine bubble_index: not supported option, use a/d instead'
  endif

!  DO i=1,N
!    print*, i, indxo(i), array(indxo(i)), sorta(i) 
!  ENDDO


!  print*, "First"
!  print*, "non-ordered:",1, array(1) 
!  print*, "ordered    :",indxo(1),sorta(1) 
!  print*,'----------------------'
!  print*, "Last"
!  print*, "non-ordered:",N, array(N) 
!  print*, "ordered    :",indxo(N),sorta(N) 
!  print*,'----------------------'
!  stop 

	
  END SUBROUTINE bubble_index
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE ibubble_index(N, array, indxo, ad)
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
	IMPLICIT NONE
	integer (kind=8), intent(IN)    :: N
	integer (kind=8), intent(INOUT) :: indxo(N)
	integer (kind=8), intent(IN)    :: array(N)
  character,intent(IN)   :: ad ! ad = 'a' ; ad = 'd'
	!Local var
	integer (kind=8)        	       :: i, j, k
	integer (kind=8)                :: sorta(N), saux
  !---------------------------------------
	!Source code
	indxo=(/ (i,i=1,N)/)
	!print*, array(1), array(N)
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
    print*, 'Subroutine bubble_index: not supported option, use a/d instead'
  endif
  
!  DO i=1,N
!    print*, i, array(i), indxo(i), array(indxo(i)), sorta(i) 
!  ENDDO

!  print*, "First"
!  print*, "non-ordered:",1, array(1) 
!  print*, "ordered    :",indxo(1),sorta(1) 
!  print*,'----------------------'
!  print*, "Last"
!  print*, "non-ordered:",N, array(N) 
!  print*, "ordered    :",indxo(N),sorta(N)
!  print*,'----------------------' 
!  stop
	
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
