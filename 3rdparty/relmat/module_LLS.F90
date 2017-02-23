!--------------------------------------------------------------------------------------------------------------------
module module_LLS
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to 
!
    interface        

        subroutine LLS_Matrix(j,k,h,molP,PerM,M, econ)
            use module_common_var
            use module_maths
            use module_phsub
            implicit none
            integer       , intent(in)   :: j,k
            type (dta_SDF), intent(in)   :: h
            type (dta_MOL), intent(in)   :: molP, PerM
            type (dta_ERR), intent(inout):: econ
            double precision, intent(out):: M(4)
        end subroutine LLS_Matrix
        
        subroutine calc_QParam(nLines, dta1, molP, PerM, econ)
            use module_common_var
            use module_error
            use module_maths
            Implicit None
            Integer,        intent(in)      :: nLines
            type (dta_SDF), intent(in)      :: dta1
            type (dta_MOL), intent(inout)   :: molP
            type (dta_MOL), intent(in)      :: PerM
            type (dta_ERR), intent(inout)   :: econ
        end subroutine calc_QParam

    end interface

END module module_LLS
!--------------------------------------------------------------------------------------------------------------------
! MATHEMATIC FUNCTIONS AND SUBROUTINES ------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE LLS_Matrix(j,k,h,molP,PerM,M,econ)
!--------------------------------------------------------------------------------------------------------------------
! "LLS_Matrix": Subroutine that fills up the matrix to be included in the LLS                   
!------------------------------------------ 
!vi := initial state (downwards transition j->k)
!CASE:  J(j) > J(k)
!M(j,k,:) = [ K1*Sum(K2(L)), 
!              -K1*Sum( K2(L) * (L*(L+1)) ), ...
!              -c2*K1*Sum( K2(L) * B0(L*(L+1)) ,...
!               K1*Sum(K2(L))];      
!------------------------------------------ 
    use module_common_var
    use module_maths
    use module_phsub
    IMPLICIT none
    ! a1, a2, a3, dc were declared in module_molecSP.SystemQparam:
    integer*8     , intent(in):: j,k
    type (dta_SDF), intent(in):: h
    type (dta_MOL), intent(in):: molP, PerM
    type (dta_ERR), intent(inout):: econ
    double precision, intent(out):: M(4)
    !internal variables:
    integer*8             :: L, incr, step
    integer*8             :: li,lf
    integer*8             :: i, iniL,endL
    integer*8             :: Ni, Ni_p, Nf, Nf_p
    double precision      :: Ji, Ji_p, Jf, Jf_p !, jmax
    double precision      :: Si, Sf
    double precision      :: dE,Jaux
    double precision      :: K1, K2, AF1, AF2
    double precision      :: T, Ptot, RT
    !-----------------------------------------
      T    = molP % Temp
      Ptot = molP % Ptot
      RT   = T0/T
      AF1  = 1.0_dp
      AF2  = 1.0_dp
    !-----------------------------------------
! molP%a1 - (molP%a2)*(E_l) - molP%a3*c2*molP%B0*E_l/T + 1.0_dp
      step = 2 ! -> Because the quadrupole-quadrupole energy potential is dominant
      !             in linear molecules.
      ! j -> Ji(j) > Ji(k)=Ji'
      Ji   = h%J(j,1); Ni = h%N(j,1); Si = h%espin(j,1)
      Jf   = h%J(j,2); Nf = h%N(j,2); Sf = h%espin(j,2)
      ! k-> 
      Ji_p = h%J(k,1); Ni_p = h%N(k,1)
      Jf_p = h%J(k,2); Nf_p = h%N(k,2)
!
        ! li and lf are the vibrational angular momentum of the vibrational
        ! available for Linear molecules in HITRAN; 
        ! (due to its symmetry though).
!
       li = h%lv2(1)
       lf = h%lv2(2) 
! 
      iniL=int(max(abs(Ji-Ji_p),abs(Jf-Jf_p)))
      If( mod(iniL,step) .ne. 0 )iniL=iniL+1
      endL=int(min((Ji+Ji_p),(Jf+Jf_p)))
!
      M(1)=0.0_dp;M(2)=0.0_dp;M(3)=0.0_dp;M(4)=0.0_dp
      if (molP%M .eq. 7) then
          if (molP%AF_ON) then ! Adiabatic factor 1:
            AF1 = AFmol_X(molP, PerM, real(Ni,dp), step)
          endif
          K1 = Kpart1_O2(Ji, Ji_p, Jf, Jf_p, &
                         Ni, Ni_p, Nf, Nf_p, AF1)
      else
          if (molP%AF_ON) then ! Adiabatic factor 1:
            AF1 = AFmol_X(molP, PerM, Ji, step)
          endif
          K1 = Kpart1(Ji_p, Jf, Jf_p, li, lf, AF1 )
      endif
      !print*, K1
      do L = iniL,endL,step
        ! Since the molecules to be considered are linear -> symmetric 
        ! and the quadrupole - quadrupole interaction for the systems mol-X
        ! are meant to be dominant the interaction potential is likely 
        ! dominated by EVEN rank contributions. 
        !
        dE = abs( real(L*(L+1),dp) )
        if (molP%AF_ON) then
          AF2 = AFmol_X(molP, PerM, real(L,dp), step)
        endif
        if (molP%M .eq. 7) then
          K2 = Kpart2_O2(L, Ji, Ji_p, Jf, Jf_p, &
                         Ni, Ni_p, Nf, Nf_p, Si, Sf, AF2, econ)
        else
          K2 = Kpart2( L, Ji, Jf, Ji_p, Jf_p, li, lf, AF2, econ)
        endif
        !
        M(1) = K2 + M(1)
        M(2) = K2*log(dE) + M(2)
        M(3) = K2*molP%B0*dE + M(3)
      enddo
      Jaux = Ji
      if (Jaux .eq. 0) Jaux = 0.5_dp
      ! Approximations:
      ! 1) Linear
      M(1) = K1 * M(1) 
      ! 2) correction 1
      !M(1) = K1 * M(1)/( (dabs(h%Sig(j)-h%Sig(k))/ molP%B0) )**(Ji_p/Jaux) ! correction 1 -> inverted results
      ! 2) correction 2
      !M(1) = K1 * M(1) / ( dabs(h%Sig(j)-h%Sig(k)) )**(Ji_p/Jaux) ! correction 2
      ! 2) correction 3
      !M(1) = K1 * M(1) / (T*( dabs(h%Sig(j)-h%Sig(k)) )**(Ji_p/Jaux)) ! correction 3 c2 (cm·K)
      ! 2) correction 4
      !M(1) = K1 * M(1) *molP%B0*molP%B0*c2/ (T*( dabs(h%Sig(j)-h%Sig(k)) )**(Ji_p/Jaux)) ! correction 4 K·cm-1
      M(2) = K1 * M(2)
      M(3) = K1 * (c2/T) * M(3)
      M(4) = M(1)


  END SUBROUTINE LLS_Matrix
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_QParam(nLines, dta1, molP, PerM, econ)
!--------------------------------------------------------------------------------------------------------------------
! "calc_QParam": gives the values of a1, a2, a3 
! after solving a Linear-Least-Square system. For that Lapack is used.
!
!------------------------------------------ 
    use module_common_var
    use module_error
    use module_maths
    Implicit None
    Integer*8,      intent(in)      :: nLines
    type (dta_SDF), intent(in)      :: dta1
    type (dta_MOL), intent(inout)   :: molP
    type (dta_MOL), intent(in)      :: PerM
    type (dta_ERR), intent(inout)   :: econ
    double Precision                :: Mlls(nLines,4)
    double Precision                :: Aux_4M(4)
    integer*8                       :: i, j, k
    integer*8                       :: jBIG, jSMALL
    integer*8                       :: indexI(nLines)
    Double Precision                :: r_kj, rD0_kj
    !Auxiliar Constants
    double Precision                :: Kaux, HWT, faH
    double precision                :: T, Ptot, RT
    !Auxiliar for LAPACK
    !
    ! Parameters:
      integer*8, Parameter   :: N = 3
      integer*8, Parameter   :: NRHS = 1
      integer*8              :: M
      integer*8              :: LDA, LDB
      !PARAMETER        ( LDA = M, LDB = M )
      integer*8, Parameter   :: LWMAX = 100
    ! Local Scalars:
      integer*8              :: INFO, LWORK
    ! Local Arrays:
      Double Precision       :: WORK( LWMAX )
      Double Precision,ALLOCATABLE,DIMENSION(:,:) :: A( :, : ), B( :, : )
!-----------------------------------------
      T    = molP % Temp
      Ptot = molP % Ptot
      RT   = T0/T
!-----------------------------------------
      M = nLines
      LDA = M
      LDB = M
      allocate ( A( LDA, N ), B( LDB, NRHS ) )
!
! LAPACK is used (installation command for mac):
! sudo port install lapack
! 
! * Compilation for "Free PGI compiler" for MAC:
! pgf90 myprog.f90 -llapack -lblas
!
! * Compilation "gfortran"
! gfortran myprog.f90 -llapack
!
!---------
! FIRST: create a zero matrix for 
    do j=1, nLines
      do k=1,4
        Mlls(j,k) = 0.0_dp
      enddo  
    enddo  
!
!
! Generate the Matrix for LLS:
    do j=1, nLines
      do k=1, nLines
        ! 
        if (j .eq. k) then
          faH = 1.0_dp
          if(T.ne.T0)faH = (RT**dta1%BHW(j)) 
          B(j,NRHS) = 2*molP%Ptot*dta1%HWT0(j)*faH
        else              
          if (isJb(dta1,j,k)) then
          ! CASE:  J(j) > J(k) (downwards transition j->k)
          ! or
          ! CASE: J(j) = J(k)
            jBIG   = j
            jSMALL = k
            r_kj = 1.0_dp 
          else
          ! CASE: J(j) < J(k)
          ! so downwards transition is (k->j)
          ! pj·<<k|W|j>> = pk·<<j|W|k>>; pk = dta1%PopuT(k); pj = dta1%PopuT(j)
            jBIG   = k
            jSMALL = j
            r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
          endif
          call LLS_Matrix(jBIG,jSMALL,dta1,molP,PerM,Aux_4M, econ)
          !
          rD0_kj = dta1%D0(k)/dta1%D0(j) 
          do i =1,4
            Mlls(j,i) = Mlls(j,i) + rD0_kj*r_kj*Aux_4M(i)
          enddo  
        endif    
      enddo
      !
    enddo
!
! **********************************************************************************
! LAPACK routine:
! ---------------
! DGELS solves overdetermined or undetermined real linear systems
! involving an M-by-N matrix A, or its transpose, using RQ or LQ factorization of A.
! Note_ it is asumed that A has full rank.
! if TRANS = 'N' and m>=n: find the least squares solution of an overdetermined system, 
! i.e., solve the least squares problem minimize || B - A*X ||
!  A * X = B,
!
!
! Definition of A, B:
    !
    A = -Mlls(1:nLines,1:3)
    !
    do i = 1,nLines
      B(i,NRHS) = B(i,NRHS) + Mlls(i,4)
    enddo    
!
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
!      CHARACTER          TRANS
!      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
!      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELS solves overdetermined or underdetermined real linear systems
!  involving an M-by-N matrix A, or its transpose, using a QR or LQ
!  factorization of A.  It is  that A has full rank.
!
!  The following options are provided:
!
!  1. If TRANS = 'N' and m >= n:  find the least squares solution of
!     an overdetermined system, i.e., solve the least squares problem
!                  minimize || B - A*X ||.
!
!  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!  an underdetermined system A * X = B.
!
!  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
!  an undetermined system A**T * X = B.
!
!  4. If TRANS = 'T' and m < n:  find the least squares solution of
!  an overdetermined system, i.e., solve the least squares problem
!               minimize || B - A**T * X ||.
!
!  Several right hand side vectors b and solution vectors x can be
!  handled in a single call; they are stored as the columns of the
!  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!  matrix X.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!       = 'N': the linear system involves A;
!       = 'T': the linear system involves A**T.
!
!  M       (input) INTEGER
!       The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!       The number of columns of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!       The number of right hand sides, i.e., the number of
!       columns of the matrices B and X. NRHS >=0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!       On entry, the M-by-N matrix A.
!       On exit,
!         if M >= N, A is overwritten by details of its QR
!                    factorization as returned by DGEQRF;
!         if M <  N, A is overwritten by details of its LQ
!                    factorization as returned by DGELQF.
!
!  LDA     (input) INTEGER
!       The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!       On entry, the matrix B of right hand side vectors, stored
!       columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!       if TRANS = 'T'.
!       On exit, if INFO = 0, B is overwritten by the solution
!       vectors, stored columnwise:
!       if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!       squares solution vectors; the residual sum of squares for the
!       solution in each column is given by the sum of squares of
!       elements N+1 to M in that column;
!       if TRANS = 'N' and m < n, rows 1 to N of B contain the
!       minimum norm solution vectors;
!       if TRANS = 'T' and m >= n, rows 1 to M of B contain the
!       minimum norm solution vectors;
!       if TRANS = 'T' and m < n, rows 1 to M of B contain the
!       least squares solution vectors; the residual sum of squares
!       for the solution in each column is given by the sum of
!       squares of elements M+1 to N in that column.
!
!  LDB     (input) INTEGER
!       The leading dimension of the array B. LDB >= MAX(1,M,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!       On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!       The dimension of the array WORK.
!       LWORK >= max( 1, MN + max( MN, NRHS ) ).
!       For optimal performance,
!       LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
!       where MN = min(M,N) and NB is the optimum block size.
!
!       If LWORK = -1, then a workspace query is assumed; the routine
!       only calculates the optimal size of the WORK array, returns
!       this value as the first entry of the WORK array, and no error
!       message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!       = 0:  successful exit
!       < 0:  if INFO = -i, the i-th argument had an illegal value
!       > 0:  if INFO =  i, the i-th diagonal element of the
!             triangular factor of A is zero, so that A does not have
!             full rank; the least squares solution could not be
!             computed.
!
!  =====================================================================
!
!  Command sentence:
!  CALL dgels( TRANS, M , N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
!  CALL dgelsd(TRANS, M , N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO ) 
!
! To use LAPACK uncomment this section -------->
!
!
!    write(*,2017),nLines
!2017 format("calling Lapack #",i3,"^2 square-matrix to invert")
!
!
!
!   Query the optimal workspace.
!
      LWORK = -1
      CALL DGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK,&
                 LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!   Solve the equations A*X = B.
!
      CALL DGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK,&
                 LWORK, INFO )
! <---------------------------------------------
! NOTE: for further information about this subroutine 
! http://www.netlib.no/netlib/lapack/double/dgels.f
! http://www.netlib.no/netlib/lapack/double/dgelsd.f
! http://www.netlib.no/netlib/lapack/double/dgglse.f
!
! For examples:
! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgels_ex.f.htm
!
! For information about LAPACK solving linear methods in general visit:
! http://www.netlib.org/lapack/lug/node26.html
!**********************************************************************************
          !print*, INFO
          molP%a1 = -B(1,NRHS)
          molP%a2 = -B(2,NRHS)
          molP%a3 = -B(3,NRHS)
          if (INFO .ne. 0) then
            call LLS_error(B(1,NRHS),B(2,NRHS),B(3,NRHS),econ)
          endif

  END SUBROUTINE calc_QParam
!--------------------------------------------------------------------------------------------------------------------