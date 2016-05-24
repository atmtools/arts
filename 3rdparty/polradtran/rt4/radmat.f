
C      This file contains a number of simple matrix manipulation routines.
C      All of the routines operate on REAL*8 matrices.
C      The size is passed as the first parameters:  N rows by M columns.
C      If there is no M parameter then a N by N square matrix is assumed.
C
C    CALL MCOPY (N,M, X, Y)        Copies matrix: Y = X
C
C    CALL MADD (N,M, X, Y, Z)      Adds matrices: Z = X + Y
C                                    Z can be X or Y
C    CALL MSUB (N,M, X, Y, Z)      Subtracts matrices: Z = X - Y
C                                    Z can be X or Y
C    CALL MSCALARMULT (N,M, C, X, Y) Scalar multiply: Y = C*X
C                                    C is real scalar; Y can be X
C    CALL MZERO (N,M, X)           Zeros all elements in X
C
C    CALL MDIAG (N, V, X)          Formats a vector (V) into a diagonal
C                                    matrix (X)
C    CALL MIDENTITY (N, X)         Make identity matrix in X
C                                                   t
C    CALL MTRANSPOSE (N,M, X, Y)   Transposes: Y = X ;   Y cannot be X
C
C    CALL MMULT (N,M,L, X, Y, Z)   Matrix multiply: Z = X*Y
C                                    X is N by M and Y is M by L.
C                                    Z cannot be X or Y
C                                                          -1
C    CALL MINVERT (N, X, Y)        Matrix inversion:  Y = X
C				     X gets LU decomposition.




      SUBROUTINE MCOPY (N, M, MATRIX1, MATRIX2)
      INTEGER  N, M, I
      REAL*8   MATRIX1(1), MATRIX2(1)

      DO 100 I = 1, N*M
          MATRIX2(I) = MATRIX1(I)
100   CONTINUE

      RETURN
      END



      SUBROUTINE MADD (N, M, MATRIX1, MATRIX2, MATRIX3)
      INTEGER  N, M, I
      REAL*8   MATRIX1(1), MATRIX2(1),  MATRIX3(1)

      DO 100 I = 1, N*M
          MATRIX3(I) = MATRIX1(I) + MATRIX2(I)
100   CONTINUE

      RETURN
      END



      SUBROUTINE MSUB (N, M, MATRIX1, MATRIX2, MATRIX3)
      INTEGER  N, M, I
      REAL*8   MATRIX1(1), MATRIX2(1),  MATRIX3(1)

      DO 100 I = 1, N*M
          MATRIX3(I) = MATRIX1(I) - MATRIX2(I)
100   CONTINUE

      RETURN
      END




      SUBROUTINE MSCALARMULT (N, M, C, MATRIX1, MATRIX2)
      INTEGER  N, M, I
      REAL*8   C, MATRIX1(1), MATRIX2(1)

      DO 100 I = 1, N*M
          MATRIX2(I) = C*MATRIX1(I)
100   CONTINUE

      RETURN
      END


      SUBROUTINE MZERO (N, M, MATRIX1)
      INTEGER  N, M, I
      REAL*8   MATRIX1(1)

      DO 100 I = 1, N*M
          MATRIX1(I) = 0.0D0
100   CONTINUE

      RETURN
      END



      SUBROUTINE MDIAG (N,  VECTOR, MATRIX)
      INTEGER  N,  I, J
      REAL*8   VECTOR(1), MATRIX(N,N)

      DO 110 I = 1, N
        DO 100 J = 1, N
          MATRIX(I,J) = 0.0
100     CONTINUE
        MATRIX(I,I) = VECTOR(I)
110   CONTINUE

      RETURN
      END



      SUBROUTINE MIDENTITY (N, MATRIX)
      INTEGER  N, I, J
      REAL*8   MATRIX(N,N)

      DO 110 I = 1, N
        DO 100 J = 1, N
          MATRIX(I,J) = 0.0
100     CONTINUE
        MATRIX(I,I) = 1.0
110   CONTINUE

      RETURN
      END



      SUBROUTINE MTRANSPOSE (N, M, MATRIX1, MATRIX2)
      INTEGER  N, M, I, J
      REAL*8   MATRIX1(N,M), MATRIX2(M,N)

      DO 100 I = 1, N
        DO 100 J = 1, M
             MATRIX2(I,J) = MATRIX1(J,I)
100     CONTINUE

      RETURN
      END



      SUBROUTINE MMULT (N, M, L, MATRIX1, MATRIX2, MATRIX3)
      INTEGER  N, M, L,  I, J, K
      REAL*8   MATRIX1(N,M), MATRIX2(M,L), MATRIX3(N,L),   SUM

      DO 200 I = 1, N
        DO 200 J = 1, L
          SUM = 0.0
          DO 100 K = 1, M
            SUM = SUM + MATRIX1(I,K)*MATRIX2(K,J)
100       CONTINUE
          MATRIX3(I,J) = SUM
200     CONTINUE
      RETURN
      END




      SUBROUTINE MINVERT (N, MATRIX1, MATRIX2)
      INTEGER  N
      REAL*8   MATRIX1(N,N), MATRIX2(N,N)
      INTEGER  NMAX
      PARAMETER (NMAX=256)
      INTEGER  I, J, INDX(NMAX), IZ
      REAL*8   DET(2), WORK(NMAX)

      IF (N .GT. NMAX) THEN
          WRITE (*,'(1X,A,I3)')
     .    'Exceeded maximum matrix size for inversion.  Max = ', NMAX
          STOP
      ENDIF
      CALL DGEFA (MATRIX1, N, N, INDX, IZ)
      IF (IZ .GT. 0) THEN
	  WRITE (*,'(1X,A,I3)')
     .    'Encountered a zero pivot at element ', IZ
	  STOP
      ENDIF
      DO 100 I = 1, N
	DO 110 J = 1, N
           MATRIX2(I,J) = MATRIX1(I,J)
110     CONTINUE
100   CONTINUE
      CALL DGEDI (MATRIX2, N, N, INDX, DET, WORK, 1 )

      RETURN
      END



      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end


      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end


      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end


      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end


      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end


      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
