c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: /srv/svn/cvs/cvsroot/arts/3rdparty/disort1.2/LINPAK.f,v 1.1 2006/02/21 16:23:28 olemke Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     The following subroutines in this file have been renamed to avoid
c     interference with BLAS/LAPACK libraries:

c          SASUM -> HSASUM
c          SDOT  -> HSDT
c          SAXPY -> HSXPY
c          ISAMAX -> HISMX
c          SSCAL -> HSSCL
c          SSWAP -> HSWP

c
c
c Call tree:
c
c    SGBCO
c       HSASUM
c       HSDT
c       HSXPY
c       SGBFA
c           HISMX
c           HSXPY
c           HSSCL
c       HSSCL
c   SGBSL
c       HSDT
c       HSXPY
c   SGECO
c       HSASUM
c       HSDT
c       HSXPY
c       SGEFA
c           HISMX
c           HSXPY
c           HSSCL
c       HSSCL
c   SGESL
c       HSDT
c       HSXPY
c   HSSWP
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      SUBROUTINE SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

c         Factors a real band matrix by Gaussian elimination
c         and estimates the condition of the matrix.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     If  RCOND  is not needed, SGBFA is slightly faster.
c     To solve  A*X = B , follow SBGCO by SGBSL.

c     input:

C        ABD     REAL(LDA, N)
c                contains the matrix in band storage.  The columns
c                of the matrix are stored in the columns of  ABD  and
c                the diagonals of the matrix are stored in rows
c                ML+1 through 2*ML+MU+1 of  ABD .
c                See the comments below for details.

C        LDA     INTEGER
c                the leading dimension of the array  ABD .
c                LDA must be .GE. 2*ML + MU + 1 .

C        N       INTEGER
c                the order of the original matrix.

C        ML      INTEGER
c                number of diagonals below the main diagonal.
c                0 .LE. ML .LT. N .

C        MU      INTEGER
c                number of diagonals above the main diagonal.
c                0 .LE. MU .LT. N .
c                more efficient if  ML .LE. MU .

c     on return

c        ABD     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                The factorization can be written  A = L*U  where
c                L  is a product of permutation and unit lower
c                triangular matrices and  U  is upper triangular.

C        IPVT    INTEGER(N)
c                an integer vector of pivot indices.

C        RCOND   REAL
c                an estimate of the reciprocal condition of  A .
c                For the system  A*X = B , relative perturbations
c                in  A  and  B  of size  epsilon  may cause
c                relative perturbations in  X  of size  epsilon/RCOND .
c                If  RCOND  is so small that the logical expression
c                           1.0 + RCOND .EQ. 1.0
c                is true, then  A  may be singular to working
c                precision.  In particular,  RCOND  is zero  if
c                exact singularity is detected or the estimate
c                underflows.

C        Z       REAL(N)
c                a work vector whose contents are usually unimportant.
c                If  A  is close to a singular matrix, then  Z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .

c     Band storage

c           If  A  is a band matrix, the following program segment
c           will set up the input.

c                   ML = (band width below the diagonal)
c                   MU = (band width above the diagonal)
c                   M = ML + MU + 1
c                   DO 20 J = 1, N
c                      I1 = MAX(1, J-MU)
c                      I2 = MIN(N, J+ML)
c                      DO 10 I = I1, I2
c                         K = I - J + M
c                         ABD(K,J) = A(I,J)
c                10    CONTINUE
c                20 CONTINUE

c           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
c           In addition, the first  ML  rows in  ABD  are used for
c           elements generated during the triangularization.
c           The total number of rows needed in  ABD  is  2*ML+MU+1 .
c           The  ML+MU by ML+MU  upper left triangle and the
c           ML by ML  lower right triangle are not referenced.

c     Example:  if the original matrix is

c           11 12 13  0  0  0
c           21 22 23 24  0  0
c            0 32 33 34 35  0
c            0  0 43 44 45 46
c            0  0  0 54 55 56
c            0  0  0  0 65 66

c      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain

c            *  *  *  +  +  +  , * = not used
c            *  * 13 24 35 46  , + = used for pivoting
c            * 12 23 34 45 56
c           11 22 33 44 55 66
c           21 32 43 54 65  *

c --------------------------------------------------------------------


c     .. Scalar Arguments ..

      INTEGER   LDA, ML, MU, N
      REAL      RCOND
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * ), Z( * )
c     ..
c     .. Local Scalars ..

      INTEGER   INFO, IS, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
      REAL      ANORM, EK, S, SM, T, WK, WKM, YNORM
c     ..
c     .. External Functions ..

      REAL      HSASUM, HSDT
      EXTERNAL  HSASUM, HSDT
c     ..
c     .. External Subroutines ..

      EXTERNAL  HSXPY, SGBFA, HSSCL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, MIN, SIGN
c     ..


c                       ** compute 1-norm of A
      ANORM  = 0.0E0
      L  = ML + 1
      IS = L + MU

      DO 10 J = 1, N

         ANORM  = MAX( ANORM, HSASUM( L,ABD( IS,J ),1 ) )

         IF( IS.GT.ML + 1 ) IS = IS - 1

         IF( J.LE.MU ) L  = L + 1

         IF( J.GE.N - ML ) L  = L - 1

   10 CONTINUE
c                                               ** factor

      CALL SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

c     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
c     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
c     trans(A) is the transpose of A.  The components of E  are
c     chosen to cause maximum local growth in the elements of W  where
c     trans(U)*W = E.  The vectors are frequently rescaled to avoid
c     overflow.

c                     ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      M  = ML + MU + 1
      JU = 0

      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( EK - Z( K ) )

            CALL HSSCL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( ABD( M,K ).NE.0.0E0 ) THEN

            WK   = WK / ABD( M, K )
            WKM  = WKM / ABD( M, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1
         JU   = MIN( MAX( JU,MU + IPVT( K ) ), N )
         MM   = M

         IF( KP1.LE.JU ) THEN

            DO 30 J = KP1, JU
               MM     = MM - 1
               SM     = SM + ABS( Z( J ) + WKM*ABD( MM,J ) )
               Z( J ) = Z( J ) + WK*ABD( MM, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM
               MM = M

               DO 40 J = KP1, JU
                  MM = MM - 1
                  Z( J ) = Z( J ) + T*ABD( MM, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )

c                         ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB
         LM = MIN( ML, N - K )

         IF( K.LT.N )
     &       Z( K ) = Z( K ) + HSDT( LM, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL HSSCL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

   60 CONTINUE


      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )

      YNORM  = 1.0E0
c                         ** solve L*V = Y
      DO 70 K = 1, N

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
         LM     = MIN( ML, N - K )

         IF( K.LT.N )
     &       CALL HSXPY( LM, T, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z(K) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z(K) )

            CALL HSSCL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

   70 CONTINUE


      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )

      YNORM  = S*YNORM

c                           ** solve  U*Z = W
      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( Z( K ) )

            CALL HSSCL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( ABD( M,K ).NE.0.0E0 ) Z( K ) = Z( K ) / ABD( M, K )
         IF( ABD( M,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         LM = MIN( K, M ) - 1
         LA = M - LM
         LZ = K - LM
         T  = -Z( K )

         CALL HSXPY( LM, T, ABD( LA,K ), 1, Z( LZ ), 1 )

   80 CONTINUE
c                              ** make znorm = 1.0

      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )

      YNORM  = S*YNORM
      IF( ANORM.NE.0.0E0 ) RCOND  = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND  = 0.0E0

      END

      SUBROUTINE SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

c         Factors a real band matrix by elimination.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     SGBFA is usually called by SBGCO, but it can be called
c     directly with a saving in time if  RCOND  is not needed.

c     Input:  same as SGBCO

c     On return:

c        ABD,IPVT    same as SGBCO

c        INFO    INTEGER
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  This is not an error
c                     condition for this subroutine, but it does
c                     indicate that SGBSL will divide by zero if
c                     called.  Use  RCOND  in SBGCO for a reliable
c                     indication of singularity.

c     (see SGBCO for description of band storage mode)

c ----------------------------------------------------------------


c     .. Scalar Arguments ..

      INTEGER   INFO, LDA, ML, MU, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, I0, J, J0, J1, JU, JZ, K, KP1, L, LM, M, MM, NM1
      REAL      T
c     ..
c     .. External Functions ..

      INTEGER   HISMX
      EXTERNAL  HISMX
c     ..
c     .. External Subroutines ..

      EXTERNAL  HSXPY, HSSCL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MAX, MIN
c     ..


      M    = ML + MU + 1
      INFO = 0
c                        ** zero initial fill-in columns
      J0 = MU + 2
      J1 = MIN( N, M ) - 1

      DO 20 JZ = J0, J1

         I0 = M + 1 - JZ

         DO 10 I = I0, ML
            ABD( I, JZ ) = 0.0E0
   10    CONTINUE

   20 CONTINUE

      JZ = J1
      JU = 0
c                       ** Gaussian elimination with partial pivoting
      NM1  = N - 1

      DO 50 K = 1, NM1

         KP1 = K + 1
c                                  ** zero next fill-in column
         JZ = JZ + 1

         IF( JZ.LE.N ) THEN

            DO 30 I = 1, ML
               ABD( I, JZ ) = 0.0E0
   30       CONTINUE

         END IF
c                                  ** find L = pivot index
         LM  = MIN( ML, N - K )
         L   = HISMX( LM + 1, ABD( M, K ), 1 ) + M - 1
         IPVT( K ) = L + K - M

         IF( ABD( L,K ).EQ.0.0E0 ) THEN
c                                      ** zero pivot implies this column
c                                      ** already triangularized
            INFO = K

         ELSE
c                                ** interchange if necessary
            IF( L.NE.M ) THEN

               T           = ABD( L, K )
               ABD( L, K ) = ABD( M, K )
               ABD( M, K ) = T
            END IF
c                                      ** compute multipliers
            T  = - 1.0E0 / ABD( M, K )

            CALL HSSCL( LM, T, ABD( M + 1,K ), 1 )

c                               ** row elimination with column indexing

            JU = MIN( MAX( JU,MU + IPVT( K ) ), N )
            MM = M

            DO 40 J = KP1, JU

               L  = L - 1
               MM = MM - 1
               T  = ABD( L, J )

               IF( L.NE.MM ) THEN

                  ABD( L, J ) = ABD( MM, J )
                  ABD( MM, J ) = T

               END IF

               CALL HSXPY( LM, T, ABD( M+1, K ), 1, ABD( MM+1, J ), 1)

   40       CONTINUE

         END IF

   50 CONTINUE


      IPVT( N ) = N
      IF( ABD( M,N ).EQ.0.0E0 ) INFO = N

      END

      SUBROUTINE SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

c         Solves the real band system
c            A * X = B  or  transpose(A) * X = B
c         using the factors computed by SBGCO or SGBFA.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     Input:

C        ABD     REAL(LDA, N)
c                the output from SBGCO or SGBFA.

C        LDA     INTEGER
c                the leading dimension of the array  ABD .

C        N       INTEGER
c                the order of the original matrix.

C        ML      INTEGER
c                number of diagonals below the main diagonal.

C        MU      INTEGER
c                number of diagonals above the main diagonal.

C        IPVT    INTEGER(N)
c                the pivot vector from SBGCO or SGBFA.

C        B       REAL(N)
c                the right hand side vector.

C        JOB     INTEGER
c                = 0         to solve  A*X = B ,
c                = nonzero   to solve  transpose(A)*X = B

c     On return

c        B       the solution vector  X

c     Error condition

c        A division by zero will occur if the input factor contains a
c        zero on the diagonal.  Technically, this indicates singularity,
c        but it is often caused by improper arguments or improper
c        setting of LDA .  It will not occur if the subroutines are
c        called correctly and if SBGCO has set RCOND .GT. 0.0
c        or SGBFA has set INFO .EQ. 0 .

c     To compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue

c --------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   JOB, LDA, ML, MU, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   K, KB, L, LA, LB, LM, M, NM1
      REAL      T
c     ..
c     .. External Functions ..

      REAL      HSDT
      EXTERNAL  HSDT
c     ..
c     .. External Subroutines ..

      EXTERNAL  HSXPY
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      M   = MU + ML + 1
      NM1 = N - 1

      IF( JOB.EQ.0 ) THEN
c                           ** solve  A * X = B

c                               ** first solve L*Y = B
         IF( ML.NE.0 ) THEN

            DO 10 K = 1, NM1

               LM = MIN( ML, N - K )
               L  = IPVT( K )
               T  = B( L )

               IF( L.NE.K ) THEN

                  B( L ) = B( K )
                  B( K ) = T

               END IF

               CALL HSXPY( LM, T, ABD( M + 1,K ), 1, B( K + 1 ), 1 )

   10       CONTINUE

         END IF

c                           ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / ABD( M, K )
            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = -B( K )

            CALL HSXPY( LM, T, ABD( LA,K ), 1, B( LB ), 1 )

   20    CONTINUE


      ELSE
c                          ** solve  trans(A) * X = B

c                                  ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = HSDT( LM, ABD( LA,K ), 1, B( LB ), 1 )
            B( K ) = ( B( K ) - T ) / ABD( M, K )

   30    CONTINUE

c                                  ** now solve trans(L)*X = Y
         IF( ML.NE.0 ) THEN

            DO 40 KB = 1, NM1

               K      = N - KB
               LM     = MIN( ML, N - K )
               B( K ) = B( K ) + HSDT( LM, ABD( M+1, K ), 1,
     &                                 B( K+1 ), 1 )
               L      = IPVT( K )

               IF( L.NE.K ) THEN

                  T    = B( L )
                  B( L ) = B( K )
                  B( K ) = T

               END IF

   40       CONTINUE

         END IF

      END IF

      END

      SUBROUTINE SGECO( A, LDA, N, IPVT, RCOND, Z )

c         Factors a real matrix by Gaussian elimination
c         and estimates the condition of the matrix.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c         If  RCOND  is not needed, SGEFA is slightly faster.
c         To solve  A*X = B , follow SGECO by SGESL.

c     On entry

c        A       REAL(LDA, N)
c                the matrix to be factored.

c        LDA     INTEGER
c                the leading dimension of the array  A .

c        N       INTEGER
c                the order of the matrix  A .

c     On return

c        A       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                The factorization can be written  A = L*U , where
c                L  is a product of permutation and unit lower
c                triangular matrices and  U  is upper triangular.

c        IPVT    INTEGER(N)
c                an integer vector of pivot indices.

c        RCOND   REAL
c                an estimate of the reciprocal condition of  A .
c                For the system  A*X = B , relative perturbations
c                in  A  and  B  of size  epsilon  may cause
c                relative perturbations in  X  of size  epsilon/RCOND .
c                If  RCOND  is so small that the logical expression
c                           1.0 + RCOND .EQ. 1.0
c                is true, then  A  may be singular to working
c                precision.  In particular,  RCOND  is zero  if
c                exact singularity is detected or the estimate
c                underflows.

C        Z       REAL(N)
c                a work vector whose contents are usually unimportant.
c                If  A  is close to a singular matrix, then  Z  is
c                an approximate null vector in the sense that
c                norm(A*Z) = RCOND*norm(A)*norm(Z) .

c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   LDA, N
      REAL      RCOND
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * ), Z( * )
c     ..
c     .. Local Scalars ..

      INTEGER   INFO, J, K, KB, KP1, L
      REAL      ANORM, EK, S, SM, T, WK, WKM, YNORM
c     ..
c     .. External Functions ..

      REAL      HSASUM, HSDT
      EXTERNAL  HSASUM, HSDT
c     ..
c     .. External Subroutines ..

      EXTERNAL  HSXPY, SGEFA, HSSCL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, SIGN
c     ..


c                        ** compute 1-norm of A
      ANORM  = 0.0E0
      DO 10 J = 1, N
         ANORM  = MAX( ANORM, HSASUM( N,A( 1,J ),1 ) )
   10 CONTINUE
c                                      ** factor

      CALL SGEFA( A, LDA, N, IPVT, INFO )

c     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
c     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E .
c     trans(A) is the transpose of A.  The components of E  are
c     chosen to cause maximum local growth in the elements of W  where
c     trans(U)*W = E.  The vectors are frequently rescaled to avoid
c     overflow.

c                        ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( EK - Z( K ) )

            CALL HSSCL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( A( K,K ).NE.0.0E0 ) THEN

            WK   = WK / A( K, K )
            WKM  = WKM / A( K, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1

         IF( KP1.LE.N ) THEN

            DO 30 J = KP1, N
               SM     = SM + ABS( Z( J ) + WKM*A( K,J ) )
               Z( J ) = Z( J ) + WK*A( K, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM

               DO 40 J = KP1, N
                  Z( J ) = Z( J ) + T*A( K, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )
c                                ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB

         IF( K.LT.N )
     &       Z( K ) = Z( K ) + HSDT( N - K, A( K+1, K ), 1, Z( K+1 ), 1)

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL HSSCL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
   60 CONTINUE


      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )
c                                 ** solve L*V = Y
      YNORM  = 1.0E0

      DO 70 K = 1, N
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

         IF( K.LT.N ) CALL HSXPY( N - K, T, A( K + 1,K ), 1, Z( K + 1 ),
     &                            1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL HSSCL( N, S, Z, 1 )

            YNORM  = S*YNORM
         END IF

   70 CONTINUE


      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )
c                                  ** solve  U*Z = V
      YNORM  = S*YNORM

      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( Z( K ) )

            CALL HSSCL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( A( K,K ).NE.0.0E0 ) Z( K ) = Z( K ) / A( K, K )

         IF( A( K,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         T  = -Z( K )

         CALL HSXPY( K - 1, T, A( 1,K ), 1, Z( 1 ), 1 )

   80 CONTINUE
c                                   ** make znorm = 1.0
      S  = 1.0E0 / HSASUM( N, Z, 1 )

      CALL HSSCL( N, S, Z, 1 )

      YNORM  = S*YNORM

      IF( ANORM.NE.0.0E0 ) RCOND = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND = 0.0E0

      END

      SUBROUTINE SGEFA( A, LDA, N, IPVT, INFO )

c         Factors a real matrix by Gaussian elimination.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     SGEFA is usually called by SGECO, but it can be called
c     directly with a saving in time if  RCOND  is not needed.
c     (time for SGECO) = (1 + 9/N) * (time for SGEFA) .

c     Input:  same as SGECO

c     On return:

c        A,IPVT  same as SGECO

c        INFO    INTEGER
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  This is not an error
c                     condition for this subroutine, but it does
c                     indicate that SGESL or SGEDI will divide by zero
c                     if called.  Use  RCOND  in SGECO for a reliable
c                     indication of singularity.

c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INFO, LDA, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * )
c     ..
c     .. Local Scalars ..

      INTEGER   J, K, KP1, L, NM1
      REAL      T
c     ..
c     .. External Functions ..

      INTEGER   HISMX
      EXTERNAL  HISMX
c     ..
c     .. External Subroutines ..

      EXTERNAL  HSXPY, HSSCL
c     ..


c                      ** Gaussian elimination with partial pivoting
      INFO = 0
      NM1  = N - 1

      DO 20 K = 1, NM1

         KP1  = K + 1
c                                            ** find L = pivot index

         L  = HISMX( N - K + 1, A( K,K ), 1 ) + K - 1
         IPVT( K ) = L

         IF( A( L,K ).EQ.0.0E0 ) THEN
c                                     ** zero pivot implies this column
c                                     ** already triangularized
            INFO = K

         ELSE
c                                     ** interchange if necessary
            IF( L.NE.K ) THEN

               T         = A( L, K )
               A( L, K ) = A( K, K )
               A( K, K ) = T

            END IF
c                                     ** compute multipliers
            T  = -1.0E0 / A( K, K )

            CALL HSSCL( N - K, T, A( K + 1,K ), 1 )

c                              ** row elimination with column indexing
            DO 10 J = KP1, N

               T  = A( L, J )

               IF( L.NE.K ) THEN

                  A( L, J ) = A( K, J )
                  A( K, J ) = T

               END IF

               CALL HSXPY( N-K, T, A( K+1, K ), 1, A( K+1, J ), 1 )

   10       CONTINUE

         END IF

   20 CONTINUE


      IPVT( N ) = N
      IF( A( N,N ) .EQ. 0.0E0 ) INFO = N

      END

      SUBROUTINE SGESL( A, LDA, N, IPVT, B, JOB )

c         Solves the real system
c            A * X = B  or  transpose(A) * X = B
c         using the factors computed by SGECO or SGEFA.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     On entry

c        A       REAL(LDA, N)
c                the output from SGECO or SGEFA.

c        LDA     INTEGER
c                the leading dimension of the array  A

c        N       INTEGER
c                the order of the matrix  A

c        IPVT    INTEGER(N)
c                the pivot vector from SGECO or SGEFA.

c        B       REAL(N)
c                the right hand side vector.

c        JOB     INTEGER
c                = 0         to solve  A*X = B ,
c                = nonzero   to solve  transpose(A)*X = B

c     On return

c        B       the solution vector  X

c     Error condition

c        A division by zero will occur if the input factor contains a
c        zero on the diagonal.  Technically, this indicates singularity,
c        but it is often caused by improper arguments or improper
c        setting of LDA.  It will not occur if the subroutines are
c        called correctly and if SGECO has set RCOND .GT. 0.0
c        or SGEFA has set INFO .EQ. 0 .

c     To compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue

c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   JOB, LDA, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   K, KB, L, NM1
      REAL      T
c     ..
c     .. External Functions ..

      REAL      HSDT
      EXTERNAL  HSDT
c     ..
c     .. External Subroutines ..

      EXTERNAL  HSXPY
c     ..


      NM1  = N - 1

      IF( JOB.EQ.0 ) THEN
c                                 ** solve  A * X = B

c                                     ** first solve  L*Y = B
         DO 10 K = 1, NM1

            L  = IPVT( K )
            T  = B( L )

            IF( L.NE.K ) THEN

               B( L ) = B( K )
               B( K ) = T

            END IF

            CALL HSXPY( N - K, T, A( K+1, K ), 1, B( K+1 ), 1 )

   10    CONTINUE
c                                    ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / A( K, K )
            T      = - B( K )

            CALL HSXPY( K-1, T, A( 1, K ), 1, B(1), 1 )

   20    CONTINUE


      ELSE
c                         ** solve  trans(A) * X = B

c                                    ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            T      = HSDT( K - 1, A( 1,K ), 1, B( 1 ), 1 )
            B( K ) = ( B( K ) - T ) / A( K, K )

   30    CONTINUE

c                                    ** now solve  trans(l)*x = y
         DO 40 KB = 1, NM1

            K      = N - KB
            B( K ) = B( K ) + HSDT( N - K, A( K+1, K ), 1, B( K+1 ), 1)
            L      = IPVT( K )

            IF( L.NE.K ) THEN

               T      = B( L )
               B( L ) = B( K )
               B( K ) = T

            END IF

   40    CONTINUE

      END IF

      END

      REAL FUNCTION HSASUM( N, SX, INCX )

c  INPUT--    N  Number of elements in vector to be summed
c            SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
c          INCX  Spacing of vector elements in SX

c  OUTPUT-- HSASUM   Sum from 0 to N-1 of  ABS(SX(1+I*INCX))
c ----------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MOD
c     ..

      HSASUM  = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN
c                                          ** non-unit increments
         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            HSASUM  = HSASUM + ABS( SX( I ) )
   10    CONTINUE

      ELSE
c                                          ** unit increments
         M  = MOD( N, 6 )

         IF( M.NE.0 ) THEN
c                             ** clean-up loop so remaining vector
c                             ** length is a multiple of 6.
            DO 20 I = 1, M
               HSASUM  = HSASUM + ABS( SX( I ) )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 6
            HSASUM  = HSASUM + ABS( SX( I ) ) + ABS( SX( I + 1 ) ) +
     &               ABS( SX( I + 2 ) ) + ABS( SX( I + 3 ) ) +
     &               ABS( SX( I + 4 ) ) + ABS( SX( I + 5 ) )
   30    CONTINUE

      END IF

      END

      SUBROUTINE HSXPY( N, SA, SX, INCX, SY, INCY )

c          Y = A*X + Y  (X, Y = vectors, A = scalar)

c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SA  Single precision scalar multiplier A
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c       SY   For I = 0 to N-1, overwrite  SY(LY+I*INCY) with
c                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
c            where LX = 1          if INCX .GE. 0,
c                     = (-INCX)*N  if INCX .LT. 0
c            and LY is defined analogously using INCY.
c ------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
      REAL      SA
c     ..
c     .. Array Arguments ..

      REAL      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SY( I ) = SY( I ) + SA*SX( I )
   10    CONTINUE

      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 4 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 4.
            DO 20 I = 1, M
               SY( I ) = SY( I ) + SA*SX( I )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 4
            SY( I ) = SY( I ) + SA*SX( I )
            SY( I + 1 ) = SY( I + 1 ) + SA*SX( I + 1 )
            SY( I + 2 ) = SY( I + 2 ) + SA*SX( I + 2 )
            SY( I + 3 ) = SY( I + 3 ) + SA*SX( I + 3 )
   30    CONTINUE


      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1
         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            SY( IY ) = SY( IY ) + SA*SX( IX )
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE

      END IF

      END

      REAL FUNCTION HSDT( N, SX, INCX, SY, INCY )

c        Single-prec dot product of vectors  X  and  Y

c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c     HSDT   Sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
c            where  LX = 1          if INCX .GE. 0,
c                      = (-INCX)*N  if INCX .LT. 0,
c            and LY is defined analogously using INCY.
c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      HSDT = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            HSDT = HSDT + SX( I )*SY( I )
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 4.
            DO 20 I = 1, M
               HSDT = HSDT + SX( I )*SY( I )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            HSDT = HSDT + SX( I )*SY( I ) + SX( I + 1 )*SY( I + 1 ) +
     &               SX( I + 2 )*SY( I + 2 ) + SX( I + 3 )*SY( I + 3 ) +
     &               SX( I + 4 )*SY( I + 4 )
   30    CONTINUE

      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            HSDT = HSDT + SX( IX )*SY( IY )
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      SUBROUTINE HSSCL( N, SA, SX, INCX )

c         Multiply vector SX by scalar SA

c  INPUT--  N  Number of elements in vector
c          SA  Single precision scale factor
c          SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
c        INCX  Spacing of vector elements in SX

c OUTPUT-- SX  Replace  SX(1+I*INCX)  with  SA * SX(1+I*INCX)
c                for I = 0 to N-1
c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, N
      REAL      SA
c     ..
c     .. Array Arguments ..

      REAL      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SX( I ) = SA*SX( I )
   10    CONTINUE


      ELSE

         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
c                           ** clean-up loop so remaining vector length
c                           ** is a multiple of 5.
            DO 20 I = 1, M
               SX( I ) = SA*SX( I )
   20       CONTINUE

         END IF
c                             ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            SX( I ) = SA*SX( I )
            SX( I + 1 ) = SA*SX( I + 1 )
            SX( I + 2 ) = SA*SX( I + 2 )
            SX( I + 3 ) = SA*SX( I + 3 )
            SX( I + 4 ) = SA*SX( I + 4 )
   30    CONTINUE

      END IF

      END

      SUBROUTINE HSSWP( N, SX, INCX, SY, INCY )

c          Interchange s.p vectors  X  and  Y, as follows:

c     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
c     where LX = 1          if INCX .GE. 0,
c              = (-INCX)*N  if INCX .LT. 0
c     and LY is defined analogously using INCY.


c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c       SX  Input vector SY (unchanged if N .LE. 0)
c       SY  Input vector SX (unchanged IF N .LE. 0)
c --------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
      REAL      STEMP1, STEMP2, STEMP3
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX
            STEMP1 = SX( I )
            SX( I ) = SY( I )
            SY( I ) = STEMP1
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 3 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 3.
            DO 20 I = 1, M
               STEMP1 = SX( I )
               SX( I ) = SY( I )
               SY( I ) = STEMP1
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 3
            STEMP1 = SX( I )
            STEMP2 = SX( I + 1 )
            STEMP3 = SX( I + 2 )
            SX( I ) = SY( I )
            SX( I + 1 ) = SY( I + 1 )
            SX( I + 2 ) = SY( I + 2 )
            SY( I ) = STEMP1
            SY( I + 1 ) = STEMP2
            SY( I + 2 ) = STEMP3
   30    CONTINUE


      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            STEMP1 = SX( IX )
            SX( IX ) = SY( IY )
            SY( IY ) = STEMP1
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      INTEGER FUNCTION HISMX( N, SX, INCX )

c INPUT--  N     Number of elements in vector of interest
c          SX    Sing-prec array, length 1+(N-1)*INCX, containing vector
c          INCX  Spacing of vector elements in SX

c OUTPUT-- HISMX   First I, I = 1 to N, to maximize
c                         ABS(SX(1+(I-1)*INCX))
c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, II
      REAL      SMAX, XMAG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS
c     ..


      IF( N.LE.0 ) THEN

         HISMX = 0

      ELSE IF( N.EQ.1 ) THEN

         HISMX = 1

      ELSE

         SMAX = 0.0
         II   = 1

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX

            XMAG = ABS( SX( I ) )

            IF( SMAX.LT.XMAG ) THEN

               SMAX   = XMAG
               HISMX = II

            END IF

            II = II + 1

   10    CONTINUE

      END IF

      END

