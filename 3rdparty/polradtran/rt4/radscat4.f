
      SUBROUTINE GET_SCAT_FILE (NSTOKES, NUMMU,
     .                          QUAD_TYPE, SCAT_FILE,
     .                          SCATTER_MATRIX,
     .                          EXTINCT_MATRIX, EMIS_VECTOR)
C       Reads in a oriented particle scattering file containing the
C       Mueller scattering matrix, extinction matrix, and emission
C       vector for a set of quadrature angles.
C       Only reads in the m=0 Fourier azimuthal angle mode.
      INTEGER NSTOKES, NUMMU
      REAL*8  SCATTER_MATRIX(NSTOKES,NUMMU,NSTOKES,NUMMU,4)
      REAL*8  EXTINCT_MATRIX(4,4,NUMMU,2)
      REAL*8  EMIS_VECTOR(4,NUMMU,2)
      CHARACTER*(*)  SCAT_FILE, QUAD_TYPE
      INTEGER NMU, NAZ, J, L, I, I1, I2, J1, J2, L1, L2, M
      REAL*8  MU, MU1, MU2
      CHARACTER*(16)  FQUAD
      CHARACTER*132   BUFFER


      OPEN (UNIT=2, FILE=SCAT_FILE, STATUS='OLD')
100   CONTINUE
          READ (2,'(A)') BUFFER
      IF (BUFFER(1:1) .EQ. 'C') GOTO 100
      BACKSPACE 2

      READ (2,*) NMU, NAZ, FQUAD
      IF (NMU .NE. NUMMU .OR. FQUAD(1:1) .NE. QUAD_TYPE(1:1)) THEN
          WRITE (*,*) 'Scattering file incompatible with quadrature.'
          STOP
      ENDIF


      READ (2,*)

      DO L1 = 1, 2
        DO J1 = 1, NUMMU
          DO L2 = 1, 2
            L = 2*(L2-1)+L1
            DO J2 = 1, NUMMU

              READ (2,*) MU1, MU2, M
              IF (M .NE. 0) STOP 'Error reading scattering file.'
              DO I2 = 1, NSTOKES
                READ (2,*) (SCATTER_MATRIX(I2,J2,I1,J1,L), I1=1,NSTOKES)
              ENDDO
              DO I2 = NSTOKES+1,4
                READ (2,*)
              ENDDO
              DO I = 1, 5*(2*NAZ)
                READ (2,*)
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO


      READ (2,*)
      DO L = 1, 2
        DO J = 1, NMU
          READ (2,*)
          DO I2 = 1, NSTOKES
            READ (2,*) (EXTINCT_MATRIX(I2,I1,J,L), I1=1,NSTOKES)
          ENDDO
          DO I2 = NSTOKES+1,4
            READ (2,*)
          ENDDO
        ENDDO
      ENDDO

      READ (2,*)
      DO L = 1, 2
        DO J = 1, NUMMU
          READ (2,*) MU, (EMIS_VECTOR(I,J,L), I=1,NSTOKES)
        ENDDO
      ENDDO

      CLOSE (2)
      RETURN
      END






      SUBROUTINE CHECK_NORM (NSTOKES, NUMMU, QUAD_WEIGHTS,
     .                       SCATTER_MATRIX,
     .                       EXTINCT_MATRIX, EMIS_VECTOR)
      INTEGER  NSTOKES, NUMMU
      REAL*8   QUAD_WEIGHTS(NUMMU)
      REAL*8   SCATTER_MATRIX(NSTOKES,NUMMU,NSTOKES,NUMMU,4)
      REAL*8   EXTINCT_MATRIX(4,4,NUMMU,2)
      REAL*8   EMIS_VECTOR(4,NUMMU,2)
      INTEGER  J1, J2, L
      REAL*8   SUM, MAXSUM, PI
      PARAMETER (PI = 3.1415926535897932384D0)

      MAXSUM = 0.0D0
      DO J1 = 1, NUMMU
        DO L = 1, 2
          SUM = EMIS_VECTOR(1,J1,L) - EXTINCT_MATRIX(1,1,J1,L)
          DO J2 = 1, NUMMU
              SUM = SUM + 2.0D0*PI*QUAD_WEIGHTS(J2)*
     .               ( SCATTER_MATRIX(1,J2,1,J1, L)
     .               + SCATTER_MATRIX(1,J2,1,J1, L+2) )
          ENDDO
          MAXSUM = DMAX1 (MAXSUM, DABS(SUM))
        ENDDO
      ENDDO
      IF (MAXSUM .GT. 1.0D-6) THEN
          WRITE (*,*) 'Scattering function not normalized:', MAXSUM
      ENDIF
      RETURN
      END
