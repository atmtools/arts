c
c...date last changed 25 Nov 2001
c  ****************************************
      SUBROUTINE ISO_82_to_85 (
     I  MOLEC, 
     O NSO82, ISO85)
c  ****************************************
c
      implicit DOUBLE PRECISION (a-h,o-z)
c
      INCLUDE 'SPECIES_2011.CMN'
      INCLUDE 'ISOTOPS.CMN'
c
c
         
         WRITE (*,'(A)')
     +   '     ENTER one isotope only'
   30    WRITE (*,'(A)') '     The choices are:'
         WRITE (*,'(3X,8(I6,5X))') (iso82(MOLEC,INS),INS=1,ISONM(MOLEC))
         READ (*,*) NSO82
c
         IF (NSO82 .eq. 0) THEN
           write(*,*) '  Zero is not valid'
             go to 30
           RETURN
         END IF
c
         DO 40 I = 1,ISONM(MOLEC)
          IF (NSO82 .EQ. iso82(MOLEC,I)) THEN
             ISO85 = I
             RETURN
          END IF
   40    CONTINUE
c
         WRITE (*,'(15x,A)') '     Incorrect choice, try again'
         GO TO 30
c
      RETURN
      END
