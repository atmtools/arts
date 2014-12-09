MODULE PARKIND1
!
!     *** Define usual kinds for strong typing ***
!     Option of Kinds modified for NEC interface Roger Saunders
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIT = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: JPIS = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)
INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)       !Standard integer type
!INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(15)     !Required for METO IBM OPS code (not standalone)

!Special integer type to be used for sensitive adress calculations
!should be *8 for a machine with 8byte adressing for optimum performance
!ifdef ADDRESS64
INTEGER, PARAMETER :: JPIA = JPIB
!#else
!INTEGER, PARAMETER :: JPIA = JPIM
!#endif
!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER :: JPRT = SELECTED_REAL_KIND(2,1)
INTEGER, PARAMETER :: JPRS = SELECTED_REAL_KIND(4,2)
INTEGER, PARAMETER :: JPRM = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300) !Standard real type
!
!     Logical Kinds
!     -------------
INTEGER, PARAMETER :: JPLM = KIND(.TRUE.)               !Standard logical type
!
END MODULE PARKIND1
