      module ffastwigxj
      implicit none

      interface

C     Initialisation

      subroutine ffastwigxj_load(filename, type)
      character (len=*) filename
      integer*4 type
      end subroutine

C     Release
      
      subroutine ffastwigxj_unload(type)
      integer*4 type
      end subroutine

C     3j
      
      subroutine ffw3jj_canon(two_jv, rx)
      integer*4 two_jv(5)
      integer*8 rx
      end subroutine

      subroutine ffw3jj_prefetch(x)
      integer*8 x
      end subroutine

      function ffw3jj_get(x)
      real*8 ffw3jj_get
      integer*8 x
      end function

      function ffw3jjl(two_jv)
      real*8 ffw3jja
      integer*4 two_jv(5)
      end function

      function ffw3jja(two_j1, two_j2, two_j3,
     *                 two_m1, two_m2)
      real*8 ffw3jja
      integer*4 two_j1, two_j2, two_j3, two_m1, two_m2
      end function

      function ffw3jja6(two_j1, two_j2, two_j3,
     *                  two_m1, two_m2, two_m3)
      real*8 ffw3jja6
      integer*4 two_j1, two_j2, two_j3, two_m1, two_m2, two_m3
      end function

C     6j

      subroutine ffw6jj_canon(two_jv, rx)
      integer*4 two_jv(6)
      integer*8 rx
      end subroutine

      subroutine ffw6jj_prefetch(x)
      integer*8 x
      end subroutine

      function ffw6jj_get(x)
      real*8 ffw6jj_get
      integer*8 x
      end function

      function ffw6jjl(two_jv)
      real*8 ffw6jja
      integer*4 two_jv(6)
      end function

      function ffw6jja(two_j1, two_j2, two_j3,
     *                 two_j4, two_j5, two_j6)
      real*8 ffw6jja
      integer*4 two_j1, two_j2, two_j3, two_j4, two_j5, two_j6
      end function

C     9j      

      subroutine ffw9jj_canon(two_jv, rx)
      integer*4 two_jv(9)
      integer*8 rx
      end subroutine

      subroutine ffw9jj_prefetch(x)
      integer*8 x
      end subroutine

      function ffw9jj_get(x)
      real*8 ffw9jj_get
      integer*8 x
      end function

      function ffw9jjl(two_jv)
      real*8 ffw9jja
      integer*4 two_jv(9)
      end function

      function ffw9jja(two_j1, two_j2, two_j3,
     *                 two_j4, two_j5, two_j6,
     *                 two_j7, two_j8, two_j9)
      real*8 ffw9jja
      integer*4 two_j1, two_j2, two_j3, two_j4, two_j5, two_j6,
     *          two_j7, two_j8, two_j9
      end function

C     -
      
      end interface

      end module ffastwigxj
