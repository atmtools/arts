      module fwigxjpf
      implicit none

      interface

C     Initialisation

      subroutine fwig_table_init(max_two_j, wigner_type)
      integer*4 max_two_j, wigner_type
      end subroutine

      subroutine fwig_temp_init(max_two_j)
      integer*4 max_two_j
      end subroutine

      subroutine fwig_thread_temp_init(max_two_j)
      integer*4 max_two_j
      end subroutine

C     Release

      subroutine fwig_table_free()
      end subroutine

      subroutine fwig_temp_free()
      end subroutine

C     3j

      function fwig3jj(two_j1, two_j2, two_j3,
     *                 two_m1, two_m2, two_m3)
      real*8 fwig3jj
      integer*4 two_j1, two_j2, two_j3, two_m1, two_m2, two_m3
      end function

C     6j

      function fwig6jj(two_j1, two_j2, two_j3,
     *                 two_j4, two_j5, two_j6)
      real*8 fwig6jj
      integer*4 two_j1, two_j2, two_j3, two_j4, two_j5, two_j6
      end function

C     9j

      function fwig9jj(two_j1, two_j2, two_j3,
     *                 two_j4, two_j5, two_j6,
     *                 two_j7, two_j8, two_j9)
      real*8 fwig9jj
      integer*4 two_j1, two_j2, two_j3, two_j4, two_j5, two_j6,
     *          two_j7, two_j8, two_j9
      end function

C     -

      end interface

      end module fwigxjpf
