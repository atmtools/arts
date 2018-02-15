
C     Copyright 2015 Haakan T. Johansson 
C	
C     This file is part of FASTWIGXJ.
C	
C     FASTWIGXJ is free software: you can redistribute it and/or modify it
C     under the terms of the GNU Lesser General Public License as
C     published by the Free Software Foundation, either version 3 of the
C     License, or (at your option) any later version.
C	
C     FASTWIGXJ is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C     Lesser General Public License for more details.
C	
C     You should have received a copy of the GNU Lesser General Public
C     License along with FASTWIGXJ.  If not, see
C     <http://www.gnu.org/licenses/>.
C     

C     Compile with: gfortran example/flookup.f -o flookup -Imod \
C       -I<path-to-wigxjpf>/mod -Llib -lfastwigxj \
C       -L<path-to-wigxjpf>/lib -lwigxjpf -lwigxjpf_quadmath

      program flookup

C     Module with function interface declarations.
      use ffastwigxj
      use fwigxjpf

C     Function prototypes (not recommended, use only when you cannot use
C     the module ffastwigxj)
C     real*8 ffw3jj_get, ffw3jjl, ffw3jja, ffw3jja6
C     real*8 ffw6jj_get, ffw6jjl, ffw6jja
C     real*8 ffw9jj_get, ffw9jjl, ffw9jja

      implicit none

C     Local variables
      real*8 val3j, val6j, val9j

      write (*,"(A)") "FASTWIGXJ fortran test program"

C     Load tables produced during build test.

      call ffastwigxj_load("test_table_18.3j",      3);
      call ffastwigxj_load("test_table_8.6j",       6);
      call ffastwigxj_load("test_hashed_3_6_12.9j", 9);

C     For fallback to WIGXJPF when symbols too large for table.

      call fwig_table_init(2*100,9)
      call fwig_temp_init(2*100)

C     Note that arguments are in two_j = 2*j.

      val3j = ffw3jja6(2*  5 , 2*  7 , 2*  5 ,
     c                 2*(-3), 2*  5 , 2*(-2))

      write (*,"(A40,F25.15)") "3J( 5   7   5; -3   5  -2):", val3j;

      val3j = ffw3jja6(2* 10 , 2* 15 , 2* 10 ,
     c                 2*(-3), 2* 12 , 2*(-9))

      write (*,"(A40,F25.15)") "3J(10  15  10; -3  12  -9):", val3j;

      val6j = ffw6jja(2*  3 , 2*  4 , 2*  2 ,
     c                2*  2,  2*  2 , 2*  3 )

      write (*,"(A40,F25.15)") "6J{ 3   4   2;  2   2   3}:", val6j;

      val6j = ffw6jja(2* 10 , 2* 15 , 2* 10 ,
     c                2*  7,  2*  7 , 2*  9 )

      write (*,"(A40,F25.15)") "6J{10  15  10;  7   7   9}:", val6j;

      val9j = ffw9jja( 1,  2,  3,
     c                 2,  3,  5,
     c                 3,  3,  6)

      write (*,"(A40,F25.15)")
     c     "9J{0.5 1 1.5; 1 1.5 2.5; 1.5 1.5 3}:", val9j;

      val9j = ffw9jja( 1,  2,  3,
     c                 4,  6,  8,
     c                 3,  6,  9)

      write (*,"(A40,F25.15)")
     c     "9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}:", val9j;

C     Print statistics.

      call ffastwigxj_print_stats();

C     Remove tables from memory.
      
      call ffastwigxj_unload(3);
      call ffastwigxj_unload(6);
      call ffastwigxj_unload(9);
      
C     Release WIGXJPF memory.

      call fwig_temp_free();
      call fwig_table_free();

      end program flookup
