
C     Copyright 2015 Haakan T. Johansson 
C	
C     This file is part of WIGXJPF.
C	
C     WIGXJPF is free software: you can redistribute it and/or modify it
C     under the terms of the GNU Lesser General Public License as
C     published by the Free Software Foundation, either version 3 of the
C     License, or (at your option) any later version.
C	
C     WIGXJPF is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C     Lesser General Public License for more details.
C	
C     You should have received a copy of the GNU Lesser General Public
C     License along with WIGXJPF.  If not, see
C     <http://www.gnu.org/licenses/>.
C     

C     Compile with: gfortran fsimple.f -o fsimple -L../lib/ -lwigxjpf

      program fsimple

C     Module with function interface declarations.
      use fwigxjpf
      
      implicit none

C     Function prototypes (not recommended, use only when you cannot use
C     the module fwigxjpf)
C     real*8 fwig3jj, fwig6jj, fwig9jj

C     Local variables
      real*8 val3j, val6j, val9j

      write (*,"(A)") "WIGXJPF fortran test program"

      call fwig_table_init(2*100,9)
      call fwig_temp_init(2*100)

C     Note that arguments are in two_j = 2*j.

      val3j = fwig3jj(2* 10 , 2* 15 , 2* 10 ,
     c                2*(-3), 2* 12 , 2*(-9))

      write (*,"(A40,F25.15)") "3J(10  15  10; -3  12  -9):", val3j;

      val6j = fwig6jj(2* 10 , 2* 15 , 2* 10 ,
     c                2*  7,  2*  7 , 2*  9 )

      write (*,"(A40,F25.15)") "6J{10  15  10;  7   7   9}:", val6j;

      val9j = fwig9jj( 1,  2,  3,
     c                 4,  6,  8,
     c                 3,  6,  9)

      write (*,"(A40,F25.15)")
     c     "9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}:", val9j;

      call fwig_temp_free();
      call fwig_table_free();

      end program fsimple
