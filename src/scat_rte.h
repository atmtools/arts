/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.
*/   

   /**
   * @file   scat_rte.h
   * @author Claudia Emde <claudia@sat.physik.uni-bremen.de>
   * @date   Thu Nov 21 12:43:28 2002
   * 
   * @brief  Radiative transfer calculation
   * 
   * 
   */



#ifndef scatrte_h
#define scatrte_h


void
stokes_vecGeneral(//WS Output and Input:
               Vector& stokes_vec,
               //WS Input:
               const Matrix& ext_mat,
               const Vector& abs_vec,
               const Vector& sca_vec,
               const Numeric& l_step,
               const Numeric& a_planck_value,
               const Index& stokes_dim);

void
stokes_vecScalar(//WS Input and Output:
	      Vector& stokes_vec,
	      //WS Input: 
	      const Matrix& ext_mat,
	      const Vector& abs_vec,
	      const Vector& sca_vec,
	      const Numeric& l_step,
	      const Numeric& a_planck_value,
	      const Index& stokes_dim);

#endif    // scatrte_h
