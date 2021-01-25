/* Copyright (C) 2021
      Jon Petersen <jon.petersen@studium.uni-hamburg.de>
      Manfred Brath  <manfred.brath@uni-hamburg.de>

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

/*!
  \file   star.cc

  \brief  This file contains functions that are adapted from TESSEM2
  code which is used to calculate surface emissivity.

  The implementation is based on the Fortran code v1.0 developed by Catherine
  Prigent and Filipe Aires in the EUMETSAT Study on surface emissivity at
  microwave and sub-millimeter frequencies project EUM/CO/14/4600001473/CJA
*/

#include "star.h"
