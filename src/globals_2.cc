/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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
   USA. */

/*!
  \file   globals_2.cc
  \brief  Global variable definitions that 
	  depend on the automatically generated header file wsv.h

  This file contains all global variable definitions that DO depend on
  the automatically generated header file wsv.h. It is necessary to
  have these in a separate file for compiler technical reasons. (With
  g++-2.95.2 it does not work to declare them constant in the same
  file where they are defined.)
  
  Maybe this file should be removed in the future. There is not much
  stuff here.  

  \author Stefan Buehler
  \date 2000-06-10 */


#include <map>
#include "arts.h"
#include "array.h"
#include "auto_wsv.h"
#include "auto_wsv_groups.h"
#include "wsv_aux.h"
#include "methods.h"

/** The workspace itself. */
WorkSpace workspace;


//                     ---------------
//--------------------< Methods Stuff >--------------------
//                     ---------------

/** The lookup information for the workspace methods. */
Array<MdRecord> md_data;

/** The map associated with md_data. */
std::map<String, Index> MdMap;


