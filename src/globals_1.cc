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

/*-----------------------------------------------------------------------
FILE:      globals_1.cc

INCLUDES:  This file contains all global variable definitions that do
	   NOT depend on the automatically generated header file
	   wsv.h. It is necessary to have these in a separate file for
	   compiler technical reasons. (With g++-2.95.2 it does not
	   work to declare them constant in the same file where they
	   are defined.)

FUNCTIONS: None

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "workspace.h"
#include "absorption.h"

//                     -----------------
//--------------------< Workspace Stuff >--------------------
//                     -----------------

/** The workspace itself. */
WorkSpace workspace;

/** The names associated with Wsv groups as strings. 
    
    \verbatim
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Must be consistent with the enum in workspace.h! Later on,
    these enums could also be generated automatically, but that would
    have to take place before the wsv_data is defined, since that
    needs these enums.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \endverbatim  */
//ARRAY<string> wsv_group_names;
// Moved to groups.cc
// FIXME: REMOVE this chunk

/** The lookup information for the workspace variables. */
ARRAY<WsvRecord> wsv_data;

/** The map assiciated with wsv_data. */
std::map<string, size_t> WsvMap;


//                     ------------------
//--------------------< Absorption Stuff >--------------------
//                     ------------------

/** The lookup information for all the different species. */
ARRAY<SpeciesRecord> species_data;

/** The map associated with species_data. */
std::map<string, size_t> SpeciesMap;

