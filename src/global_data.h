/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \file   globals_data.h
  \brief  Global variable declarations

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef global_data_h
#define global_data_h

#include <map>
#include "array.h"
#include "methods.h"
#include "agenda_record.h"


// Needed for global_data::species_data
class SpeciesRecord;
// Needed for global_data::lineshape_data
class LineshapeRecord;
// Needed for global_data::lineshape_norm_data
class LineshapeNormRecord;

namespace global_data {

//                     ---------------
//--------------------<     Methods   >--------------------
//                     ---------------

//! Lookup information for workspace methods.
/*!  
 This is the original data, corresponding directly to what is in
 methods.cc. Later, supergeneric methods are expanded for all groups
 to produce md_data.
 
 Defined in methods.cc.
 */
extern const Array<MdRecord> md_data_raw;

//! Lookup information for workspace methods.
/*!
 This is the data with expanded supergeneric methods. That means,
 e.g., instead of supergeneric method Copy(Any,Any) there will be
 Copy(Vector,Vector), Copy(Matrix,Matrix), etc..
 
 Defined in methods_aux.cc.
 */
extern const Array<MdRecord> md_data;

//! The map associated with md_data.
/**
 Defined in methods_aux.cc.
 */
extern const map<String, Index> MdMap;

//! The map associated with md_data_raw.
/**
 Defined in methods_aux.cc.
 */
extern const map<String, Index> MdRawMap;

//! The lookup information for the agendas.
/**
 Defined in agendas.cc.
 */
extern const Array<AgRecord> agenda_data;

//! The map associated with agenda_data.
/**
 Defined in agenda_record.cc.
 */
extern const map<String, Index> AgendaMap;

//! The names associated with Wsv groups as Strings.
/**
 See function define_wsv_group_names for more information.
 
 Defined in groups.cc.
 */
extern const ArrayOfString wsv_group_names;

//! The map associated with wsv_group_names.
/**
 Defined in groups.cc.
 */
extern const map<String, Index> WsvGroupMap;

//! Species Data
/**
 Defined in species_data.cc.
 */
extern const Array<SpeciesRecord> species_data;

//! The lookup data for the different lineshapes.
/**
 Defined in lineshapes.cc.
 */
extern const Array<LineshapeRecord> lineshape_data;

//! The lookup data for the different normalization factors to the lineshapes.
/**
 Defined in lineshapes.cc.
 */
extern const Array<LineshapeNormRecord> lineshape_norm_data;

} /* namespace global_data */

#endif /* global_data_h */
