/*!
  \file   globals_data.h
  \brief  Global variable declarations

  \author Oliver Lemke
  \date 2013-04-25 */

#ifndef global_data_h
#define global_data_h

#include "agenda_record.h"
#include "array.h"
#include "groups.h"
#include "methods.h"
#include "workspace_memory_handler.h"
#include <map>

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
 See function define_wsv_groups for more information.
 
 Defined in groups.cc.
 */
extern const ArrayOfGroupRecord wsv_groups;

//! The map associated with wsv_groups.
/**
 Defined in groups.cc.
 */
extern const map<String, Index> WsvGroupMap;

/** The workspace memory handler
 * Defined in workspace_ng.cc.
 */
extern WorkspaceMemoryHandler workspace_memory_handler;
} /* namespace global_data */

#endif /* global_data_h */
