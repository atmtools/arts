#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "workspace_agendas.h"

struct WorkspaceDimensionRecord {
  std::string desc;
};

/*! The dimensions that workspace variables may name.
 *
 * A dimension is a size that several workspace variables share.  The variables
 * name their dimensions in WorkspaceVariableInternalRecord::dims and the groups
 * say how to read them in WorkspaceGroupRecord::dim_size.  Nothing here knows
 * how to read a size, because that depends on the type of the variable and not
 * on the dimension itself.
 */
const std::unordered_map<std::string, WorkspaceDimensionRecord>& internal_workspace_dimensions();

/*! The size checks of an agenda, as generated from the dimensions of its variables.
 *
 * Every output variable that names a dimension that some other variable of the
 * agenda also names is checked against that other variable.  Inputs are
 * preferred as the thing to check against, since the agenda is not supposed to
 * change the size of what it was given.
 *
 * Returns nothing if the agenda has turned output size verification off.
 */
std::vector<StringVectorAgendaHelper> agenda_output_size_checks(const WorkspaceAgendaInternalRecord& ag);

/*! Documents the shape of a workspace variable, or nothing if it names no dimensions. */
std::string variable_dimension_docs(const std::string& name);
