#pragma once

#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "workspace_agendas.h"
#include "workspace_methods.h"

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

struct SizeCheck {
  //! Fails the check when false
  std::string test;

  //! What the check means, for the error message and the documentation
  std::string constraint;

  //! What to report when the check fails, as label and the expression to read
  std::vector<std::pair<std::string, std::string>> printables;
};

/*! How a generator writes the name of a workspace variable where it generates code.
 *
 * The workspace variables are not always in scope under their own name.  Agendas
 * and the python bindings have them as named variables, but a generated method
 * body reaches them through the workspace, so it writes something else.
 */
using DimAccess = std::function<std::string(const std::string&)>;

//! Writes the name as it stands, for generators that have the variables in scope
std::string dim_access_by_name(const std::string& name);

/*! Checks that the targets agree with the references about the dimensions they share.
 *
 * The references are tried in order, so put the variables that are trusted to
 * carry the right size first.  A target is never checked against itself.
 */
std::vector<SizeCheck> size_checks(const std::vector<std::vector<std::string>>& references,
                                   const std::vector<std::string>&              targets,
                                   const DimAccess&                             access = dim_access_by_name);

//! The checks of an agenda, verifying its outputs once it has run
std::vector<SizeCheck> agenda_output_size_checks(const WorkspaceAgendaInternalRecord& ag);

//! The checks of a method, verifying what it is given before it runs
std::vector<SizeCheck> method_input_size_checks(const WorkspaceMethodInternalRecord& wsmr,
                                                const DimAccess&                     access = dim_access_by_name);

/*! The checks of a method, verifying what it created once it has run.
 *
 * Only the outputs that the method did not also take as input are checked here,
 * since the others have already been checked before the method ran.
 */
std::vector<SizeCheck> method_output_size_checks(const WorkspaceMethodInternalRecord& wsmr,
                                                 const DimAccess&                     access = dim_access_by_name);

/*! Throws unless every workspace variable names its dimensions sensibly.
 *
 * Checked for all of them rather than only for those that some method or agenda
 * happens to verify, so that naming a dimension wrongly is always caught.
 */
void check_workspace_dimensions();

/*! Writes the C++ that performs a check and reports what failed.
 *
 * This is how every generator turns a check into code, so that a size failure
 * reads the same whether it came from an agenda, a method, or python.
 */
std::string size_check_code(const SizeCheck& check, std::string_view indent);

//! Writes the checks as a documentation rubric, or nothing if there are none
std::string size_check_docs(const std::vector<SizeCheck>& checks);

//! Documents the shape of a workspace variable, or nothing if it names no dimensions
std::string variable_dimension_docs(const std::string& name);
