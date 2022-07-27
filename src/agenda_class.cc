/* Copyright (C) 2002-2012
   Stefan Buehler <sbuehler@ltu.se>
   Oliver Lemke <olemke@core-dump.info>

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
  \file   agenda_class.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Thu Mar 14 08:49:33 2002
  
  \brief  Implementation of agendas.
*/

#include <algorithm>
#include <iterator>
#include <memory>
#include <ostream>

#include "agenda_class.h"
#include "agenda_record.h"  // only for test functions
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "global_data.h"
#include "messages.h"
#include "methods.h"
#include "workspace_ng.h"

MRecord::MRecord(const std::shared_ptr<Workspace>& ws)
    : moutput(), minput(), msetvalue(), mtasks(ws) { /* Nothing to do here. */
}

MRecord::MRecord(const Index id,
                 ArrayOfIndex output,
                 ArrayOfIndex input,
                 const TokVal& setvalue,
                 const Agenda& tasks,
                 bool internal)
    : mid(id),
      moutput(std::move(output)),
      minput(std::move(input)),
      msetvalue(setvalue),
      mtasks(tasks),
      minternal(internal) {
}

Agenda::Agenda(const std::shared_ptr<Workspace>& workspace)
    : ws(workspace), mname(), mml(), moutput_push(), moutput_dup() {
}

//! Appends methods to an agenda
/*!
  This function appends a workspace method to the agenda. It currently only
  supports appending WSMs which have no generic input or output, and which are
  defined as Set methods which take only one keyword.
   
  The keyword value has to be a string, which for no value should be of length
  zero.
   
  \param methodname    The name of the WSM
  \param keywordvalue  The value of the keyword

  \author Mattias Ekstrom
  \date   2005-01-05
*/
void Agenda::append(const String& methodname, const TokVal& keywordvalue) {
  using global_data::MdMap;

  const auto i2 = MdMap.find(methodname);
  ARTS_ASSERT(i2 != MdMap.end());
  Index id = i2->second;

  using global_data::md_data;
  ArrayOfIndex output = md_data[id].Out();

  // Find explicit method id in MdMap.
  ArrayOfIndex input = md_data[id].InOnly();

  mml.push_back(MRecord(id, output, input, keywordvalue, Agenda(ws)));
  mchecked = false;
}

//! Checks consistency of an agenda.
/*! 
  Checks that the input used by the agenda and the output produced by the
  actual methods corresponds to what is desired in the lookup data.
*/
void Agenda::check(Workspace& ws_in, const Verbosity& verbosity) {
  // Check that this agenda has a default workspace
  if (not ws.get()) {
    mchecked = false;
    return;
  }

  // Make external data visible
  using global_data::agenda_data;
  using global_data::AgendaMap;

  // First we have to find the lookup information for this agenda. We
  // use AgendaMap for this.

  auto mi = AgendaMap.find(mname);

  // Find return end() if the string is not found. This means we deal with
  // agenda defined in the control file and therefore we don't check its
  // consistency. Custom agendas can't be executed and we delay the check
  // until it is copied to a predefined agenda.
  if (mi == AgendaMap.end()) {
    mchecked = false;
    return;
  }

  const AgRecord& this_data = agenda_data[mi->second];

  // Ok, we have the lookup data now.

  // Check that the output produced by the actual methods in the
  // agenda corresponds to what is desired in the lookup data:
  for (Index i = 0; i < this_data.Out().nelem(); ++i) {
    // The WSV for which to check: 
    Index this_wsv = this_data.Out()[i];

    if (!is_output(this_wsv)) {
      ostringstream os;
      os << "The agenda " << mname << " must generate the output WSV "
         << (*ws_in.wsv_data_ptr)[this_wsv].Name() << ",\n"
         << "but it does not. It only generates:\n";
      for (Index j = 0; j < (*ws_in.wsv_data_ptr).nelem(); ++j)
        if (is_output(j)) os << (*ws_in.wsv_data_ptr)[j].Name() << "\n";
      throw runtime_error(os.str());
    }
  }

  // Check that the input used by the actual methods in the
  // agenda corresponds to what is desired in the lookup data:
  for (Index i = 0; i < this_data.In().nelem(); ++i) {
    // The WSV for which to check:
    Index this_wsv = this_data.In()[i];

    if (!is_input(ws_in, this_wsv)) {
      ostringstream os;
      os << "The agenda " << mname << " must use the input WSV "
         << (*ws_in.wsv_data_ptr)[this_wsv].Name() << ",\n"
         << "but it does not. It only uses:\n";
      for (Index j = 0; j < (*ws_in.wsv_data_ptr).nelem(); ++j)
        if (is_input(ws_in, j)) os << (*ws_in.wsv_data_ptr)[j].Name() << "\n";
      throw runtime_error(os.str());
    }
  }

  set_outputs_to_push_and_dup(verbosity);

  mchecked = true;
}

//! Execute an agenda.
/*! 
  This executes the methods specified in tasklist on the given
  workspace. It also checks for errors during the method execution and
  stops the program if an error has occured. 
*/
void Agenda::execute(Workspace& ws_in) const {
  ARTS_USER_ERROR_IF(not correct_workspace(ws_in),
                     mname,
                     " is on another original workspace.  DEBUG: Pointers are at: ",
                     ws_in.original_workspace,
                     " vs ",
                     ws -> original_workspace)

  if (!mchecked) {
    ostringstream os;
    os << "Agenda *" << mname << "* hasn't been checked for consistency yet."
       << endl
       << "This check is usually done by AgendaSet or AgendaAppend." << endl
       << "There are two possible causes for this:" << endl
       << "1) You're trying to execute an agenda that has been created in"
       << endl
       << "   the controlfile with AgendaCreate. This is not allowed. You have"
       << endl
       << "   to use *Copy* to store it into one of the predefined agendas and"
       << endl
       << "   execute that one." << endl
       << "2) Developer error: If you have written code that modifies an Agenda"
       << endl
       << "   directly (changing its name or altering its method list), it's up"
       << endl
       << "   to you to call Agenda::check in your code after your modifications.";
    throw runtime_error(os.str());
  }

  // An empty Agenda name indicates that something going wrong here
  ARTS_ASSERT(mname != "");

  // The method description lookup table:
  using global_data::md_data;

  // The array holding the pointers to the getaway functions:
  extern void (*getaways[])(Workspace&, const MRecord&);

  const Index wsv_id_verbosity = ws_in.WsvMap_ptr->at("verbosity");
  ws_in.duplicate(wsv_id_verbosity);

  Verbosity& averbosity =
      *(static_cast<Verbosity*>(ws_in[wsv_id_verbosity].get()));

  averbosity.set_main_agenda(is_main_agenda());

  ArtsOut1 aout1(averbosity);
  {
    //    ostringstream os;  // disabled for performance reasons
    //    os << "Executing " << name() << "\n"
    //       << "{\n";
    //    aout1 << os.str();
    aout1 << "Executing " << name() << "\n"
          << "{\n";
  }

  for (Index i = 0; i < mml.nelem(); ++i) {
    const Verbosity& verbosity =
        *static_cast<Verbosity*>(ws_in[wsv_id_verbosity].get());
    CREATE_OUT1;
    CREATE_OUT3;

    // Runtime method data for this method:
    const MRecord& mrr = mml[i];
    // Method data for this method:
    const MdRecord& mdd = md_data[mrr.Id()];

    try {
      {
        if (mrr.isInternal()) {
          out3 << "- " + mdd.Name() + "\n";
        } else {
          out1 << "- " + mdd.Name() + "\n";
        }
      }

      {  // Check if all input variables are initialized:
        const ArrayOfIndex& v(mrr.In());
        for (Index s = 0; s < v.nelem(); ++s) {
          if ((s != v.nelem() - 1 || !mdd.SetMethod()) &&
              !ws_in.is_initialized(v[s]))
            throw runtime_error(
                "Method " + mdd.Name() +
                " needs input variable: " + (*ws_in.wsv_data_ptr)[v[s]].Name());
        }
      }

      {  // Check if all output variables which are also used as input
        // are initialized
        const ArrayOfIndex& v = mdd.InOut();
        for (Index s = 0; s < v.nelem(); ++s)
          if (!ws_in.is_initialized(mrr.Out()[v[s]]))
            throw runtime_error("Method " + mdd.Name() +
                                " needs input variable: " +
                                (*ws_in.wsv_data_ptr)[mrr.Out()[v[s]]].Name());
      }

      // Call the getaway function:
      getaways[mrr.Id()](ws_in, mrr);

    } catch (const std::bad_alloc& x) {
      aout1 << "}\n";

      ostringstream os;
      os << "Memory allocation error in method: " << mdd.Name() << '\n'
         << "For memory intensive jobs it could help to limit the\n"
         << "number of threads with the -n option.\n"
         << x.what();

      throw runtime_error(os.str());
    } catch (const std::exception& x) {
      aout1 << "}\n";

      ostringstream os;
      os << "Run-time error in method: " << mdd.Name() << '\n' << x.what();

      throw runtime_error(os.str());
    }
  }

  aout1 << "}\n";

  ws_in.pop(wsv_id_verbosity);
}

//! Retrieve indexes of all input and output WSVs
/*!
  Builds arrays of WSM output variables which need to be
  duplicated or pushed on the WSV stack before the agenda
  is executed.
*/
void Agenda::set_outputs_to_push_and_dup(const Verbosity& verbosity) {
  using global_data::agenda_data;
  using global_data::AgendaMap;
  using global_data::md_data;
  using global_data::MdMap;
  using global_data::WsvGroupMap;

  set<Index> inputs;
  set<Index> outputs;
  set<Index> outs2push;
  set<Index> outs2dup;

  const Index WsmAgendaExecuteIndex = MdMap.find("AgendaExecute")->second;
  const Index WsmAgendaExecuteExclIndex =
      MdMap.find("AgendaExecuteExclusive")->second;
  const Index WsmDeleteIndex = MdMap.find("Delete")->second;
  const Index WsvAgendaGroupIndex = WsvGroupMap.find("Agenda")->second;

  for (auto & method : mml) {
    // Collect output WSVs
    const ArrayOfIndex& gouts = method.Out();

    // Put the outputs into a new set to sort them. Otherwise
    // set_intersection and set_difference screw up.
    set<Index> souts;
    souts.insert(gouts.begin(), gouts.end());

    // Collect generic input WSVs
    const ArrayOfIndex& gins = method.In();
    inputs.insert(gins.begin(), gins.end());

    /* Special case: For the Delete WSM add its input to the list
       * of output variables to force a duplication of those variables.
       * It avoids deleting variables outside the agenda's scope.
       */
    if (method.Id() == WsmDeleteIndex) {
      souts.insert(gins.begin(), gins.end());
    }

    /* Special case: For WSMs that execute generic input Agendas
         it is necessary for proper scoping to also add the output and input variables of the
         agenda to our output and input variable lists.
       */
    if (method.Id() == WsmAgendaExecuteIndex ||
        method.Id() == WsmAgendaExecuteExclIndex) {
      for (Index j = 0; j < md_data[method.Id()].GInType().nelem(); j++) {
        if (md_data[method.Id()].GInType()[j] == WsvAgendaGroupIndex) {
          const String& agenda_name = (*ws->wsv_data_ptr)[gins[j]].Name();
          const auto agenda_it =
              AgendaMap.find(agenda_name);
          // The executed agenda must not be a user created agenda
          if (agenda_it == AgendaMap.end()) {
            ostringstream os;
            os << "Manual execution of the agenda \"" << agenda_name
               << "\" is not supported.";
            throw std::runtime_error(os.str());
          }
          const ArrayOfIndex& agouts = agenda_data[agenda_it->second].Out();
          souts.insert(agouts.begin(), agouts.end());
          const ArrayOfIndex& agins = agenda_data[agenda_it->second].In();
          inputs.insert(agins.begin(), agins.end());
        }
      }
    }

    // Add all outputs of this WSM to global list of outputs
    outputs.insert(souts.begin(), souts.end());

    // Find out all output WSVs of current WSM which were
    // already used as input. We have to place a copy of them on
    // the WSV stack.
    set_intersection(souts.begin(),
                     souts.end(),
                     inputs.begin(),
                     inputs.end(),
                     insert_iterator<set<Index> >(outs2dup, outs2dup.begin()));

    // Get a handle on the InOut list for the current method:
    const ArrayOfIndex& mdinout = md_data[method.Id()].InOut();
    const ArrayOfIndex& output = method.Out();

    for (Index j = 0; j < mdinout.nelem(); j++) {
      inputs.insert(output[mdinout[j]]);
    }
  }

  // Find all outputs which are not in the list of WSVs to duplicate
  set_difference(outputs.begin(),
                 outputs.end(),
                 outs2dup.begin(),
                 outs2dup.end(),
                 insert_iterator<set<Index> >(outs2push, outs2push.begin()));

  const AgRecord& agr = agenda_data[AgendaMap.find(name())->second];
  const ArrayOfIndex& aout = agr.Out();
  const ArrayOfIndex& ain = agr.In();

  // We have to build a new set of agenda input and output because the
  // set_difference function only works properly on sorted input.
  set<Index> saout;
  set<Index> sain;

  saout.insert(aout.begin(), aout.end());
  sain.insert(ain.begin(), ain.end());

  moutput_push.clear();
  moutput_dup.clear();

  // Remove the WSVs which are agenda input from the list of
  // output variables for which we have to create an new
  // entry on the stack. This is already done for agenda inputs.
  set<Index> outs2push_without_agenda_input;
  set_difference(
      outs2push.begin(),
      outs2push.end(),
      sain.begin(),
      sain.end(),
      insert_iterator<set<Index> >(outs2push_without_agenda_input,
                                   outs2push_without_agenda_input.begin()));

  // Same for agenda output variables.
  set_difference(
      outs2push_without_agenda_input.begin(),
      outs2push_without_agenda_input.end(),
      saout.begin(),
      saout.end(),
      insert_iterator<ArrayOfIndex>(moutput_push, moutput_push.begin()));

  // Remove the WSVs which are agenda input from the list of
  // output variables for which we have to create a duplicate
  // on the stack. This is already done for agenda inputs.
  set<Index> outs2dup_without_agenda_input;
  set_difference(
      outs2dup.begin(),
      outs2dup.end(),
      sain.begin(),
      sain.end(),
      insert_iterator<set<Index> >(outs2dup_without_agenda_input,
                                   outs2dup_without_agenda_input.begin()));

  // Same for agenda output variables.
  set_difference(
      outs2dup_without_agenda_input.begin(),
      outs2dup_without_agenda_input.end(),
      saout.begin(),
      saout.end(),
      insert_iterator<ArrayOfIndex>(moutput_dup, moutput_dup.begin()));

  // Special case: Variables which are defined in the agenda only
  // as output but are used first as input in one of the WSMs
  // For those the current WSV value must be copied to the agenda
  // input variable
  set<Index> saout_only;

  set_difference(saout.begin(),
                 saout.end(),
                 sain.begin(),
                 sain.end(),
                 insert_iterator<set<Index> >(saout_only, saout_only.begin()));

  ArrayOfIndex agenda_only_out_wsm_in;
  set_intersection(outs2dup.begin(),
                   outs2dup.end(),
                   saout_only.begin(),
                   saout_only.end(),
                   insert_iterator<ArrayOfIndex>(
                       agenda_only_out_wsm_in, agenda_only_out_wsm_in.begin()));

  // Special case: Variables which are defined in the agenda only
  // as input but are used as output in one of the WSMs
  // For those the current WSV value must be copied to the agenda
  // input variable
  set<Index> sain_only;

  set_difference(sain.begin(),
                 sain.end(),
                 saout.begin(),
                 saout.end(),
                 insert_iterator<set<Index> >(sain_only, sain_only.begin()));

  ArrayOfIndex agenda_only_in_wsm_out;
  set_intersection(outs2push.begin(),
                   outs2push.end(),
                   sain_only.begin(),
                   sain_only.end(),
                   insert_iterator<ArrayOfIndex>(
                       agenda_only_in_wsm_out, agenda_only_in_wsm_out.begin()));

  CREATE_OUT3;

  out3 << "  [Agenda::pushpop]                 : " << name() << "\n";
  out3 << "  [Agenda::pushpop] - # Funcs in Ag : " << mml.nelem() << "\n";
  out3 << "  [Agenda::pushpop] - AgOut         : ";
  PrintWsvNames(out3, *ws, aout);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - AgIn          : ";
  PrintWsvNames(out3, *ws, ain);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - All WSM output: ";
  PrintWsvNames(out3, *ws, outputs);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - All WSM input : ";
  PrintWsvNames(out3, *ws, inputs);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - Output WSVs push     : ";
  PrintWsvNames(out3, *ws, moutput_push);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - Output WSVs dup      : ";
  PrintWsvNames(out3, *ws, moutput_dup);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - Ag inp dup    : ";
  PrintWsvNames(out3, *ws, agenda_only_in_wsm_out);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - Ag out dup    : ";
  PrintWsvNames(out3, *ws, agenda_only_out_wsm_in);
  out3 << "\n";
}

//! Check if given variable is agenda input.
/*! 
  A variable is agenda input if it is an input variable to any of the
  methods making up the agenda. 

  \param[in] var The workspace variable to check.

  \return True if var is an input variable of this agenda.
*/
bool Agenda::is_input(Workspace& ws_in, Index var) const {
  // Make global method data visible:
  using global_data::agenda_data;
  using global_data::AgendaMap;
  using global_data::md_data;
  using global_data::MdMap;
  using global_data::WsvGroupMap;

  const Index WsmAgendaExecuteIndex = MdMap.find("AgendaExecute")->second;
  const Index WsmAgendaExecuteExclIndex =
      MdMap.find("AgendaExecuteExclusive")->second;

  // Make sure that var is the index of a valid WSV:
  ARTS_ASSERT(0 <= var);
  ARTS_ASSERT(var < (*ws_in.wsv_data_ptr).nelem());

  // Determine the index of WsvGroup Agenda
  const Index WsvAgendaGroupIndex = WsvGroupMap.find("Agenda")->second;

  // Loop all methods in this agenda:
  for (Index i = 0; i < nelem(); ++i) {
    // Get a handle on this methods runtime data record:
    const MRecord& this_method = mml[i];

    // Is var a generic input?
    {
      // Get a handle on the Input list:
      const ArrayOfIndex& input = this_method.In();

      for (Index j = 0; j < input.nelem(); ++j) {
        if (var == input[j]) return true;
      }

      // Get a handle on the InOut list for the current method:
      const ArrayOfIndex& mdinout = md_data[this_method.Id()].InOut();
      const ArrayOfIndex& output = this_method.Out();

      for (Index j = 0; j < mdinout.nelem(); j++) {
        if (var == output[mdinout[j]]) return true;
      }

      /* Special case: For WSMs that execute generic input Agendas
         it is necessary to check their inputs.
       */
      if (this_method.Id() == WsmAgendaExecuteIndex ||
          this_method.Id() == WsmAgendaExecuteExclIndex) {
        for (Index j = 0; j < md_data[this_method.Id()].GInType().nelem();
             j++) {
          if (md_data[this_method.Id()].GInType()[j] == WsvAgendaGroupIndex) {
            const String& agenda_name = (*ws_in.wsv_data_ptr)[input[j]].Name();
            const auto agenda_it =
                AgendaMap.find(agenda_name);
            // The executed agenda must not be a user created agenda
            if (agenda_it == AgendaMap.end()) {
              ostringstream os;
              os << "Manual execution of the agenda \"" << agenda_name
                 << "\" is not supported.";
              throw std::runtime_error(os.str());
            }
            const ArrayOfIndex& agins = agenda_data[agenda_it->second].In();
            for (Index k = 0; k < agins.nelem(); ++k) {
              if (var == agins[k]) return true;
            }
          }
        }
      }
    }
  }

  // Ok, that means var is no input at all.
  return false;
}

//! Check if given variable is agenda output.
/*! 
  A variable is agenda output if it is an output variable to any of the
  methods making up the agenda. 

  \param var The workspace variable to check.

  \return True if var is an output variable of this agenda.
*/
bool Agenda::is_output(Index var) const {
  // Make global method data visible:
  using global_data::agenda_data;
  using global_data::AgendaMap;
  using global_data::md_data;
  using global_data::MdMap;
  using global_data::WsvGroupMap;

  const Index WsmAgendaExecuteIndex = MdMap.find("AgendaExecute")->second;
  const Index WsmAgendaExecuteExclIndex =
      MdMap.find("AgendaExecuteExclusive")->second;

  // Determine the index of WsvGroup Agenda
  const Index WsvAgendaGroupIndex = WsvGroupMap.find("Agenda")->second;

  // Loop all methods in this agenda:
  for (Index i = 0; i < nelem(); ++i) {
    // Get a handle on this methods runtime data record:
    const MRecord& this_method = mml[i];

    // Is var a generic output?
    {
      // Get a handle on the Output list:
      const ArrayOfIndex& output = this_method.Out();

      for (Index j = 0; j < output.nelem(); ++j) {
        if (var == output[j]) return true;
      }

      /* Special case: For WSMs that execute generic input Agendas
         it is necessary to check their outputs.
        */
      if (this_method.Id() == WsmAgendaExecuteIndex ||
          this_method.Id() == WsmAgendaExecuteExclIndex) {
        for (Index j = 0; j < md_data[this_method.Id()].GInType().nelem();
             j++) {
          if (md_data[this_method.Id()].GInType()[j] == WsvAgendaGroupIndex) {
            const String& agenda_name =
                (*ws->wsv_data_ptr)[this_method.In()[j]].Name();
            const auto agenda_it =
                AgendaMap.find(agenda_name);
            // The executed agenda must not be a user created agenda
            if (agenda_it == AgendaMap.end()) {
              ostringstream os;
              os << "Manual execution of the agenda \"" << agenda_name
                 << "\" is not supported.";
              throw std::runtime_error(os.str());
            }
            const ArrayOfIndex& agouts = agenda_data[agenda_it->second].Out();
            for (Index k = 0; k < agouts.nelem(); k++)
              if (agouts[k] == var) return true;
          }
        }
      }
    }
  }

  // Ok, that means var is no output at all.
  return false;
}

//! Set agenda name.
/*! 
  This sets the private member mname to the given string. 

  \param nname The name for the agenda.
*/
void Agenda::set_name(const String& nname) {
  mname = nname;
  mchecked = false;
}

//! Agenda name.
/*! 
  Returns the private member mname.

  \return The name of this agenda.
*/
String Agenda::name() const { return mname; }

//! Check if method is in Agenda.
/*!
  This function checks if the method with the given name is
  called by this agenda.

  \param methodname     Name of method to look for.

  \return True if method is part of Agenda.

  \author Oliver Lemke
  \date   2013-03-18
*/
bool Agenda::has_method(const String& methodname) const {
  using global_data::md_data;

  bool found = false;
  for (auto it = mml.begin();
       !found && it != mml.end();
       it++) {
    if (md_data[it->Id()].Name() == methodname) found = true;
  }

  return found;
}

void Agenda::set_methods(const Array<MRecord>& ml) {
  mml = ml;
  mchecked = false;
}

//! Print an agenda.
/*!
  This prints an agenda, by printing the individual methods, just as
  they would appear in the controlfile.

  \param os     Output stream.
  \param indent How many characters of indentation.

  \author Stefan Buehler
  \date   2002-12-02
*/
void Agenda::print(ostream& os, const String& indent) const {
  for (Index i = 0; i < mml.nelem(); ++i) {
    // Print member methods with 3 characters more indentation:
    mml[i].print(os, indent);
  }
}

//! Resize the method list.
/*!
  Resizes the agenda's method list to n elements
 */
void Agenda::resize(Index n) { mml.resize(n); }

//! Return the number of agenda elements.
/*!  
  This is needed, so that we can find out the correct size for
  resize, befor we do a copy.

  \return Number of agenda elements.
*/
Index Agenda::nelem() const { return mml.nelem(); }

//! Append a new method to end of list.
/*! 
  This is used by the parser to fill up the agenda.

  \param n New method to add.
*/
void Agenda::push_back(const MRecord& n) {
  mml.push_back(n);
  mchecked = false;
}

Agenda& Agenda::operator=(const Agenda& x) {
  ws = x.ws;
  mml = x.mml;
  mname = x.mname;
  moutput_push = x.moutput_push;
  moutput_dup = x.moutput_dup;
  mchecked = x.mchecked;
  return *this;
}

//! Output operator for Agenda.
/*! 
  This is useful for debugging.
  
  \param os Output stream.
  \param a The Agenda to write.

  \return Output stream.

  \author Stefan Buehler
  \date   2002-12-02
*/
ostream& operator<<(ostream& os, const Agenda& a) {
  // Print agenda as it would apear in a controlfile.
  a.print(os, "    ");
  return os;
}

//--------------------------------
//     Functions for MRecord:
//--------------------------------

//! Print an MRecord.
/*!
  Since the MRecord contains all runtime information for one method,
  the best way to print it is exactly as it would appear in the
  controlfile. 

  This has to work in a recursive way, since the method can be an
  agenda method, which includes other methods, which can be agenda
  methods, ...

  Therefore, the indentation is increased more and more for recursive
  calls. 

  At the moment, this is used just for debugging.

  \param os     Output stream.
  \param indent How many characters of indentation.

  \author Stefan Buehler
  \date   2002-12-02
*/
void MRecord::print(ostream& os, const String& indent) const {
  using global_data::md_data;

  // Get a handle on the right record:
  const MdRecord& tmd = md_data[Id()];

  os << indent << tmd.Name();

  // Is this a generic method? -- Then we need round braces.
  if (0 != tmd.GOutType().nelem() + tmd.GInType().nelem()) {
    // First entry needs no leading comma:
    bool first = true;

    os << "(";

    for (Index i = 0; i < Out().nelem(); ++i) {
      if (first)
        first = false;
      else
        os << ",";

      os << (*mtasks.wsptr() -> wsv_data_ptr)[Out()[i]].Name();
    }

    for (Index i = 0; i < In().nelem(); ++i) {
      if (first)
        first = false;
      else
        os << ",";

      os << (*mtasks.wsptr() -> wsv_data_ptr)[In()[i]].Name();
    }

    os << ")";
  }

  if (0 != Tasks().nelem()) {
    os << " {\n";
    Tasks().print(os, indent + "    ");
    os << indent << "}\n";
  } else
    os << "\n";
}

MRecord& MRecord::operator=(const MRecord& x) {
  mid = x.mid;

  msetvalue = x.msetvalue;

  moutput = x.moutput;

  minput = x.minput;

  mtasks = x.mtasks;

  return *this;
}

void MRecord::ginput_only(ArrayOfIndex& ginonly) const {
  ginonly = minput;  // Input
  for (auto j = moutput.begin(); j < moutput.end(); ++j)
    for (auto k = ginonly.begin(); k < ginonly.end(); ++k)
      if (*j == *k) {
        //              erase_vector_element(vi,k);
        k = ginonly.erase(k) - 1;
        // We need the -1 here, otherwise due to the
        // following increment we would miss the element
        // behind the erased one, which is now at the
        // position of the erased one.
      }
}

//! Output operator for MRecord.
/*! 
  This is useful for debugging.
  
  \param os Output stream.
  \param a The method runtime data record to write.

  \return Output stream.

  \author Stefan Buehler
  \date   2002-12-02
*/
ostream& operator<<(ostream& os, const MRecord& a) {
  a.print(os, "");
  return os;
}
