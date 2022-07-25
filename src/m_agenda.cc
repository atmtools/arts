/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \file   m_agenda.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Mar 15 09:12:46 2002
  
  \brief  Workspace methods for Agenda.
*/

#include <algorithm>
#include <map>
#include "agenda_class.h"
#include "agenda_record.h"
#include "debug.h"
#include "global_data.h"
#include "messages.h"
#include "workspace_ng.h"
#include "wsv_aux.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void AgendaExecute(Workspace& ws,
                   // WS Generic Input:
                   const Agenda& this_agenda,
                   const Verbosity& verbosity) {
  CREATE_OUT3;
  out3 << "  Manual agenda execution\n";

  using global_data::agenda_data;
  using global_data::AgendaMap;

  if (!this_agenda.checked()) {
    std::ostringstream os;
    if (this_agenda.name().nelem() == 0)
      os << "This agenda is uninitialized. We don't even know its name.";
    else
      os << "*" << this_agenda.name() << "* is uninitialized.";
    os << " Use *AgendaSet* or *@arts_agenda* to add methods to it.";
    throw std::runtime_error(os.str());
  }

  const AgRecord& agr = agenda_data[AgendaMap.find(this_agenda.name())->second];

  // Duplicate input-only arguments of the agenda as they might be
  // changed inside the agenda.
  const ArrayOfIndex& ain = agr.In();
  const ArrayOfIndex& aout = agr.Out();

  // Put the input and outputs into new sets to sort them. Otherwise
  // set_difference screws up.
  set<Index> sain;
  set<Index> saout;
  sain.insert(ain.begin(), ain.end());
  saout.insert(aout.begin(), aout.end());

  set<Index> in_only;
  set_difference(sain.begin(),
                 sain.end(),
                 saout.begin(),
                 saout.end(),
                 insert_iterator<set<Index> >(in_only, in_only.begin()));
  for (Index it : in_only) {
    ws.duplicate(it);
  }

  const ArrayOfIndex& outputs_to_push = this_agenda.get_output2push();
  const ArrayOfIndex& outputs_to_dup = this_agenda.get_output2dup();

  for (Index it : outputs_to_push) {
    if (ws.is_initialized(it))
      ws.duplicate(it);
    else
      ws.emplace(it);
  }

  for (auto it : outputs_to_dup) ws.duplicate(it);

  String msg;
  bool agenda_failed = false;

  try {
    this_agenda.execute();
  } catch (const std::exception& e) {
    agenda_failed = true;
    msg = e.what();
  }

  for (auto it : outputs_to_push) ws.pop(it);

  for (auto it : outputs_to_dup) ws.pop(it);

  for (auto it : in_only) ws.pop(it);

  ARTS_USER_ERROR_IF(agenda_failed,
                     "Run-time error in agenda: ",
                     this_agenda.name(),
                     '\n',
                     msg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfAgendaExecute(Workspace& ws,
                          // WS Generic Input:
                          const Index& agenda_array_index,
                          const ArrayOfAgenda& agenda_array,
                          const Verbosity& verbosity) {
  if (agenda_array_index < 0 || agenda_array_index >= agenda_array.nelem()) {
    std::ostringstream os;
    os << "Agenda index " << agenda_array_index
       << " out of bounds. 0 <= index < " << agenda_array.nelem();
    throw std::runtime_error(os.str());
  }
  AgendaExecute(ws, agenda_array[agenda_array_index], verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AgendaExecuteExclusive(Workspace& ws,
                            // WS Generic Input:
                            const Agenda& this_agenda,
                            const Verbosity& verbosity) {
  CREATE_OUT3;
  out3 << "  Manual, exclusive agenda execution\n";

#pragma omp critical(AgendaExecuteExclusive_region)
  AgendaExecute(ws, this_agenda, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AgendaSet(Workspace&,
               // WS Generic Output:
               Agenda& output_agenda,
               // WS Generic Output Names:
               const String& agenda_name,
               // Agenda from controlfile:
               const Agenda& input_agenda,
               const Verbosity& verbosity) {
  output_agenda = input_agenda;
  output_agenda.set_name(agenda_name);

  output_agenda.check(verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfAgendaAppend(Workspace&,
                         // WS Generic Output:
                         ArrayOfAgenda& out,
                         // WS Generic Names:
                         const String& agenda_name,
                         // Agenda from controlfile:
                         const Agenda& input_agenda,
                         const Verbosity& verbosity) {
  out.push_back(input_agenda);

  Agenda& appended_agenda = out[out.nelem() - 1];
  appended_agenda.set_name(agenda_name);

  appended_agenda.check(verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AgendaAppend(Workspace&,
                  // WS Generic Output:
                  Agenda& output_agenda,
                  // WS Generic Output Names:
                  const String& output_agenda_name,
                  // WS Generic Input:
                  const Agenda& in_agenda _U_,
                  // WS Generic Input Names:
                  const String& in_agenda_name,
                  // Agenda from controlfile:
                  const Agenda& input_agenda,
                  const Verbosity& verbosity) {
  if (output_agenda_name != in_agenda_name) {
    ostringstream os;
    os << "Output and input agenda must be the same!" << endl
       << "*" << output_agenda_name << "* and *" << in_agenda_name << "* "
       << "are not.";
    throw runtime_error(os.str());
  }

  Array<MRecord> methods = output_agenda.Methods();
  for (Index i = 0; i < input_agenda.Methods().nelem(); i++)
    methods.push_back(input_agenda.Methods()[i]);

  output_agenda.set_methods(methods);
  output_agenda.check(verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Arts(Workspace&,
          // Agenda from controlfile:
          const Agenda&,
          const Verbosity& verbosity) {
  CREATE_OUT0;
  arts_exit_with_error_message(
      "The 'Arts' method is obsolete. Arts1 controlfiles are no longer supported.",
      out0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Arts2(Workspace& ws,
           // Agenda from controlfile:
           const Agenda& input_agenda,
           const Verbosity& verbosity) {
  Verbosity* v = static_cast<Verbosity*>(ws[ws.WsvMap.find("verbosity") -> second].get());

  // If the verbosity in the current workspace and the verbosity parameter point
  // to the same variable in memory, that means we were called
  // from inside a controlfile by the user. That is not permitted.
  if (v == &verbosity) {
    CREATE_OUT0;
    arts_exit_with_error_message("The 'Arts2' method can't be called directly.",
                                 out0);
  }
  (*v) = verbosity;
  input_agenda.execute();
}
