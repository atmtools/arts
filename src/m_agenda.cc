/* Copyright (C) 2002-2008 Stefan Buehler <sbuehler@ltu.se>

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

#include <map>
#include "agenda_class.h"
#include "messages.h"
#include "wsv_aux.h"
#include "agenda_record.h"
#include "workspace_ng.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void AgendaExecute(Workspace& ws,
                   // WS Generic Input:
                   const Agenda& this_agenda,
                   const Verbosity& verbosity)
{
  CREATE_OUT3
  out3 << "  Manual agenda execution\n";
  this_agenda.execute(ws);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AgendaSet(Workspace& ws,
               // WS Generic Output:
               Agenda& output_agenda,
               // WS Generic Output Names:
               const String& agenda_name,
               // Agenda from controlfile:
               const Agenda& input_agenda,
               const Verbosity& verbosity)
{
  output_agenda = input_agenda;
  output_agenda.set_name(agenda_name);

  output_agenda.check(ws, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AgendaAppend(Workspace& ws,
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
                  const Verbosity& verbosity)
{
  if (output_agenda_name != in_agenda_name)
    {
      ostringstream os;
      os << "Output and input agenda must be the same!" << endl
        << "*" << output_agenda_name << "* and *" << in_agenda_name << "* "
        << "are not.";
      throw runtime_error (os.str());
    }

  Array<MRecord> methods = output_agenda.Methods();
  for (Index i = 0; i < input_agenda.Methods().nelem(); i++)
    methods.push_back(input_agenda.Methods()[i]);

  output_agenda.set_methods (methods);
  output_agenda.check(ws, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Arts(Workspace&,
           // Agenda from controlfile:
           const Agenda&,
           const Verbosity&)
{
  arts_exit_with_error_message("The 'Arts' method is obsolete. Arts1 controlfiles are no longer supported.");
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Arts2(Workspace& ws,
           // Agenda from controlfile:
           const Agenda& input_agenda,
           const Verbosity& verbosity)
{
  Verbosity *v = (Verbosity*)ws[get_wsv_id("verbosity")];

  // If the verbosity in the current workspace and the verbosity parameter point
  // to the same variable in memory, that means we were called
  // from inside a controlfile by the user. That is not permitted.
  if (v == &verbosity)
  {
    arts_exit_with_error_message("The 'Arts2' method can't be called directly.");
  }
  (*v) = verbosity;
  input_agenda.execute(ws);
}

