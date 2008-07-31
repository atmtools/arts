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
                   const Agenda& this_agenda)
{
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
               const Agenda& input_agenda)
{
  // Make external data visible
  extern const Array<AgRecord>  agenda_data;
  extern const map<String, Index> AgendaMap;

  output_agenda = input_agenda;
  output_agenda.set_name(agenda_name);

  // First we have to find the lookup information for this agenda. We
  // use AgendaMap for this.

  map<String, Index>::const_iterator mi =
    AgendaMap.find( output_agenda.name() );

  // Find return end() if the string is not found. This means that the
  // lookup data for this agenda is missing!
  assert( mi != AgendaMap.end() );

  const AgRecord& this_data = agenda_data[mi->second];

  // Ok, we have the lookup data now.

  // Check that the output produced by the actual methods in the
  // agenda corresponds to what is desired in the lookup data:
  for ( Index i=0; i<this_data.Output().nelem(); ++i )
    {
      // The WSV for which to check:
      Index this_wsv = this_data.Output()[i];

      if ( !output_agenda.is_output(this_wsv) )
        {
          ostringstream os;
          os << "The agenda " << output_agenda.name()
             << " must generate the output WSV "
             << Workspace::wsv_data[this_wsv].Name() << ",\n"
             << "but it does not. It only generates:\n";
          for ( Index j=0; j<Workspace::wsv_data.nelem(); ++j )
            if ( output_agenda.is_output(j) )
              os << Workspace::wsv_data[j].Name() << "\n";
          throw runtime_error (os.str());
        }
    }

  // Check that the input used by the actual methods in the
  // agenda corresponds to what is desired in the lookup data:
  for ( Index i=0; i<this_data.Input().nelem(); ++i )
    {
      // The WSV for which to check:
      Index this_wsv = this_data.Input()[i];

      if ( !output_agenda.is_input(ws, this_wsv) )
        {
          ostringstream os;
          os << "The agenda " << output_agenda.name()
             << " must use the input WSV "
             << Workspace::wsv_data[this_wsv].Name() << ",\n"
             << "but it does not. It only uses:\n";
          for ( Index j=0; j<Workspace::wsv_data.nelem(); ++j )
            if ( output_agenda.is_input(ws, j) )
              os << Workspace::wsv_data[j].Name() << "\n";
          throw runtime_error (os.str());
        }
    }

  output_agenda.set_outputs_to_push_and_dup ();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Arts(Workspace& ws,
          // Agenda from controlfile:
          const Agenda&    input_agenda)
{
  input_agenda.execute(ws);
}

