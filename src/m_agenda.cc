/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri Mar 15 09:12:46 2002
  
  \brief  Workspace methods for Agenda.
*/

#include "agenda.h"

void AgendaSet(// WS Generic Output:
	       Agenda& output_agenda,
	       // WS Generic Output Names:
	       const String& agenda_name,
	       // Agenda from controlfile:
	       const Agenda& input_agenda)
{
  output_agenda.resize(input_agenda.nelem());
  output_agenda = input_agenda;
  output_agenda.set_name(agenda_name);
}

void Main(// Agenda from controlfile:
	  const Agenda& input_agenda)
{
  input_agenda.execute();
}
