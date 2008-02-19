/* Copyright (C) 2000-2007 Stefan Buehler <sbuehler@ltu.se>

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

/* Header file for functions and classes related to parsing the
   control text. */ 

#ifndef parser_old_h
#define parser_old_h

#include <map>
#include "parser.h"
#include "agenda_class.h"

void parse_main(Agenda& tasklist, SourceText& text);

void parse_agenda(Agenda& tasklist,
                  SourceText& text);


#endif
