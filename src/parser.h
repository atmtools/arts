/* Copyright (C) 2008 Oliver Lemke <olemke@core-dump.info>

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

#ifndef parser_h
#define parser_h

#include <map>
#include "agenda_class.h"
#include "sourcetext.h"


class ArtsParser {
public:
  ArtsParser(Agenda& tasklist, String controlfile);

  void parse_tasklist ();

private:
  void parse_main();

  void parse_agenda(Agenda& tasklist);

  void parse_method(Index& id, 
                    Array<TokVal>& values,
                    ArrayOfIndex& output,
                    ArrayOfIndex& input,
                    Agenda&       tasks,
                    String&       include_file,
                    bool no_eot=false);

  bool is_whitespace(const char c);

  void eat_whitespace();

  void eat_whitespace_from_string(String& str, size_t& pos);

  void read_name(String& name);

  void assertain_character(char c);

  void parse_String(String& res);

  void read_integer(String& res);

  void read_numeric(String& res);

  void parse_integer(Index& n);

  void parse_numeric(Numeric& n);

  void parse_Stringvector(ArrayOfString& res);

  void parse_intvector(ArrayOfIndex& res);

  void parse_numvector(Vector& res);

  bool parse_numvector_from_string (Vector& res, String& str);

  bool parse_stringarray_from_string (ArrayOfString& res, String& str);

  Agenda& mtasklist;

  String mcfile;

  SourceText msource;
};

#endif /* parser_h */

