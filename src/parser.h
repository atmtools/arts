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
#include "methods.h"


class ArtsParser {
public:
  ArtsParser(Agenda& tasklist, String controlfile);

  void parse_tasklist ();

private:
  void parse_main();

  void parse_agenda(Agenda& tasklist);

  void parse_method(Index& id, 
                    Array<TokVal>& values,
                    ArrayOfIndex&  output,
                    ArrayOfIndex&  input,
                    Agenda&        tasks,
                    ArrayOfIndex&  auto_vars,
                    Array<TokVal>& auto_vars_values,
                    String&        include_file,
                    bool no_eot=false);

  void parse_input(const MdRecord*      mdd,
                         ArrayOfIndex&  input,
                         ArrayOfIndex&  auto_vars,
                         Array<TokVal>& auto_vars_values,
                         bool&          first);

  void parse_output(const MdRecord* mdd, ArrayOfIndex& output, bool& first);

  void parse_output_and_input(const MdRecord*&     mdd,
                                    Index&         id,
                                    String&        methodname,
                                    Array<TokVal>& values,
                                    ArrayOfIndex&  output,
                                    ArrayOfIndex&  input,
                                    ArrayOfIndex&  auto_vars,
                                    Array<TokVal>& auto_vars_values);

  void parse_keywords(const MdRecord*      mdd,
                            Array<TokVal>& values,
                      const bool           found_curly_brace);

  void parse_keywords2(const MdRecord*      mdd,
                             Array<TokVal>& values,
                             bool&          first);

  void tasklist_insert_set_delete(const ArrayOfIndex&  auto_vars,
                                  const Array<TokVal>& auto_vars_values,
                                  const Index          method_type,
                                        Agenda&        tasklist);

  bool is_whitespace(const char c);

  void eat_whitespace();

  void eat_whitespace_from_string(String& str, size_t& pos);

  void read_name(String& name);

  Index read_name_or_value(      String&        name,
                                 ArrayOfIndex&  auto_vars,
                                 Array<TokVal>& auto_vars_values,
                           const String&        default_name,
                           const MdRecord*      mdd,
                           const Index          group);

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

  Index mcfile_version;
};

#endif /* parser_h */

