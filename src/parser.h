/* Copyright (C) 2012 Oliver Lemke <olemke@core-dump.info>

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

#include "agenda_class.h"
#include "methods.h"
#include "sourcetext.h"
#include <map>

class ArtsParser {
 public:
  ArtsParser(Agenda& tasklist, String controlfile, const Verbosity& verbosity);

  void parse_tasklist();

 private:
  using NamedArgument = struct {
    String name;
    Index line;
    Index column;
  };

  using NamedArguments = vector<NamedArgument>;

  void find_named_arguments(vector<NamedArgument>& named_args);

  void skip_to_next_argument();

  void at_end_of_argument(const String& argname);

  void get_argument_index_by_name(Index& arg_index,
                                  NamedArguments& named_args,
                                  const String& name);

  void parse_main();

  void parse_agenda(Agenda& tasklist, const String& agenda_name);

  void parse_method(Index& id,
                    ArrayOfIndex& output,
                    ArrayOfIndex& input,
                    Agenda& tasks,
                    ArrayOfIndex& auto_vars,
                    Array<TokVal>& auto_vars_values,
                    String& include_file,
                    bool no_eot = false);

  void parse_generic_input(const MdRecord*& mdd,
                           Index& id,
                           String& methodname,
                           ArrayOfIndex& input,
                           ArrayOfIndex& auto_vars,
                           Array<TokVal>& auto_vars_values,
                           bool& first,
                           bool& still_supergeneric,
                           String& supergeneric_args,
                           Index& supergeneric_index,
                           NamedArguments& named_args,
                           bool call_by_name);

  void parse_generic_output(const MdRecord*& mdd,
                            Index& id,
                            String& methodname,
                            ArrayOfIndex& output,
                            bool& first,
                            bool& still_supergeneric,
                            String& supergeneric_args,
                            Index& supergeneric_index,
                            NamedArguments& named_args,
                            bool call_by_name);

  void parse_specific_input(const MdRecord* mdd,
                            ArrayOfIndex& input,
                            ArrayOfIndex& auto_vars,
                            Array<TokVal>& auto_vars_values,
                            bool& first,
                            NamedArguments& named_args,
                            bool call_by_name);

  void parse_specific_output(const MdRecord* mdd,
                             ArrayOfIndex& output,
                             bool& first,
                             NamedArguments& named_args,
                             bool call_by_name);

  void parse_method_args(const MdRecord*& mdd,
                         Index& id,
                         String& methodname,
                         ArrayOfIndex& output,
                         ArrayOfIndex& input,
                         ArrayOfIndex& auto_vars,
                         Array<TokVal>& auto_vars_values);

  String set_gin_to_default(const MdRecord* mdd,
                            ArrayOfIndex& auto_vars,
                            Array<TokVal>& auto_vars_values,
                            Index keyword_index);

  void tasklist_insert_set_delete(const ArrayOfIndex& auto_vars,
                                  const Array<TokVal>& auto_vars_values,
                                  const Index method_type,
                                  Agenda& tasklist);

  bool is_whitespace(const char c);

  void eat_whitespace();

  void eat_whitespace_from_string(String& str, size_t& pos);

  void read_name(String& name);

  Index read_name_or_value(String& name,
                           ArrayOfIndex& auto_vars,
                           Array<TokVal>& auto_vars_values,
                           const String& default_name,
                           const MdRecord* mdd,
                           const Index group);

  void assertain_character(char c);

  void parse_String(String& res);

  void read_integer(String& res);

  void read_numeric(String& res);

  void parse_integer(Index& n);

  void parse_numeric(Numeric& n);

  void parse_Stringvector(ArrayOfString& res);

  void parse_intvector(ArrayOfIndex& res);

  void parse_numvector(Vector& res);

  void parse_matrix(Matrix& res);

  bool parse_intvector_from_string(ArrayOfIndex& res, String& str);

  bool parse_numvector_from_string(Vector& res, String& str);

  bool parse_stringarray_from_string(ArrayOfString& res, String& str);

  Agenda& mtasklist;

  std::shared_ptr<Workspace> ws;

  String mcfile;

  SourceText msource;

  Index mcfile_version;

  const Verbosity& verbosity;
};

#endif /* parser_h */
