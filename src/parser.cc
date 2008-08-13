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

#include <iostream>
#include "parser.h"
#include "arts.h"
#include "exceptions.h"
#include "file.h"
#include "methods.h"
#include "parameters.h"
#include "wsv_aux.h"
#include "workspace_ng.h"


/** Constructs a new parser.

    \param[out] tasklist    Method list read from the controlfile.
    \param[in]  controlfile Path to the controlfile.

    \author Oliver Lemke
*/
ArtsParser::ArtsParser(Agenda& tasklist, String controlfile)
  : mtasklist (tasklist), mcfile (controlfile), mcfile_version (1)
{
  msource.AppendFile (mcfile);
}


/** Public interface to the main function of the parser.

    \author Oliver Lemke
*/
void ArtsParser::parse_tasklist ()
{
  parse_main ();
}


/** The main function of the parser. This will parse the entire
    text. FIXME: Add more documentation here.

    \author Stefan Buehler, Oliver Lemke
*/
void ArtsParser::parse_main()
{
  try 
    {
      extern const Array<MdRecord> md_data;

      // For method ids:
      Index id;         
      // For keyword parameter values:
      ArrayOfIndex values;
      // Output workspace variables (for generic methods):
      ArrayOfIndex output;              
      // Input workspace variables (for generic methods):
      ArrayOfIndex input;
      // For include statements, holding the include file's name
      String include_file;

      ArrayOfIndex  auto_vars;
      Array<TokVal> auto_vars_values;

      out3 << "\nParsing control text:\n";

      msource.Init();
      eat_whitespace();

      parse_method(id,values,output,input,mtasklist,auto_vars,auto_vars_values,
                   include_file,true);
          
      if ( "Arts" != md_data[id].Name() &&  "Arts2" != md_data[id].Name() )
        {
          out0 << "The outermost method must be Arts!\n"
               << "(But it seems to be " << md_data[id].Name() << ".)\n";
          arts_exit ();
        }
    }
  catch (const Eot x)
    {
      // Unexpected end of the source text:
      out0 << "Unexpected end of control script.\n";
      out0 << "File: " << x.file() << '\n';
      out0 << "Line: " << x.line() << '\n';
      arts_exit ();
    }
  catch (const UnexpectedChar x)
    {
      // Unexpected Character:
      out0 << "Unexpected character:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      arts_exit ();
    }
  catch (const IllegalLinebreak x)
    {
      // A line break in an illegal position:
      out0 << "Illegal Line break:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      arts_exit ();
    }
  catch (const UnknownMethod x)
    {
      // Method unknown:
      // [**This should give a hint on how to obtain a list of allowed 
      // methods.]
      out0 << "Unknown Method:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      arts_exit ();
    }
  catch (const UnknownWsv x)
    {
      // Workspace variable unknown:
      // [**This should give a hint on how to obtain a list of allowed 
      // Wsvs.]
      out0 << "Unknown workspace variable:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      arts_exit ();
    }
  catch (const WrongWsvGroup x)
    {
      // Workspace variable unknown:
      // [**This should give a hint on how to obtain a list of Wsvs in 
      // this group.
      out0 << "Workspace variable belongs to the wrong group:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      arts_exit ();
    }
  catch (const UnexpectedKeyword x)
    {
      // Keyword unknown:
      // [**This should give a hint on how to obtain a list of allowed 
      // keywords.]
      out0 << "Unknown keyword:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      arts_exit ();
    }
  catch (const ParseError x)
    {
      // General Parse Error (parent of all the above):
      out0 << "Parse error:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      arts_exit ();
    }
}


/** Parse the Contents of text as ARTS control input. 
  
    This method is used to parse the list of methods given in the
    curly braces of an agenda method. So the end is marked by a
    closing curly brace ahead.

    \param[out] tasklist The ids and keyword parameter values for the
                         methods to run.
  
    \see eat_whitespace
    \see parse_method
          
    \author Stefan Buehler
*/
void ArtsParser::parse_agenda( Agenda& tasklist )
{
  extern const Array<MdRecord> md_data;

  // For method ids:
  Index id;             
  // For keyword parameter values:
  ArrayOfIndex values;
  // Output workspace variables:
  ArrayOfIndex output;          
  // Input workspace variables:
  ArrayOfIndex input;
  // For Agenda, if there is any:
  Agenda tasks;
  // For include statements, holding the include file's name
  String include_file;

  ArrayOfIndex  auto_vars;
  Array<TokVal> auto_vars_values;

  eat_whitespace();

  while ( '}' != msource.Current() )
    {
      parse_method(id,values,output,input,tasks,
                   auto_vars,auto_vars_values,include_file);

      // If parse_method found an include statement it returnes -1 for the
      // method id
      if (id == -1)
        {
          // Command line parameters which give us the include search path.
          extern const Parameters parameters;

          if (!find_file (include_file, ".arts"))
            {
              ostringstream os;
              os << "Cannot find include file " << include_file
                << ".\n";
              os << "File: "   << msource.File() << '\n';
              os << "Search path was: . " << parameters.includepath
                << "\n";
              throw runtime_error (os.str());
            }

          out3 << "- Including control file " << include_file << "\n";

          ArtsParser include_parser (tasks, include_file);
          include_parser.parse_tasklist();

          for (Index i = 0; i < tasks.nelem(); i++)
            tasklist.push_back (tasks.Methods()[i]);
        }
      else
        {
          if (md_data[id].SetMethod())
            {
              // Append task to task list:      
              tasklist.push_back(MRecord(id,values,output,input,
                                         auto_vars_values[0],tasks));
            }
          else
            {
              tasklist_insert_set_delete(auto_vars, auto_vars_values,
                                         0, tasklist);

              // Append task to task list:      
              tasklist.push_back(MRecord(id,values,output,input,TokVal(),tasks));

              tasklist_insert_set_delete(auto_vars, auto_vars_values,
                                         1, tasklist);
            }

          {
            // Everything in this block is just to generate some
            // informative output.  

            out3 << "- " << md_data[id].Name() << "\n";

            /* FIXME OLE */
            /*for ( Index j=0 ; j<values.nelem() ; ++j )
              {
                out3 << "   " 
                  << md_data[id].Keywords()[j] << ": "
                  << values[j] << '\n';
              }*/

            // Output workspace variables for generic methods:
            if ( 0 < md_data[id].GOutType().nelem()
                 + md_data[id].GInType().nelem() )
              {
                out3 << "   Output: ";
                for ( Index j=0 ; j<output.nelem() ; ++j )
                  {
                    out3 << Workspace::wsv_data[output[j]].Name() << " ";
                  }
                out3 << "\n";

                out3 << "   Input: ";
                for ( Index j=0 ; j<input.nelem() ; ++j )
                  {
                    out3 << Workspace::wsv_data[input[j]].Name() << " ";
                  }
                out3 << "\n";
              }
          }
        }

      eat_whitespace();
    } 
}


/** Parse the Contents of text as ARTS control input. 

  Either values or tasks will be empty.

  \param[out]    id           Method id.
  \param[out]    values       Keyword parameter values for this method.
  \param[out]    output       Output workspace variables (for generic methods).
  \param[out]    input        Input workspace variables (for generic methods).
  \param[out]    tasks        A list of other methods.
  \param[out]    auto_vars    Indexes of automatically created variables.
  \param[out]    auto_vars_values    Values of automatically created variables.
  \param[out]    include_file The input to parse.
  \param[in]     no_eot       Suppress throwing an error on EOT after the
                              closing curly brace.
   
  \see read_name
  \see eat_whitespace
  \see assertain_character
  \see parse_String
  \see parse_integer
  \see parse_numeric
  \see parse_Stringvector
  \see parse_intvector
  \see parse_numvector
   
  \exception UnknownMethod
  \exception UnknownWsv
  \exception WrongWsvGroup
  \exception UnexpectedKeyword

  \author Stefan Buehler
*/
void ArtsParser::parse_method(Index& id, 
                              ArrayOfIndex&  values,
                              ArrayOfIndex&  output,
                              ArrayOfIndex&  input,
                              Agenda&        tasks,
                              ArrayOfIndex&  auto_vars,
                              Array<TokVal>& auto_vars_values,
                              String&        include_file,
                              bool no_eot)
{
  String methodname;            // We need this out here, since it is
                                // set once and later modified.

  const MdRecord* mdd;          // Handle on the method record. Needed here,
                                // because it is modified.

  bool found_curly_brace = false;

  // Clear all output variables:
  id = 0;
  values.resize( 0 );
  output.resize( 0 );
  input.resize(  0 );
  tasks.resize(  0 );
  auto_vars.resize(  0 );
  auto_vars_values.resize(  0 );
  include_file = "";
  
  read_name(methodname);

  if (methodname == "INCLUDE")
    {
      String s;
      eat_whitespace ();
      parse_String (include_file);

      id = -1;

      return;
    }
  else
    {
      if (methodname == "Arts2")
        {
          mcfile_version = 2;
          methodname = "Arts";
        }

      eat_whitespace();

      parse_output_and_input(mdd, id, methodname, values, output, input,
                             auto_vars, auto_vars_values);

      eat_whitespace();

      // Now look for the curly braces:
      if (msource.Current() == '{')
        {
          msource.AdvanceChar();
          eat_whitespace();
          found_curly_brace = true;
        }

      // There are two kind of methods, agenda methods, which have other
      // methods in the body, and normal methods, expecting keywords and
      // values. Let's take the agenda case first...
      if ( mdd->AgendaMethod() )
        {
          out3 << "- " << mdd->Name() << "\n";
          out3 << "{\n";
          parse_agenda(tasks);
          out3 << "}\n";
        }

      // For controlfiles v1 parse the keywords (or set them to default
      // values)
      if (mcfile_version == 1 && !mdd->AgendaMethod())
        {
          parse_keywords(mdd, values, auto_vars, auto_vars_values,
                         found_curly_brace);
        } 

      // Curly braces in non-agenda methods are not valid in v2 controlfiles
      if (mcfile_version == 2 && !mdd->AgendaMethod() && found_curly_brace)
        {
          ostringstream os;
          os << "Expected method name , but got `" << msource.Current() << "'.";
          throw UnexpectedChar( os.str(),
                                msource.File(),
                                msource.Line(),
                                msource.Column() );os << "" << endl;
        }
    }

  // Now look for the closing curly braces.  We have to catch Eot,
  // because after a method description may be a good place to end
  // the control file.
  if (found_curly_brace)
    {
      try
        {
          assertain_character('}');
        }
      catch (const Eot x)
        {
          // Re-throw the error if the no_eot flag is not set:
          if (!no_eot) throw Eot(x);
        }
    }
}


//! Parse keywords.
/** Parses the keyword parameters from a v1 controlfile.
 
  \param[in]  mdd               Method object.
  \param[out] values            Parsed keyword values.
  \param[in]  found_curly_brace Tells this method if the parser found an
                                opening curly brace before. If not, we just
                                fill all keywords with their default values.

  \author Oliver Lemke
  \date   2008-03-05
*/
void ArtsParser::parse_keywords(const MdRecord* mdd,
                                ArrayOfIndex&   values,
                                ArrayOfIndex&   auto_vars,
                                Array<TokVal>&  auto_vars_values,
                                const bool      found_curly_brace)
{
  // Now we have to deal with two different cases: Keywords with
  // parameters, or (optionally) only a parameter without a keyword
  // for methods that have only a single argument.
  //
  // We can distinguish the two cases if we check whether the current
  // character is a letter. (If the parameter is specified directly it
  // must be either a number, a +- sign or a quotation mark)
  //
  // KEYWORDS THAT START WITH A NUMBER WILL BREAK THIS CODE!!

  if (!found_curly_brace)
    {
      bool all_kw_have_defaults = true;
      for (Index kw = 0;
           all_kw_have_defaults && kw < mdd->Keywords().nelem();
           ++kw)
        {
          if (mdd->Defaults()[kw] == NODEF)
            all_kw_have_defaults = false;
        }

      if (!all_kw_have_defaults)
        {
          ostringstream os;
          os << "Expected '{', but got `" << msource.Current() << "'.";
          throw UnexpectedChar( os.str(),
                                msource.File(),
                                msource.Line(),
                                msource.Column() );
        }
    }
  values.resize (mdd->Keywords().nelem());

  // Use this array to remember which keywords have been set
  ArrayOfIndex initialized_keywords (mdd->Keywords().nelem(), 0);

  while (found_curly_brace && msource.Current () != '}')
    {
      Index keyword_index;
      if (!isalpha(msource.Current()) && 1==mdd->Keywords().nelem())
        {
          keyword_index = 0;
          initialized_keywords[keyword_index] = 1;
          // Parameter specified directly, without a keyword. This is only
          // allowed for single parameter methods!

          // We don't have to do anything here.
        }
      else  
        {      // Look for the keywords and read the parameters:

          bool found_keyword = false;
          String keyname;
          read_name(keyname);

          // Is the keyname the expected keyname?
          for (keyword_index = 0; keyword_index < mdd->Keywords().nelem();
               ++keyword_index)
            {
              if ( keyname == mdd->Keywords()[keyword_index] )
                {
                  found_keyword = true;
                  if (initialized_keywords[keyword_index])
                    {
                      ostringstream os;
                      os << "Keyword " << mdd->Keywords()[keyword_index]
                        << " was already set.\n";
                      throw ParseError (os.str (),
                                        msource.File(),
                                        msource.Line(),
                                        msource.Column());
                    }
                  else
                    {
                      initialized_keywords[keyword_index] = 1;
                    }
                  break;
                }
            }
          if ( !found_keyword )
            {
              throw UnexpectedKeyword( keyname,
                                       msource.File(),
                                       msource.Line(),
                                       msource.Column());
            }

          eat_whitespace();

          // Look for '='
          assertain_character('=');
          eat_whitespace();
        }

        {
          // Now parse the key value.
          String wsvname;
          ostringstream os;
          os << keyword_index;
          read_name_or_value(wsvname,
                             auto_vars, auto_vars_values,
                             "keyword" + os.str() + "_"
                             + mdd->Keywords()[keyword_index],
                             mdd, mdd->Types()[keyword_index]);

            {
              // Find Wsv id:
              const map<String, Index>::const_iterator wsvit =
                Workspace::WsvMap.find(wsvname);
              if ( wsvit == Workspace::WsvMap.end() )
                throw UnknownWsv( wsvname,
                                  msource.File(),
                                  msource.Line(),
                                  msource.Column() );
              values[keyword_index] = wsvit->second;
            }

      }

      eat_whitespace();

      // Check:
      //      cout << "Value: " << mdd->Values()[i] << '\n';
    }

  // Check if all keywords are set
  for (Index i = 0; i < initialized_keywords.nelem(); ++i)
    {
      if (!initialized_keywords[i])
        {
          set_keyword_to_default(mdd, values, auto_vars, auto_vars_values, i);
        }
    }
}


//! Set keyword to default value.
/** Sets the value of the keyword with the given index to its default value.
 
  \param[in]      mdd           Method object.
  \param[out]     values        Keyword values.
  \param[in]      keyword_index Index of the keyword which should be set to
                                its default value.

  \author Oliver Lemke
  \date   2008-07-16
*/
void ArtsParser::set_keyword_to_default(const MdRecord* mdd,
                                        ArrayOfIndex&   values,
                                        ArrayOfIndex&   auto_vars,
                                        Array<TokVal>&  auto_vars_values,
                                        Index           keyword_index)
{
  if (mdd->Defaults()[keyword_index] != NODEF)
    {
      TokVal tv;
      // Now parse the key value. This can be:
      // String, Index, Numeric, ArrayOfString, ArrayOfIndex, Vector
      bool failed = false;
      if (mdd->Types()[keyword_index] == get_wsv_group_id ("String")) 
        {
          tv = mdd->Defaults()[keyword_index];
        }
      else if (mdd->Types()[keyword_index] == get_wsv_group_id ("Index")) 
        {
          Index n;
          istringstream is(mdd->Defaults()[keyword_index]);
          is >> n;
          tv = n;
          if (is.bad () || is.fail ())
            failed = true;
        }
      else if (mdd->Types()[keyword_index] == get_wsv_group_id ("Numeric")) 
        {
          Numeric n;
          istringstream is(mdd->Defaults()[keyword_index]);
          is >> n;
          tv = n;
          if (is.bad () || is.fail ())
            failed = true;
        }
      else if (mdd->Types()[keyword_index] == get_wsv_group_id ("ArrayOfString")) 
        {
          ArrayOfString v;
          String s = mdd->Defaults()[keyword_index];
          if (!parse_stringarray_from_string(v, s))
            {
              failed = true;
            }
          tv = v;
        }
      else if (mdd->Types()[keyword_index] == get_wsv_group_id ("ArrayOfIndex")) 
        {
          ostringstream os;
          os << "Default values for keywords with type "
            << "ArrayOfIndex are not supported.\n"
            << "Either remove the default value for keyword '"
            << mdd->Keywords()[keyword_index] << "' in workspace method *"
            << mdd->Name() << "* in methods.cc or discuss this "
            << "issue on the arts-dev mailing list.\n";
          throw ParseError (os.str(),
                            msource.File(),
                            msource.Line(),
                            msource.Column());
        }
      else if (mdd->Types()[keyword_index] == get_wsv_group_id ("Vector")) 
        {
          Vector v;
          String s = mdd->Defaults()[keyword_index];
          if (!parse_numvector_from_string(v, s))
            {
              failed = true;
            }
          tv = v;
        }
      else
        {
          failed = true;
        }

      Index wsvid;
      String name;

      {
          ostringstream os;
          os << keyword_index;

          name = "auto_" + mdd->Name() + "_" + "keyword" + os.str() + "_"
            + mdd->Keywords()[keyword_index];
      }

      map<String, Index>::const_iterator wsvit =
        Workspace::WsvMap.find(name);
      if (wsvit == Workspace::WsvMap.end())
        {
          wsvid = Workspace::add_wsv(
            WsvRecord(name.c_str(),
                      "Automatically allocated variable.",
                      mdd->Types()[keyword_index],
                      true));
        }
      else
        {
          wsvid = wsvit->second;
        }

      auto_vars.push_back(wsvid);
      auto_vars_values.push_back(tv);
      values[keyword_index] = wsvid;

      if (failed)
        {
          ostringstream os;
          os << "Failed to assign default value for keyword '"
            << mdd->Keywords()[keyword_index] << "'.\n"
            << "Check definition of workspace method *"
            << mdd->Name() << "* in methods.cc.\n";
          throw ParseError (os.str(),
                            msource.File(),
                            msource.Line(),
                            msource.Column());
        }
    }
  else
    {
      ostringstream os;
      os << "Keyword '" << mdd->Keywords()[keyword_index]
        << "' omitted but no default value found.\n"
        << "Check definition of workspace method *"
        << mdd->Name() << "* in methods.cc.\n";
      throw ParseError (os.str(),
                        msource.File(),
                        msource.Line(),
                        msource.Column());
    }
}


//! Parse keywords.
/** Parses the keyword parameters from a v2 controlfile.
 
  \param[in]      mdd     Method object.
  \param[out]     values  Parsed keyword values.
  \param[in,out]  first   Tells this method if there were any other arguments
                          passed to this method before. If so, we have to look
                          for a comma first.

  \author Oliver Lemke
  \date   2008-03-05
*/
void ArtsParser::parse_keywords2(const MdRecord* mdd,
                                 ArrayOfIndex&   values,
                                 ArrayOfIndex&   auto_vars,
                                 Array<TokVal>&  auto_vars_values,
                                 bool&           first)
{
  Index keyword_index = 0;

  eat_whitespace();

  values.resize(mdd->Keywords().nelem());

  while (msource.Current() != ')' && keyword_index < mdd->Keywords().nelem())
    {
      if (!first)
        {
          assertain_character(',');
          eat_whitespace();
        }
      else
        first = false;

      // If there is a comma or a closing brace means no value was
      // specified and we use the default value instead (if there is one)
      if (msource.Current() == ',' || msource.Current() == ')')
        {
          set_keyword_to_default(mdd, values, auto_vars, auto_vars_values,
                                 keyword_index);
        }
      else
        {
          String wsvname;
          ostringstream os;
          os << keyword_index;
          read_name_or_value(wsvname,
                             auto_vars, auto_vars_values,
                             "keyword" + os.str() + "_"
                             + mdd->Keywords()[keyword_index],
                             mdd, mdd->Types()[keyword_index]);

          {
            // Find Wsv id:
            const map<String, Index>::const_iterator wsvit =
              Workspace::WsvMap.find(wsvname);
            if ( wsvit == Workspace::WsvMap.end() )
              throw UnknownWsv( wsvname,
                                msource.File(),
                                msource.Line(),
                                msource.Column() );
            values[keyword_index] = wsvit->second;
          }

        }

      eat_whitespace();
      keyword_index++;
    }

  if (keyword_index < mdd->Keywords().nelem())
    {
      ostringstream os;
      os << "Expected keyword '" << mdd->Keywords()[keyword_index]
        << "', but got `" << msource.Current() << "'.";
      throw UnexpectedChar( os.str(),
                            msource.File(),
                            msource.Line(),
                            msource.Column() );
    }
  else if (msource.Current() != ')')
    {
      ostringstream os;
      os << "Expected ')', but got '" << msource.Current() << "'.";
      throw UnexpectedChar( os.str(),
                            msource.File(),
                            msource.Line(),
                            msource.Column() );
    }

  msource.AdvanceChar();
}


//! Parse method's argument list.
/** This function parses a method's argument list for v1 and v2 controlfiles.

  \param[in]  mdd        Method
  \param[out] id         Method index in md_data.
  \param[out] methodname Name of the WSM
  \param[out] values     Keyword values
  \param[out] output     Output WSVs
  \param[out] input      Input WSVs
  \param[out] auto_vars  Indexes of automatically created variables.
  \param[out] auto_vars_values    Values of automatically created variables.

  \author Oliver Lemke
  \date   2008-03-05
*/
void ArtsParser::parse_output_and_input(const MdRecord*& mdd,
                                        Index&           id,
                                        String&          methodname,
                                        ArrayOfIndex&    values,
                                        ArrayOfIndex&    output,
                                        ArrayOfIndex&    input,
                                        ArrayOfIndex&    auto_vars,
                                        Array<TokVal>&   auto_vars_values)
{
  extern const map<String, Index> MdRawMap;
  extern const Array<MdRecord> md_data;
  extern const Array<MdRecord> md_data_raw;
  extern const ArrayOfString wsv_group_names;
  extern const map<String, Index> MdMap;

  Index wsvid;                  // Workspace variable id, is used to
                                // access data in wsv_data.

  bool still_supergeneric=true; // Flag that our MdRecord still is
                                // from md_data_raw, not from
                                // md_data. 

                                // Find method raw id in raw map:
  const map<String, Index>::const_iterator i = MdRawMap.find(methodname);
  if ( i == MdRawMap.end() )
    throw UnknownMethod(methodname,
                        msource.File(),
                        msource.Line(),
                        msource.Column());

  id = i->second;

  // Get a convenient handle on the data record for this method. We
  // have to use a pointer here, not a reference, because we later
  // want to change where mdd is pointing!
  mdd = &md_data_raw[id];

  // Is this a supergeneric method? If not, take the record in
  // md_data, rather than in md_data_raw:
  if ( !mdd->Supergeneric() )
    {
      // Find explicit method id in MdMap:
      const map<String, Index>::const_iterator i2 = MdMap.find(methodname);
      assert ( i2 != MdMap.end() );
      id = i2->second;            

      mdd = &md_data[id];

      still_supergeneric = false;
    }

  // For generic methods in v1 controlfiles the output and input workspace
  // variables have to be parsed (given in round brackets).
  // For v2 controlfiles, we also have to parse the keywords which are now
  // part of the argument list in parenthesis.
  if ( (mdd->GOutType().nelem() + mdd->GInType().nelem() > 0)
       || (mcfile_version == 2
           && 0 < mdd->GOutType().nelem() + mdd->GInType().nelem()))
    {
      String wsvname;
      bool first = true;        // To skip the first comma.

      assertain_character('(');
      eat_whitespace();

      if (mcfile_version == 2 && mdd->Output().nelem())
        {
          // Parse all output variables.
          parse_output(mdd, output, first);
        }
      else
        {
          // To keep v1 controlfiles compatible, we have to add the static
          // outputs to the output variable list and just hardwire them to the
          // WSVs given in the method's definition in methods.cc.
          ArrayOfIndex vo=mdd->Output();
          for (ArrayOfIndex::const_iterator outs=vo.begin(); outs<vo.end(); ++outs)
            {
              output.push_back (*outs);
            }
        }

      // Parse all generic output variables
      for ( Index j=0 ; j<mdd->GOutType().nelem() ; ++j )
        {
          if (first)
            first = false;
          else
            {
              assertain_character(',');
              eat_whitespace();
            }

          read_name(wsvname);

          {
            wsvid = -1;
            // Find Wsv id:
            map<String, Index>::const_iterator wsvit =
              Workspace::WsvMap.find(wsvname);
            if ( wsvit == Workspace::WsvMap.end() )
              {
                if (still_supergeneric)
                  {
                    ostringstream os;
                    os << "This might be either a typo or you have to create "
                      << "the variable\nby calling TYPECreate(" << wsvname
                      << ") first. Replace TYPE with the\n"
                      << "WSV group your variable should belong to.";

                    throw UnknownWsv( os.str(),
                                      msource.File(),
                                      msource.Line(),
                                      msource.Column() );
                  }
                else
                  {
                    if (mdd->Name().length() > 6
                        && mdd->Name().substr (mdd->Name().length() - 6)
                        != "Create")
                      {
                        ostringstream os;
                        os << "This might be either a typo or you have to create "
                          << "the variable\nby calling "
                          << wsv_group_names[mdd->GOutType()[j]]
                          << "Create( " << wsvname
                          << " ) first.\n";

                        throw UnknownWsv( os.str(),
                                          msource.File(),
                                          msource.Line(),
                                          msource.Column() );
                      }
                    else
                      {
                        wsvid = Workspace::add_wsv(
                          WsvRecord(wsvname.c_str(),
                                    "Automatically allocated variable.",
                                    mdd->GOutType()[j],
                                    true));
                      }
                  }
              }

            if (wsvid == -1)
              wsvid = wsvit->second;
          }

          // If this is a supergeneric method, now is the time to find
          // out the actual group of the argument(s)!
          if ( still_supergeneric )
            {
//            cout << "wsvid = " << wsvid << "\n";
//            cout << "wsv_group_names[Workspace::wsv_data[wsvid].Group()] = "
//                 << wsv_group_names[Workspace::wsv_data[wsvid].Group()] << "\n";
              ostringstream os;
              os << mdd->Name() << "_sg_"
                 << wsv_group_names[Workspace::wsv_data[wsvid].Group()];
              methodname = os.str();

              // Find explicit method id in MdMap:
              const map<String, Index>::const_iterator mdit = MdMap.find(methodname);
              assert ( mdit != MdMap.end() );
              id = mdit->second;         
              
              mdd = &md_data[id];

              still_supergeneric = false;

//            cout << "Adjusted id=" << id << '\n';   
//            cout << "Adjusted Method: " << mdd->Name() << '\n';
            }

          // Now we have explicitly the method record for the right
          // group. From now on no special treatment of supergeneric
          // methods should be necessary.

          // Check that this Wsv belongs to the correct group:
          if ( Workspace::wsv_data[wsvid].Group() != mdd->GOutType()[j] )
            {
            throw WrongWsvGroup( wsvname+" is not "+
                        wsv_group_names[mdd->GOutType()[j]]+", it is "+ 
                        wsv_group_names[Workspace::wsv_data[wsvid].Group()],
                                 msource.File(),
                                 msource.Line(),
                                 msource.Column() );
            }

          // Add this one to the list of output workspace variables:
          output.push_back(wsvid);
          
          eat_whitespace();
        }

      if (mcfile_version == 2 && mdd->Input().nelem())
        {
          // Parse all input variables
          parse_input(mdd, input, auto_vars, auto_vars_values, first);
          first = false;
        }
      else
        {
          // To keep v1 controlfiles compatible, we have to add the static
          // inputs to the input variable list and just hardwire them to the
          // WSVs given in the method's definition in methods.cc.
          ArrayOfIndex vi;
          mdd->input_only(vi);
          for (ArrayOfIndex::const_iterator ins=vi.begin(); ins<vi.end(); ++ins)
            {
              input.push_back (*ins);
            }
        }

      // Then parse all generic input variables
      for ( Index j=0 ; j<mdd->GInType().nelem() ; ++j )
        {
          if (first)
            first = false;
          else
            {
              assertain_character(',');
              eat_whitespace();
            }

          {
            ostringstream os;
            os << j;
            read_name_or_value(wsvname, auto_vars, auto_vars_values,
                               "generic" + os.str(),
                               mdd, mdd->GInType()[j]);
          }

          {
            // Find Wsv id:
            const map<String, Index>::const_iterator wsvit =
              Workspace::WsvMap.find(wsvname);
            if ( wsvit == Workspace::WsvMap.end() )
              throw UnknownWsv( wsvname,
                                msource.File(),
                                msource.Line(),
                                msource.Column() );

            wsvid = wsvit->second;
          }

          // Is the method data record still supergeneric? This could
          // be the case if there are no output arguments, only input
          // arguments. In that case, let's find out the actual group!
          if ( still_supergeneric )
            {
              ostringstream os;
              os << mdd->Name() << "_sg_"
                 << wsv_group_names[Workspace::wsv_data[wsvid].Group()];
              methodname = os.str();

              // Find explicit method id in MdMap:
              const map<String, Index>::const_iterator mdit =
                MdMap.find(methodname);
              assert ( mdit != MdMap.end() );
              id = mdit->second;         

              mdd = &md_data[id];

              still_supergeneric = false;
            }

          // Now we have explicitly the method record for the right
          // group. From now on no special treatment of supergeneric
          // methods should be necessary.

          // Check that this Wsv belongs to the correct group:
          if ( Workspace::wsv_data[wsvid].Group() != mdd->GInType()[j] )
            throw WrongWsvGroup( wsvname+" is not "+
                        wsv_group_names[mdd->GInType()[j]]+", it is "+ 
                        wsv_group_names[Workspace::wsv_data[wsvid].Group()],
                                 msource.File(),
                                 msource.Line(),
                                 msource.Column() );

          // Add this one to the list of input variables:
          input.push_back(wsvid);
          
          eat_whitespace();
        }

      if (mcfile_version == 2)
        {
          // For v2 controlfiles we also parse the keywords
          parse_keywords2(mdd, values, auto_vars, auto_vars_values, first);
        }
      else
        {
          assertain_character(')');
          eat_whitespace();
        }
    }
  else
    {
      // Even if the method has no GInType or GOutType the user can
      // pass all inputs and outputs in parenthesis
      if (msource.Current() == '(')
        {
          bool first = true;
          msource.AdvanceChar();
          eat_whitespace();
          if (mcfile_version == 2)
            {
              if (mdd->Output().nelem())
                {
                  parse_output(mdd, output, first);
                  if (mdd->Input().nelem())
                    parse_input(mdd, input, auto_vars, auto_vars_values, first);
                }
              else if (mdd->Input().nelem())
                parse_input(mdd, input, auto_vars, auto_vars_values, first);

              parse_keywords2(mdd, values, auto_vars, auto_vars_values, first);
            }

          if (mcfile_version == 1)
            {
              assertain_character(')');
              eat_whitespace();
            }
        }
      else
        {
          // If the parenthesis were omitted we still have to add the implicit
          // outputs and inputs to the methods input and output variable lists
          ArrayOfIndex vo=mdd->Output();
          for (ArrayOfIndex::const_iterator outs=vo.begin();
               outs<vo.end(); ++outs)
            {
              output.push_back (*outs);
            }
          ArrayOfIndex vi;
          mdd->input_only(vi);
          for (ArrayOfIndex::const_iterator ins=vi.begin(); ins<vi.end(); ++ins)
            {
              input.push_back (*ins);
            }

          if (mcfile_version == 2)
            {
              // Make sure all keywords have default values, otherwise the
              // user has to specify them in the controlfile.
              bool all_kw_have_defaults = true;
              values.resize (mdd->Keywords().nelem());
              for (Index kw = 0;
                   all_kw_have_defaults && kw < mdd->Keywords().nelem();
                   ++kw)
                {
                  if (mdd->Defaults()[kw] == NODEF)
                    all_kw_have_defaults = false;
                  else
                    set_keyword_to_default(mdd, values,
                                           auto_vars, auto_vars_values, kw);
                }

              if (!all_kw_have_defaults)
                {
                  ostringstream os;
                  os << "Not all keywords of this method have default values, you have "
                    << "to specify them!";
                  throw UnexpectedChar( os.str(),
                                        msource.File(),
                                        msource.Line(),
                                        msource.Column() );
                }
            }

        }
    }
}


/** Parse the input WSVs for current method from the controlfile.

  \param[in]  mdd    Pointer to the current WSM
  \param[out] input  Indexes of input variables for the WSM
  \param[in]  first  If set to false, there must be a comma before the first WSV
                     in the controlfile

  \author Oliver Lemke
  \date   2008-03-05
  */
void ArtsParser::parse_input(const MdRecord*      mdd,
                                   ArrayOfIndex&  input,
                                   ArrayOfIndex&  auto_vars,
                                   Array<TokVal>& auto_vars_values,
                                   bool&          first)
{
  extern const ArrayOfString wsv_group_names;

  // There are two lists of parameters that we have to read.
  ArrayOfIndex  vo=mdd->Output();  // Output 
  ArrayOfIndex  vi;                // Input

  Index wsvid;  // Workspace variable id, is used to
                // access data in wsv_data.

  mdd->input_only(vi);             // Input only variables

  for (ArrayOfIndex::const_iterator ins=vi.begin(); ins<vi.end(); ++ins)
    {
      String wsvname;

      if (first)
        first = false;
      else
        {
          try
            {
              assertain_character(',');
            }
          catch (UnexpectedChar)
            {
              ostringstream os;
              os << "Expected input WSV *" << Workspace::wsv_data[*ins].Name() << "*";
            }
          eat_whitespace();
        }

      read_name_or_value(wsvname, auto_vars, auto_vars_values,
                         Workspace::wsv_data[*ins].Name(),
                         mdd, Workspace::wsv_data[*ins].Group());

      {
        // Find Wsv id:
        const map<String, Index>::const_iterator wsvit =
          Workspace::WsvMap.find(wsvname);
        if ( wsvit == Workspace::WsvMap.end() )
          throw UnknownWsv( wsvname,
                            msource.File(),
                            msource.Line(),
                            msource.Column() );

        wsvid = wsvit->second;
      }

      // Check that this Wsv belongs to the correct group:
      if ( Workspace::wsv_data[wsvid].Group() != Workspace::wsv_data[*ins].Group() )
        {
          throw WrongWsvGroup( wsvname+" is not "+
                               wsv_group_names[Workspace::wsv_data[*ins].Group()]+", it is "+ 
                               wsv_group_names[Workspace::wsv_data[wsvid].Group()],
                               msource.File(),
                               msource.Line(),
                               msource.Column() );
        }

      input.push_back(wsvid);
    }
}


/** Parse the output WSVs for current method from the controlfile.

  \param[in]  mdd     Pointer to the current WSM
  \param[out] output  Indexes of output variables for the WSM
  \param[in]  first   If set to false, there must be a comma before the first WSV
                      in the controlfile

  \author Oliver Lemke
  \date   2008-03-05
  */
void ArtsParser::parse_output(const MdRecord* mdd, ArrayOfIndex& output, bool& first)
{
  extern const ArrayOfString wsv_group_names;

  ArrayOfIndex  vo=mdd->Output();

  Index wsvid;                  // Workspace variable id, is used to
                                // access data in wsv_data.

  for (ArrayOfIndex::const_iterator outs=vo.begin(); outs<vo.end(); ++outs)
    {
      String wsvname;

      if (first)
        first = false;
      else
        {
          try
            {
              assertain_character(',');
            }
          catch (UnexpectedChar)
            {
              ostringstream os;
              os << "Expected output WSV *" << Workspace::wsv_data[*outs].Name() << "*";
              throw ParseError( os.str(),
                                msource.File(),
                                msource.Line(),
                                msource.Column() );
            }
          eat_whitespace();
        }

      read_name(wsvname);

      {
        wsvid = -1;
        // Find Wsv id:
        map<String, Index>::const_iterator wsvit =
          Workspace::WsvMap.find(wsvname);
        if ( wsvit == Workspace::WsvMap.end() )
          {
            if (mcfile_version == 1)
              {
                throw UnknownWsv( wsvname,
                                  msource.File(),
                                  msource.Line(),
                                  msource.Column() );
              }
            else
              {
                if (mdd->Name().length() > 6
                    && mdd->Name().substr (mdd->Name().length() - 6)
                    != "Create")
                  {
                    ostringstream os;
                    os << "This might be either a typo or you have to create "
                      << "the variable\nby calling "
                      << wsv_group_names[Workspace::wsv_data[*outs].Group()]
                      << "Create( " << wsvname
                      << " ) first.\n";

                    throw UnknownWsv( os.str(),
                                      msource.File(),
                                      msource.Line(),
                                      msource.Column() );
                  }
                else
                  {
                    wsvid = Workspace::add_wsv(
                      WsvRecord(wsvname.c_str(),
                                "Automatically allocated variable.",
                                Workspace::wsv_data[*outs].Group(),
                                true));
                  }
              }
          }

        if (wsvid == -1)
          wsvid = wsvit->second;
      }

      // Check that this Wsv belongs to the correct group:
      if ( Workspace::wsv_data[wsvid].Group() != Workspace::wsv_data[*outs].Group() )
        {
          throw WrongWsvGroup( wsvname+" is not "+
                               wsv_group_names[Workspace::wsv_data[*outs].Group()]+", it is "+ 
                               wsv_group_names[Workspace::wsv_data[wsvid].Group()],
                               msource.File(),
                               msource.Line(),
                               msource.Column() );
        }

      output.push_back(wsvid);
    }
}


/** Insert Set and Delete methods for automatically allocated output WSVs.

    This function inserts either a bunch of Set or Delete methods for
    implicitly allocated output WSVs. This needs to be done if the controlfile contains
    a value instead of an output variable name.

    \see ArtsParser::parse_agenda

    \param[in]  auto_vars        Indexes of automatically created variables.
    \param[in]  auto_vars_values Values of automatically created variables.
    \param[in]  method_type      0 = insert Set method, 1 = insert Delete method.
    \param[out] tasklist         Agenda to which the methods should be appended.
*/
void ArtsParser::tasklist_insert_set_delete(const ArrayOfIndex&  auto_vars,
                                            const Array<TokVal>& auto_vars_values,
                                            const Index          method_type,
                                                  Agenda&        tasklist)
{
  extern const map<String, Index> MdMap;
  extern const ArrayOfString wsv_group_names;

  for (Index i=0; i<auto_vars.nelem(); i++)
    {
      map<String, Index>::const_iterator mdit;
      Index         init_mdid;
      TokVal        auto_keyword_value;
      ArrayOfIndex  auto_output_var;
      ArrayOfIndex  auto_input_var;
      Agenda        auto_tasks;

      const Index auto_group = Workspace::wsv_data[auto_vars[i]].Group();
      if (auto_group != get_wsv_group_id("Index")
          && auto_group != get_wsv_group_id("Numeric")
          && auto_group != get_wsv_group_id("ArrayOfString")
          && auto_group != get_wsv_group_id("String")
          && auto_group != get_wsv_group_id("Vector"))
        {
          ostringstream os;
          os << "Passing a "
            << wsv_group_names[Workspace::wsv_data[auto_vars[i]].Group()]
            << " constant to a WSM is not supported!";
          throw ParseError (os.str(),
                            msource.File(),
                            msource.Line(),
                            msource.Column());
        }

      String method_name;
      switch (method_type)
        {
        case 0:
          auto_keyword_value = auto_vars_values[i];
          auto_output_var.push_back(auto_vars[i]);
          method_name = wsv_group_names[
            Workspace::wsv_data[auto_vars[i]].Group()] + "Set";
          break;
        case 1:
          auto_input_var.push_back(auto_vars[i]);
          method_name = "Delete_sg_" + wsv_group_names[Workspace::wsv_data[auto_vars[i]].Group()];
          break;
        default:
          throw runtime_error("Invalid method_type");
        }

      mdit = MdMap.find(method_name);
      assert ( mdit != MdMap.end() );
      init_mdid = mdit->second;         

      tasklist.push_back(MRecord(init_mdid, ArrayOfIndex(),
                                 auto_output_var, auto_input_var,
                                 auto_keyword_value,
                                 auto_tasks));
    }
}


/** Returns true if this character is considered whitespace. This
    includes the comment sign `\#'. This function is used by other
    functions to test for delimiting whitespace. 

    The whitespace cases implemented here must be consistent with
    eat_whitespace! 
    \see eat_whitespace  */
bool ArtsParser::is_whitespace(const char c)
{
  switch (c) 
    {
    case ' ':
    case '\r':
    case '\t':
    case '#':
      return true;
      break;
    }

  return false;
}


/** Eats whitespace. Comments are a special case of
    whitespace. Everything from the `\#' to the end of the line is
    eaten. 
  
    The whitespace cases implemented here must be consistent with
    is_whitespace! 
    \see is_whitespace    

    \exception UnexpectedChar Non-whitespace character encountered. */
void ArtsParser::eat_whitespace()
{
  char dummy;

  while (is_whitespace(dummy=msource.Current())) 
    {
      switch (dummy)
        {
        case ' ':
        case '\r':
        case '\t':
          msource.AdvanceChar();
          break;
        case '#':
          msource.AdvanceLine();
          break;
        default:
          {
            ostringstream os;
            os << "Expected whitespace, but got `" << dummy << "'.";
            throw UnexpectedChar( os.str(),
                                  msource.File(),
                                  msource.Line(),
                                  msource.Column() ); 
            break;
          }
        }
    }
}


/** Eats whitespace from a String.

  \param[in]     str String.
  \param[in,out] pos Current position in the String.
 */
void ArtsParser::eat_whitespace_from_string(String& str, size_t& pos)
{
  while (pos < str.length() && is_whitespace(str[pos]))
    pos++;
}


/** Reads name of method, keyword, or workspace variable.
  
    These names may consist only of letters (case matters!), numbers, and
    underscores.
    Line break or any other character ends the name.
  
    Whitespace has to have been eaten before. Scanns source for the
    name, starting at position specified by line and column.
    
    \param[out] name Method, keyword or WSV name
*/
void ArtsParser::read_name(String& name)
{
  bool stop = false;
  name = "";

  if (!isalpha(msource.Current()))
    {
      ostringstream os;
      os << "Workspace variable names must start with a letter!";
      throw ParseError( os.str(),
                        msource.File(),
                        msource.Line(),
                        msource.Column() );
    }

  while (!stop) 
    {
      char dummy = msource.Current();

      if ( isalnum(dummy) || '_'==dummy )
        {
          name += dummy;
          // AdvanceChar sets LineBreak if a line break occured.
          msource.LineBreak() = false;
          msource.AdvanceChar();
          if ( msource.LineBreak() ) stop = true;
        }
      else
        {
          stop = true;
        }
    }

  //  cout << "Name: " << name << '\n';
}


/** Reads name of a workspace variable or a value.
  
    These names may consist only of letters (case matters!), numbers, and
    underscores.
    Line break or any other character ends the name.
  
    Whitespace has to have been eaten before. Scanns source for the
    name, starting at position specified by line and column.

    \param[out] name  WSV name or value
    \param[in]  group Expected WSV group index
    
    \return -1 If a WSV name was found or the index of the newly created
               WSV
*/
Index ArtsParser::read_name_or_value(      String&        name,
                                           ArrayOfIndex&  auto_vars,
                                           Array<TokVal>& auto_vars_values,
                                     const String&        default_name,
                                     const MdRecord*      mdd,
                                     const Index          group)
{
  name = "";

  if (isalpha(msource.Current()))
    {
      read_name(name);
      return -1;
    }

  if (group == get_wsv_group_id ("Any")) 
    {
      ostringstream os;
      os << "Passing constants as supergeneric parameters is not supported.";
      throw ParseError (os.str (),
                        msource.File(),
                        msource.Line(),
                        msource.Column());
    }

  // If a value was given instead of a variable name, we create
  // a new variable in the workspace and fill it with the given
  // value

  Index wsvid;

  name = "auto_" + mdd->Name() + "_" + default_name;
  map<String, Index>::const_iterator wsvit = Workspace::WsvMap.find(name);
  if (wsvit == Workspace::WsvMap.end())
    {
      wsvid = Workspace::add_wsv(WsvRecord(name.c_str(),
                                           "Automatically allocated variable.",
                                           group,
                                           true));
    }
  else
    {
      wsvid = wsvit->second;
    }

  auto_vars.push_back(wsvid);

  // Now parse the value. This can be:
  // String_, Index_, Numeric_, Array_String_, Array_Index_, Vector_
  if (group == get_wsv_group_id("String"))
    {
      String dummy;
      parse_String(dummy);
      auto_vars_values.push_back(dummy);
    }
  else if (group == get_wsv_group_id("Index"))
    {
      Index n;
      parse_integer(n);
      auto_vars_values.push_back(n);
    }
  else if (group == get_wsv_group_id("Numeric"))
    {
      Numeric n;
      parse_numeric(n);
      auto_vars_values.push_back(n);
    }
  else if (group == get_wsv_group_id("ArrayOfString"))
    {
      ArrayOfString dummy;
      parse_Stringvector(dummy);
      auto_vars_values.push_back(dummy);
    }
  else if (group == get_wsv_group_id("ArrayOfIndex"))
    {
      ArrayOfIndex dummy;
      parse_intvector(dummy);
      auto_vars_values.push_back(dummy);
    }
  else if (group == get_wsv_group_id("Vector"))
    {
      Vector dummy;
      parse_numvector(dummy);
      auto_vars_values.push_back(dummy);
    }
  else
    {
      extern const ArrayOfString wsv_group_names;
      ostringstream os;
      os << "Unsupported parameter type: " << wsv_group_names[group];
      throw runtime_error (os.str());
    }

  return wsvid;
}


/** Make sure that the current character is equal to c and go to the
    next character.
  
    \exception UnexpectedChar The character is not right. */
void ArtsParser::assertain_character(char c)
{
  if ( c != msource.Current() )
    {
      ostringstream os;
      os << "Expected `" << c << "', but got `" << msource.Current() << "'.";
      throw UnexpectedChar( os.str(),
                            msource.File(),
                            msource.Line(),
                            msource.Column() ); 
    }
  
  msource.AdvanceChar();
}


/** Reads a String, complete with quotation marks. Whitespace has to
    have been eaten before, that is, the current character must be
    the quotation mark (").
    Quotation marks inside Strings are currently not possible.
  
    Line breaks inside Strings are not allowed. 

    \param[out] res Output string.

    \exception IllegalLinebreak An illegal linebreak has occured. */
void ArtsParser::parse_String(String& res)
{
  bool stop = false;
  res = "";

  msource.LineBreak() = false;
  assertain_character('"');
  if ( msource.LineBreak() )
    throw IllegalLinebreak( "Line break before end of String.",
                            msource.File(),
                            msource.Line(),
                            msource.Column() ); 

  while (!stop) 
    {
      char dummy = msource.Current();
      if ( dummy != '"' )
        {
          res += dummy;
          msource.AdvanceChar();

          if ( msource.LineBreak() )
            throw IllegalLinebreak( "Line break before end of String.",
                                    msource.File(),
                                    msource.Line(),
                                    msource.Column() ); 
        }
      else
        {
          stop = true;
          msource.AdvanceChar();
        }
    }
}


/** Reads an integer. Whitespace has to
    have been eaten before, that is, the current character must be
    a number or `+' or `-'.
  
    Whitespace or line breaks terminate the scanning!
    There are no whitespaces allowed anywhere, consisten with ANSI C
    scanf. 

    \param[out] res Output string containing the integer.

    \exception IllegalLinebreak An illegal linebreak has occured. 
    \exception UnexpectedChar Unexpected character encountered. */
void ArtsParser::read_integer(String& res)
{
  bool stop = false;
  res = "";
  char dummy;
  msource.LineBreak() = false;

  dummy = msource.Current();
  if ( '+' == dummy || '-' == dummy )
    {
      res += dummy;
      msource.AdvanceChar();
      if ( msource.LineBreak() )
        throw IllegalLinebreak( "Line break after sign.",
                                msource.File(),
                                msource.Line(),
                                msource.Column() ); 
    }

  if (!isdigit(msource.Current()))
    {
      ostringstream os;
      os << "Expected digit, but got `" << msource.Current() << "'.";
      throw UnexpectedChar(os.str(),
                           msource.File(),
                           msource.Line(),
                           msource.Column());
    }

  while (!stop) 
    {
      char chtmp = msource.Current();
      if ( isdigit(chtmp) )
        {
          res += chtmp;
          msource.AdvanceChar();
          if ( msource.LineBreak() ) stop = true;
        }
      else
        {
          stop = true;
        }
    }
}


/** Reads a floating point number. Whitespace has to
    have been eaten before, that is, the current character must be
    a number or `+' or `-'.
  
    Example numbers: 23.,  1.0,  -.3,  -3.3e5,  +3e8,  1.0E-9
  
    Illegal numbers: ., 3e, e3, 2e-
  
    Whitespace is not allowed inside the number.
    Line breaks or whitespace terminates the scanning. 

    \param[out] res Output string containing the numeric.

    \exception IllegalLinebreak Illegal line break.
    \exception ParseError Cannot parse this as a number. */
void ArtsParser::read_numeric(String& res)
{
  bool stop;
  res = "";
  char dummy;
  msource.LineBreak() = false;

  // To make sure that there is at least one digit:
  bool found_digit = false;

  // Check if there is a sign:
  dummy = msource.Current();
  if ( '+' == dummy || '-' == dummy )
    {
      res += dummy;
      msource.AdvanceChar();
      if ( msource.LineBreak() )
        throw IllegalLinebreak( "Linebreak after sign.",
                                msource.File(),
                                msource.Line(),
                                msource.Column() ); 
    }

  // There could be some digits here:
  stop = false;
  while (!stop) 
    {
      char chtmp = msource.Current();
      if ( isdigit(chtmp) )
        {
          found_digit = true;
          res += chtmp;
          msource.AdvanceChar();
          if ( msource.LineBreak() ) return; // Line break ends scanning immediately.
        }
      else
        {
          stop = true;
        }
    }

  // Next there can be a decimal point
  if ( '.' == msource.Current() )
    {
      res += ".";
      msource.AdvanceChar();
      if ( msource.LineBreak() )
        {
          if (found_digit)
            {
              // Line break ends scanning immediately, if we have
              // already found at least one digit.
              return;
            }
          else
            {
              throw IllegalLinebreak("Expected at least one digit.",
                                     msource.File(),
                                     msource.Line(),
                                     msource.Column());
            }
        }

      // ... followed by optional more digits
      stop = false;
      while (!stop) 
        {
          char chtmp = msource.Current();
          if ( isdigit(chtmp) )
            {
              found_digit = true;
              res += chtmp;
              msource.AdvanceChar();
              if ( msource.LineBreak() ) return; // Line break ends scanning immediately.
            }
          else
            {
              stop = true;
            }
        }    
    }

  // At this point, we must have found at least one digit.
  if (!found_digit)
    throw ParseError("Expected at least one digit.",
                     msource.File(),
                     msource.Line(),
                     msource.Column());

  // Now there could be a `e' or `E':
  dummy = msource.Current();
  if ( 'e' == dummy || 'E' == dummy )
    {
      res += dummy;
      msource.AdvanceChar();
      if ( msource.LineBreak() )
        throw IllegalLinebreak( "Linebreak after e/E.",
                                msource.File(),
                                msource.Line(),
                                msource.Column() );
      
      // Now there must be an integer (with optional sign)
      {
        String s;
        read_integer(s);
        res += s;
      }
    }
}


/** Use a String stream to parse an integer number.
 
  \param[out] n Parsed integer.
*/
void ArtsParser::parse_integer(Index& n)
{
  String res;
  read_integer(res);
  istringstream is(res);
  is >> n;
}


/** Use a String stream to parse a floating point number.
 
  \param[out] n Parsed numeric.
 */
void ArtsParser::parse_numeric(Numeric& n)
{
  String res;
  read_numeric(res);
  istringstream is(res);
  is >> n;
}


/** Read a vector of Strings. This looks as follows in the control
    file: ["String1","String2"]

    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
  
    Quotation marks inside Strings are currently not possible.
  
    Line breaks are allowed before and after each String. Line breaks
    inside Strings are not allowed. 
   
    \see parse_String */
void ArtsParser::parse_Stringvector(ArrayOfString& res)
{
  bool first = true;            // To skip the first comma.
  res.resize(0);                // Clear the result vector (just in case).

  // Make sure that the current character really is `[' and proceed.
  assertain_character('[');
  // There might have occured a linebreak, which is fine.
  
  eat_whitespace();

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( ']' != msource.Current() )
    {
      String dummy;

      if (first)
        first = false;
      else
        {
          assertain_character(',');
          eat_whitespace();
        }

      parse_String(dummy);
      res.push_back(dummy);
      eat_whitespace();
    }

  msource.AdvanceChar();
}


/** Read a vector of integers. This looks as follows in the control
    file: [123,5,334]
    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
    
    Line breaks are allowed before and after each number. Line breaks
    inside numbers are not allowed. 
   
    \see parse_integer */
void ArtsParser::parse_intvector(ArrayOfIndex& res)
{
  bool first = true;            // To skip the first comma.
  res.resize(0);                        // Clear the result vector (just in case).

  // Make sure that the current character really is `[' and proceed.
  assertain_character('[');
  // There might have occured a linebreak, which is fine.
  
  eat_whitespace();

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( ']' != msource.Current() )
    {
      Index dummy;

      if (first)
        first = false;
      else
        {
          assertain_character(',');
          eat_whitespace();
        }

      parse_integer(dummy);
      res.push_back(dummy);
      eat_whitespace();
    }

  msource.AdvanceChar();
}


/** Read a vector of Numerics. This looks as follows in the control
    file: [1.3, 5, 3.4]
    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
  
    Line breaks are allowed before and after each number. Line breaks
    inside numbers are not allowed. 
   
    \see parse_numeric */
void ArtsParser::parse_numvector(Vector& res)
{
  bool first = true;            // To skip the first comma.

  // We need a temporary Array<Numeric>, so that we can use push_back
  // to store the values. FIXME: Need also constructor for Vector from
  // Array<Numeric>.
  Array<Numeric> tres;

  // Make sure that the current character really is `[' and proceed.
  assertain_character('[');
  // There might have occured a linebreak, which is fine.
  
  eat_whitespace();

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( ']' != msource.Current() )
    {
      Numeric dummy;

      if (first)
        first = false;
      else
        {
          assertain_character(',');
          eat_whitespace();
        }

      parse_numeric(dummy);
      tres.push_back(dummy);
      eat_whitespace();
    }

  // Copy tres to res:
  res.resize(tres.nelem());
  for (int i = 0; i < tres.nelem (); i++)
    {
      res[i] = tres[i];
    }

  msource.AdvanceChar();
}


/** Read a vector of Numerics from a String. This looks as follows: [1.3, 5]
    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
  
    \see parse_numeric */
bool ArtsParser::parse_numvector_from_string (Vector& res, String& str)
{
  bool first = true;            // To skip the first comma.
  size_t pos = 0;

  // We need a temporary Array<Numeric>, so that we can use push_back
  // to store the values.
  Array<Numeric> tres;

  eat_whitespace_from_string (str, pos);

  // Make sure that the current character really is `[' and proceed.
  if (str[pos] != '[')
    {
      throw runtime_error ("No opening bracket\n");
    }

  pos++;

  eat_whitespace_from_string (str, pos);

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( pos < str.length() && str[pos] != ']'  )
    {
      if (first)
        first = false;
      else
        {
          if (str[pos] != ',')
            {
              return false;
            }
          pos++;
          eat_whitespace_from_string (str, pos);
        }

      Numeric dummy;
      istringstream is (str.substr(pos));
      is >> dummy;
      if (is.bad () || is.fail ())
        return false;
      tres.push_back(dummy);
      while (pos < str.length()
             && (isdigit(str[pos]) || str[pos] == '-' || str[pos] == '.'
                 || str[pos] == 'e'))
        pos++;
      eat_whitespace_from_string (str, pos);
    }

  // Copy tres to res:
  res.resize(tres.nelem());
  for (int i = 0; i < tres.nelem (); i++)
    {
      res[i] = tres[i];
    }

  return true;
}


/** Read an of Strings from a String. This looks as follows:
    [ "String1", "String2"]
    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
  
    \see parse_numeric */
bool ArtsParser::parse_stringarray_from_string (ArrayOfString& res, String& str)
{
  bool first = true;            // To skip the first comma.
  size_t pos = 0;

  // We need a temporary Array<Numeric>, so that we can use push_back
  // to store the values.
  ArrayOfString tres;

  eat_whitespace_from_string (str, pos);

  // Make sure that the current character really is `[' and proceed.
  if (str[pos] != '[')
    {
      throw runtime_error ("No opening bracket\n");
    }

  pos++;

  eat_whitespace_from_string (str, pos);

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( pos < str.length() && str[pos] != ']'  )
    {
      if (first)
        first = false;
      else
        {
          if (str[pos] != ',')
            {
              return false;
            }
          pos++;
          eat_whitespace_from_string (str, pos);
        }

      if (str[pos] != '"')
        {
          throw runtime_error ("Expected quotes\n");
        }

      pos++;

      String dummy;
      while ( pos < str.length() && str[pos] != '"'  )
        {
          dummy += str[pos];
          pos++;
        }

      if (pos == str.length() || str[pos] != '"')
        return false;

      tres.push_back(dummy);

      eat_whitespace_from_string (str, pos);
    }

  // Copy tres to res:
  res.resize(tres.nelem());
  for (int i = 0; i < tres.nelem (); i++)
    {
      res[i] = tres[i];
    }

  return true;
}

