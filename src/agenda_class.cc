/* Copyright (C) 2002-2008
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

#include <ostream>
#include <algorithm>
#include <iterator>

#include "arts.h"
#include "agenda_class.h"
#include "agenda_record.h" // only for test functions
#include "methods.h"
#include "wsv_aux.h"
#include "messages.h"
#include "workspace_ng.h"
#include "arts_omp.h"
#include "auto_md.h"


/** Print the error message and exit. */
void give_up(const String& message)
{
  out0 << message << '\n';
  arts_exit();
}


//! Appends methods to an agenda
/*!
  This function appends a workspace method to the agenda. It currently only
  supports appending WSMs which have no generic input or output, and which are
  defined as Set methods which take only one keyword.
   
  The keyword value has to be a string, which for no value should be of length
  zero.
   
  \param ws            Workspace reference
  \param methodname    The name of the WSM
  \param keywordvalue  The value of the keyword

  \author Mattias Ekstrom
  \date   2005-01-05
*/
void Agenda::append(Workspace& ws _U_,
                    const String& methodname,
                    const TokVal& keywordvalue)
{
  extern const map<String, Index> MdMap;

  const map<String, Index>::const_iterator i2 = MdMap.find(methodname);
  assert ( i2 != MdMap.end() );
  Index id = i2->second;            
          
  extern const Array<MdRecord> md_data;
  ArrayOfIndex output = md_data[id].Out();

  // Find explicit method id in MdMap.
  ArrayOfIndex input = md_data[id].In();
  Agenda dummy;
  dummy.resize(0);
  
  // Append the method
  dummy.push_back (MRecord(id,output,input,keywordvalue,dummy));
  AgendaAppend(ws, *this, mname, *this, mname, dummy);
}


//! Checks consistency of an agenda.
/*! 
  Checks that the input used by the agenda and the output produced by the
  actual methods corresponds to what is desired in the lookup data.
*/
void Agenda::check(Workspace& ws)
{
    // Make external data visible
  extern const Array<AgRecord>  agenda_data;
  extern const map<String, Index> AgendaMap;

  // First we have to find the lookup information for this agenda. We
  // use AgendaMap for this.

  map<String, Index>::const_iterator mi =
    AgendaMap.find( mname );

  // Find return end() if the string is not found. This means that the
  // lookup data for this agenda is missing!
  assert( mi != AgendaMap.end() );

  const AgRecord& this_data = agenda_data[mi->second];

  // Ok, we have the lookup data now.

  // Check that the output produced by the actual methods in the
  // agenda corresponds to what is desired in the lookup data:
  for ( Index i=0; i<this_data.Out().nelem(); ++i )
    {
      // The WSV for which to check:
      Index this_wsv = this_data.Out()[i];

      if ( !is_output(this_wsv) )
        {
          ostringstream os;
          os << "The agenda " << mname
             << " must generate the output WSV "
             << Workspace::wsv_data[this_wsv].Name() << ",\n"
             << "but it does not. It only generates:\n";
          for ( Index j=0; j<Workspace::wsv_data.nelem(); ++j )
            if ( is_output(j) )
              os << Workspace::wsv_data[j].Name() << "\n";
          throw runtime_error (os.str());
        }
    }

  // Check that the input used by the actual methods in the
  // agenda corresponds to what is desired in the lookup data:
  for ( Index i=0; i<this_data.In().nelem(); ++i )
    {
      // The WSV for which to check:
      Index this_wsv = this_data.In()[i];

      if ( !is_input(ws, this_wsv) )
        {
          ostringstream os;
          os << "The agenda " << mname
             << " must use the input WSV "
             << Workspace::wsv_data[this_wsv].Name() << ",\n"
             << "but it does not. It only uses:\n";
          for ( Index j=0; j<Workspace::wsv_data.nelem(); ++j )
            if ( is_input(ws, j) )
              os << Workspace::wsv_data[j].Name() << "\n";
          throw runtime_error (os.str());
        }
    }

  set_outputs_to_push_and_dup ();

  mchecked = true;
}


//! Execute an agenda.
/*! 
  This executes the methods specified in tasklist on the given
  workspace. It also checks for errors during the method execution and
  stops the program if an error has occured. 
*/
void Agenda::execute(Workspace& ws) const
{
  
  if (!mchecked)
    {
      ostringstream os;
      os << "Agenda *" << mname << "* hasn't been checked for consistency yet."
        << endl
        << "This check is usually done by AgendaSet or AgendaAppend."
        << endl
        << "However, if you have written code that modifies an Agenda directly"
        << endl
        << "(changing its name or altering its method list), it's up to you to"
        << endl
        << "call Agenda::check after your modifications.";
      throw runtime_error(os.str());
    }

  // If (and only if) this agenda is the main agenda, then we set the
  // thread local global variable in_main_agenda to true, otherwise to
  // false. The variable in_main_agenda is declared in messages.h,
  // defined in messages.cc, and initialized in main.cc.
  bool original_in_main_agenda = in_main_agenda;
  in_main_agenda = is_main_agenda();

  // An empty Agenda name indicates that something going wrong here
  assert (mname != "");

  // The method description lookup table:
  extern const Array<MdRecord> md_data;

  // The array holding the pointers to the getaway functions:
  extern void (*getaways[])(Workspace&, const MRecord&);

  {
    ostringstream os;
    os << "Executing " << name() << "\n"
       << "{\n";
    out1 << os.str();
  }

  for (Index i=0; i<mml.nelem(); ++i)
    {
      // Runtime method data for this method:
      const MRecord&  mrr = mml[i];
      // Method data for this method:
      const MdRecord& mdd = md_data[mrr.Id()];
      
      try
        {
          {
            ostringstream os;
            if ((mdd.SetMethod() && mrr.Out().nelem()
                && Workspace::wsv_data[mrr.Out()[0]].Name().substr(0, 5)
                   == "auto_")
                || (mdd.Name() == "Delete" && mrr.In().nelem()
                && Workspace::wsv_data[mrr.In()[0]].Name().substr(0, 5)
                   == "auto_"))
              {
                ostringstream os;
                os << "- " << mdd.Name() << "\n";
                out3 << os.str();
              }
            else
              {
                ostringstream os;
                os << "- " << mdd.Name() << "\n";
                out1 << os.str();
              }
          }
        
          { // Check if all input variables are initialized:
            const ArrayOfIndex& v(mrr.In());
            for (Index s=0; s<v.nelem(); ++s)
              if ((s != v.nelem()-1 || !mdd.SetMethod())
                  && !ws.is_initialized(v[s])  )
                give_up("Method "+mdd.Name()+" needs input variable: "+
                        Workspace::wsv_data[v[s]].Name());
          }

          { // Check if all output variables which are also used as input
            // are initialized
            const ArrayOfIndex& v = mdd.InOut();
            for (Index s=0; s<v.nelem(); ++s)
              if (!ws.is_initialized(mrr.Out()[v[s]]) )
                give_up("Method "+mdd.Name()+" needs input variable: "+
                        Workspace::wsv_data[mrr.Out()[v[s]]].Name());
          }

          // Call the getaway function:
          getaways[mrr.Id()]( ws, mrr );

        }
      catch (runtime_error x)
        {
          out1 << "}\n";

          // We have to restore the original content of
          // in_main_agenda, otherwise no output will be visible in
          // case the exception is caught higher up and execution
          // continues.
          in_main_agenda = original_in_main_agenda;

          ostringstream os;
          os << "Run-time error in method: " << mdd.Name() << '\n'
               << x.what();

          throw runtime_error(os.str());
        }
    }

  out1 << "}\n";

  // Restore the original content of in_main_agenda:
  in_main_agenda = original_in_main_agenda;
}


//! Retrieve indexes of all input and output WSVs
/*!
  Builds arrays of WSM output variables which need to be
  duplicated or pushed on the WSV stack before the agenda
  is executed.
*/
void Agenda::set_outputs_to_push_and_dup ()
{
  extern const Array<MdRecord>  md_data;

  set<Index> inputs;
  set<Index> outputs;
  set<Index> outs2push;
  set<Index> outs2dup;

  for (Array<MRecord>::const_iterator method = mml.begin ();
       method != mml.end (); method++)
    {
      // Collect output WSVs
      const ArrayOfIndex& outs  = md_data[method->Id()].Out();
      const ArrayOfIndex& gouts = method->Out();

      // Put the outputs into a new set to sort them. Otherwise
      // set_intersection and set_difference screw up.
      set<Index> souts;
      souts.insert ( outs.begin (), outs.end ());
      souts.insert ( gouts.begin (), gouts.end ());

      // Collect generic input WSVs
      const ArrayOfIndex& gins = method->In();
      inputs.insert (gins.begin (), gins.end ());

      /* Special case: For the Delete WSM add its input to the list
       * of output variables to force a duplication of those variables.
       * It avoids deleting variables outside the agenda's scope.
       */
      if (md_data[method->Id()].Name() == "Delete")
        {
          souts.insert ( gins.begin (), gins.end ());
        }

      // Collect input WSVs
      const ArrayOfIndex& ins = md_data[method->Id()].In();
      inputs.insert (ins.begin (), ins.end ());

      // Add all outputs of this WSM to global list of outputs
      outputs.insert (souts.begin (), souts.end ());

      // Find out all output WSVs of current WSM which were
      // already used as input. We have to place a copy of them on
      // the WSV stack.
      set_intersection (souts.begin (), souts.end (),
                        inputs.begin (), inputs.end (),
                        insert_iterator< set<Index> >(outs2dup,
                                                      outs2dup.begin ()));

    }

  // Find all outputs which are not in the list of WSVs to duplicate
  set_difference (outputs.begin (), outputs.end (),
                  outs2dup.begin (), outs2dup.end (),
                  insert_iterator< set<Index> >(outs2push,
                                                outs2push.begin ()));

  extern const map<String, Index> AgendaMap;
  extern const Array<AgRecord> agenda_data;

  const AgRecord& agr = agenda_data[AgendaMap.find (name ())->second];
  const ArrayOfIndex& aout = agr.Out();
  const ArrayOfIndex& ain = agr.In();

  // We have to build a new set of agenda input and output because the
  // set_difference function only works properly on sorted input.
  set<Index> saout;
  set<Index> sain;

  saout.insert ( aout.begin (), aout.end ());
  sain.insert ( ain.begin (), ain.end ());

  moutput_push.clear ();
  moutput_dup.clear ();

  // Remove the WSVs which are agenda input from the list of
  // output variables for which we have to create an new
  // entry on the stack. This is already done for agenda inputs.
  set<Index> outs2push_without_agenda_input;
  set_difference (outs2push.begin (), outs2push.end (),
                  sain.begin (), sain.end (),
                  insert_iterator< set<Index> >(outs2push_without_agenda_input,
                                                outs2push_without_agenda_input.begin ()));

  // Same for agenda output variables.
  set_difference (outs2push_without_agenda_input.begin (),
                  outs2push_without_agenda_input.end (),
                  saout.begin (), saout.end (),
                  insert_iterator<ArrayOfIndex>(moutput_push,
                                                moutput_push.begin ()));

  // Remove the WSVs which are agenda input from the list of
  // output variables for which we have to create a duplicate
  // on the stack. This is already done for agenda inputs.
  set<Index> outs2dup_without_agenda_input;
  set_difference (outs2dup.begin (), outs2dup.end (),
                  sain.begin (), sain.end (),
                  insert_iterator< set<Index> >(outs2dup_without_agenda_input,
                                                outs2dup_without_agenda_input.begin ()));

  // Same for agenda output variables.
  set_difference (outs2dup_without_agenda_input.begin (),
                  outs2dup_without_agenda_input.end (),
                  saout.begin (), saout.end (),
                  insert_iterator<ArrayOfIndex>(moutput_dup,
                                                moutput_dup.begin ()));

  // Special case: Variables which are defined in the agenda only
  // as output but are used first as input in one of the WSMs
  // For those the current WSV value must be copied to the agenda
  // input variable
  set<Index> saout_only;

  set_difference (saout.begin (), saout.end (),
                  sain.begin (), sain.end (),
                  insert_iterator< set<Index> >(saout_only,
                                                saout_only.begin ()));

  ArrayOfIndex agenda_only_out_wsm_in;
  set_intersection (outs2dup.begin (), outs2dup.end (),
                    saout_only.begin (), saout_only.end (),
                    insert_iterator<ArrayOfIndex>(agenda_only_out_wsm_in,
                                                  agenda_only_out_wsm_in.begin ()));

  // Special case: Variables which are defined in the agenda only
  // as input but are used as output in one of the WSMs
  // For those the current WSV value must be copied to the agenda
  // input variable
  set<Index> sain_only;

  set_difference (sain.begin (), sain.end (),
                  saout.begin (), saout.end (),
                  insert_iterator< set<Index> >(sain_only,
                                                sain_only.begin ()));

  ArrayOfIndex agenda_only_in_wsm_out;
  set_intersection (outs2push.begin (), outs2push.end (),
                    sain_only.begin (), sain_only.end (),
                    insert_iterator<ArrayOfIndex>(agenda_only_in_wsm_out,
                                                  agenda_only_in_wsm_out.begin ()));

  out3 << "  [Agenda::pushpop]                 : " << name() << "\n";
  out3 << "  [Agenda::pushpop] - # Funcs in Ag : " << mml.nelem () << "\n";
  out3 << "  [Agenda::pushpop] - AgOut         : ";
  PrintWsvNames (out3, aout);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - AgIn          : ";
  PrintWsvNames (out3, ain);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - All WSM output: ";
  PrintWsvNames (out3, outputs);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - WSVs push     : ";
  PrintWsvNames (out3, moutput_push);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - WSVs dup      : ";
  PrintWsvNames (out3, moutput_dup);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - Ag inp dup    : ";
  PrintWsvNames (out3, agenda_only_in_wsm_out);
  out3 << "\n";
  out3 << "  [Agenda::pushpop] - Ag out dup    : ";
  PrintWsvNames (out3, agenda_only_out_wsm_in);
  out3 << "\n";

  if (agenda_only_in_wsm_out.nelem ())
    {
      ostringstream err;
      err << "At least one variable is only defined as input\n"
        << "in agenda " << name () << ", but\n"
        << "used as output in a WSM called by the agenda!!!\n"
        << "This is not allowed.\n"
        << "Variable(s): ";
      PrintWsvNames (err, agenda_only_in_wsm_out);
      throw runtime_error (err.str ());
    }
}

//! Check if given variable is agenda input.
/*! 
  A variable is agenda input if it is an input variable to any of the
  methods making up the agenda. 

  \param[in,out] ws Current Workspace
  \param[in] var The workspace variable to check.

  \return True if var is an input variable of this agenda.
*/
bool Agenda::is_input(Workspace& ws, Index var) const
{
  // Make global method data visible:
  extern const Array<MdRecord>  md_data;
  extern const ArrayOfString wsv_group_names;

  // Make sure that var is the index of a valid method:
  assert( 0<=var );
  assert( var<md_data.nelem() );

  // Determine the index of WsvGroup Agenda
  Index WsvAgendaGroupIndex = 0;
  for (Index i = 0; !WsvAgendaGroupIndex && i < wsv_group_names.nelem (); i++)
    {
      if (wsv_group_names[i] == "Agenda")
        WsvAgendaGroupIndex = i;
    }

  // Loop all methods in this agenda:
  for ( Index i=0; i<nelem(); ++i )
    {
      // Get a handle on this methods runtime data record:
      const MRecord& this_method = mml[i];
      
      // Is var a specific input?
      {
        // Get a handle on the Input list for the current method:
        const ArrayOfIndex& input = md_data[this_method.Id()].In();

        for ( Index j=0; j<input.nelem(); ++j )
          {
            if ( var == input[j] ) return true;
          }
      }

      // Is var a generic input?
      {
        // Get a handle on the Input list:
        const ArrayOfIndex& input = this_method.In();

        for ( Index j=0; j<input.nelem(); ++j )
          {
            if ( var == input[j] ) return true;
          }
      }

      // If a General Input variable of this method (e.g. AgendaExecute)
      // is of type Agenda, check its input recursively for matches
      for ( Index j = 0; j < md_data[this_method.Id ()].GInType().nelem(); j++)
        {
          if (md_data[this_method.Id ()].GInType()[j] == WsvAgendaGroupIndex)
            {
              Agenda *AgendaFromGeneralInput =
                (Agenda *)ws[this_method.In()[j]];

              if ((*AgendaFromGeneralInput).is_input(ws, var))
                {
                  return true;
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
bool Agenda::is_output(Index var) const
{
  // Loop all methods in this agenda:
  for ( Index i=0; i<nelem(); ++i )
    {
      // Get a handle on this methods runtime data record:
      const MRecord& this_method = mml[i];
      
      // Is var a specific output?
      {
        // Make global method data visible:
        extern const Array<MdRecord>  md_data;

        // Get a handle on the Output list for the current method:
        const ArrayOfIndex& output = md_data[this_method.Id()].Out();

        for ( Index j=0; j<output.nelem(); ++j )
          {
            if ( var == output[j] ) return true;
          }
      }

      // Is var a generic output?
      {
        // Get a handle on the Output list:
        const ArrayOfIndex& output = this_method.Out();

        for ( Index j=0; j<output.nelem(); ++j )
          {
            if ( var == output[j] ) return true;
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
void Agenda::set_name(const String& nname)
{
  mname = nname;
  mchecked = false;
}

//! Agenda name.
/*! 
  Returns the private member mname.

  \return The name of this agenda.
*/
String Agenda::name() const
{
  return mname;
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
void Agenda::print( ostream& os,
                    const String& indent ) const
{
  for ( Index i=0; i<mml.nelem(); ++i )
    {
      // Print member methods with 3 characters more indentation:
      mml[i].print(os,indent);
    }
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
ostream& operator<<(ostream& os, const Agenda& a)
{
  // Print agenda as it would apear in a controlfile.
  a.print(os, "");
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
void MRecord::print( ostream& os,
                     const String& indent ) const
{
  extern const Array<MdRecord>  md_data;

  // Get a handle on the right record:
  const MdRecord tmd = md_data[Id()];

  os << indent << tmd.Name();

  // Is this a generic method? -- Then we need round braces.
  if ( 0 != tmd.GOutType().nelem()+tmd.GInType().nelem() )
    {
      // First entry needs no leading comma:
      bool first=true;

      os << '(';

      for (Index i=0; i<Out().nelem(); ++i)
        {
          if (first) first=false;
          else os << ",";

          os << Workspace::wsv_data[Out()[i]];
        }

      for (Index i=0; i<In().nelem(); ++i)
        {
          if (first) first=false;
          else os << ",";

          os << Workspace::wsv_data[In()[i]];
        }

      os << ')';
    }

  os << "{\n";

  if ( 0 != Tasks().nelem() )
    {
      Tasks().print(os,indent+"   ");
    }

  os << indent << "}";
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
ostream& operator<<(ostream& os, const MRecord& a)
{
  a.print(os,"");
  return os;
}

