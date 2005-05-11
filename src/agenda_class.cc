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
  \file   agenda_class.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Mar 14 08:49:33 2002
  
  \brief  Implementation of agendas.
*/

#include "arts.h"
#include "agenda_class.h"
#include "methods.h"
#include "wsv_aux.h"
#include "messages.h"
#include "auto_wsv.h"
#include "workspace_ng.h"

//! Append a new method to end of list.
/*! 
  This is used by the parser to fill up the agenda.

  \param n New method to add.
*/
void Agenda::push_back(MRecord n)
{
  mml.push_back(n);
}

/** Print the error message and exit. */
void give_up(const String& message)
{
  out0 << message << '\n';
  arts_exit ();
}

//! Execute an agenda.
/*! 
  This executes the methods specified in tasklist on the given
  workspace. It also checks for errors during the method execution and
  stops the program if an error has occured.

  The workspace itself is made visible by an external declaration.
*/
void Agenda::execute(bool silent) const
{
  // The workspace:
  extern Workspace workspace;

  // The method description lookup table:
  extern const Array<MdRecord> md_data;

  // The workspace variable lookup table:
  extern const Array<WsvRecord> wsv_data;
  
  // The array holding the pointers to the getaway functions:
  extern const void (*getaways[])(Workspace&, const MRecord&);

    // The messages level. We will manipulate it in this function, if
  // silent execution is desired.
  extern Messages messages;

  // Backup for the original message level:
  Messages messages_original( messages );

  // Manipulate the message level, to allow silent execution:
  if (silent)
    {
      // Level 4 means that the output should remain visible even for
      // silent execution.  
      if ( messages.screen < 4 )
        messages.screen = 0;
      if ( messages.file < 4 )
        messages.file   = 0;
    }

  out1 << "Executing " << name() << "\n"
       << "{\n";

//   for (Index i=0; i<mml.nelem(); ++i)
//     {
//       const MRecord&  mrr = mml[i];
//       cout << "id, input: " << mrr.Id() << ", ";
//       print_vector(mrr.Input());
//       cout << "id, output: " << mrr.Id() << ", ";
//       print_vector(mrr.Output());
//     }

  for (Index i=0; i<mml.nelem(); ++i)
    {
      // Runtime method data for this method:
      const MRecord&  mrr = mml[i];
      // Method data for this method:
      const MdRecord& mdd = md_data[mrr.Id()];
      
      try
        {
          out1 << "- " << mdd.Name() << '\n';
        
          { // Check if all specific input variables are occupied:
            const ArrayOfIndex& v(mdd.Input());
            for (Index s=0; s<v.nelem(); ++s)
              if (!workspace.is_occupied(v[s]))
                give_up("Method "+mdd.Name()+" needs input variable: "+
                        wsv_data[v[s]].Name());
          }

          { // Check if all generic input variables are occupied:
            const ArrayOfIndex& v(mrr.Input());
            //      cout << "v.nelem(): " << v.nelem() << endl;
            for (Index s=0; s<v.nelem(); ++s)
              if (!workspace.is_occupied(v[s]))
                give_up("Generic Method "+mdd.Name()+" needs input variable: "+
                        wsv_data[v[s]].Name());
          }

          // The output is flagged as occupied before the actual
          // function call, so that output of a method can already be used by
          // an agenda executed by the method. This is important for
          // example in the case of absorption, where an agenda is
          // used to compute the lineshape, which needs some input
          // data pased by the calling method, such as line-width, etc.. 

          { // Flag the specific output variables as occupied:
            //OLE: FIXME
          //  const ArrayOfIndex& v(mdd.Output());
          //  for (Index s=0; s<v.nelem(); ++s) workspace.set(v[s]);
          }

          { // Flag the generic output variables as occupied:
            //OLE: FIXME
           // const ArrayOfIndex& v(mrr.Output());
           // for (Index s=0; s<v.nelem(); ++s) workspace.set(v[s]);
          }

          // Call the getaway function:
          getaways[mrr.Id()]
            ( workspace, mrr );

        }
      catch (runtime_error x)
        {
          out0 << "Run-time error in method: " << mdd.Name() << '\n'
               << x.what() << '\n';
          arts_exit ();
        }
      catch (logic_error x)
        {
          out0 << "Logic error in method: " << mdd.Name() << '\n'
               << x.what() << '\n';
          arts_exit ();
        }
    }

  out1 << "}\n";

  // Restore the original message level:
  if (silent)
    {
      messages = messages_original;
    }
}

//! Set size to n.
void Agenda::resize(Index n)
{
  mml.resize(n);
}

//! Return the number of agenda elements.
/*!  
  This is needed, so that we can find out the correct size for
  resize, befor we do a copy.

  \return Number of agenda elements.
*/
Index Agenda::nelem() const
{
  return mml.nelem();
}

//! Assignment operator.
/*! 
  Size of target must have been adjusted before, otherwise an
  assertion fails. 
*/
Agenda& Agenda::operator=(const Agenda& x)
{
//   cout << name() << ": " << mml.nelem() << "\n"
//        << x.name() << ": " << x.mml.nelem() << "\n\n";
  assert( mml.nelem() == x.mml.nelem() );
  mml = x.mml;
  return *this;
}

//! Check if given variable is agenda input.
/*! 
  A variable is agenda input if it is an input variable to any of the
  methods making up the agenda. 

  \param var The workspace variable to check.

  \return True if var is an input variable of this agenda.
*/
bool Agenda::is_input(Index var) const
{
  // Make global method data visible:
  extern const Array<MdRecord>  md_data;

  // Make sure that var is the index of a valid method:
  assert( 0<=var );
  assert( var<md_data.nelem() );

  // Loop all methods in this agenda:
  for ( Index i=0; i<nelem(); ++i )
    {
      // Get a handle on this methods runtime data record:
      const MRecord& this_method = mml[i];
      
      // Is var a specific input?
      {
        // Get a handle on the Input list for the current method:
        const ArrayOfIndex& input = md_data[this_method.Id()].Input();

        for ( Index j=0; j<input.nelem(); ++j )
          {
            if ( var == input[j] ) return true;
          }
      }

      // Is var a generic input?
      {
        // Get a handle on the Input list:
        const ArrayOfIndex& input = this_method.Input();

        for ( Index j=0; j<input.nelem(); ++j )
          {
            if ( var == input[j] ) return true;
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
        const ArrayOfIndex& output = md_data[this_method.Id()].Output();

        for ( Index j=0; j<output.nelem(); ++j )
          {
            if ( var == output[j] ) return true;
          }
      }

      // Is var a generic output?
      {
        // Get a handle on the Output list:
        const ArrayOfIndex& output = this_method.Output();

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

//! Assignment operator for MRecord.
/*! 
  This is necessary, because it is used implicitly if agendas (which
  contain an array of MRecord) are copied. The default assignment
  operator generated by the compiler does not do the right thing!

  This became clear due to a bug when agendas were re-defined in the
  controlfile, which was discoverd by Patrick.

  The problem is that MRecord contains some arrays. The copy semantics
  for Array require the target Array to have the right size. But if we
  overwrite an old MRecord with a new one, we want all arrays to be
  overwritten. We don't care about their old size.

  \param x The other MRecord to assign.

  \return The freshly assigned MRecord.

  \author Stefan Buehler
  \date   2002-12-02
*/
MRecord& MRecord::operator=(const MRecord& x)
{
  mid = x.mid;

  mvalues.resize(x.mvalues.nelem());
  mvalues = x.mvalues;

  moutput.resize(x.moutput.nelem());
  moutput = x.moutput;

  minput.resize(x.minput.nelem());
  minput = x.minput;

  mtasks.resize(x.mtasks.nelem());
  mtasks = x.mtasks;

  return *this;
}

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
  extern const Array<WsvRecord> wsv_data;
  extern const Array<MdRecord>  md_data;

  // Get a handle on the right record:
  const MdRecord tmd = md_data[Id()];

  os << indent << tmd.Name();

  // Is this a generic method? -- Then we need round braces.
  if ( 0 != tmd.GOutput().nelem()+tmd.GInput().nelem() )
    {
      // First entry needs no leading comma:
      bool first=true;

      os << '(';

      for (Index i=0; i<Output().nelem(); ++i)
        {
          if (first) first=false;
          else os << ",";

          os << wsv_data[Output()[i]];
        }

      for (Index i=0; i<Input().nelem(); ++i)
        {
          if (first) first=false;
          else os << ",";

          os << wsv_data[Input()[i]];
        }

      os << ')';
    }

  os << "{\n";

  // Is this an agenda method?
  if ( 0 != Tasks().nelem() )
    {
      // Assert that the keyword list really is empty:
      assert ( 0 == tmd.Keywords().nelem() );

      Tasks().print(os,indent+"   ");
    }
  else
    {
      // We must have a plain method.

      // Print the keywords:

      // Determine the length of the longest keyword:
      Index maxsize = 0;
      for (Index i=0; i<tmd.Keywords().nelem(); ++i)
        if ( tmd.Keywords()[i].nelem() > maxsize )
          maxsize = tmd.Keywords()[i].nelem();

      // The number of actual parameters must match the number of
      // keywords: 
      assert( tmd.Keywords().nelem() == Values().nelem() );

      for (Index i=0; i<tmd.Keywords().nelem(); ++i)
        {
          os << indent << "   " << setw(maxsize)
             << tmd.Keywords()[i] << " = "
             << Values()[i] << "\n";
        }
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
