/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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
  \file   agenda.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Mar 14 08:49:33 2002
  
  \brief  Implementation of agendas.
*/

#include "agenda.h"
#include "methods.h"
#include "wsv_aux.h"
#include "messages.h"
#include "auto_wsv.h"

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
  exit(1);
}

//! Execute an agenda.
/*! 
  This executes the methods specified in tasklist on the given
  workspace. It also checks for errors during the method execution and
  stops the program if an error has occured.

  The workspace itself is made visible by an external declaration.
*/
void Agenda::execute() const
{
  // The workspace:
  extern WorkSpace workspace;

  // The method description lookup table:
  extern const Array<MdRecord> md_data;

  // The workspace variable lookup table:
  extern const Array<WsvRecord> wsv_data;
  
  // The array holding the pointers to the getaway functions:
  extern const void (*getaways[])(WorkSpace&, const MRecord&);

  out3 << "\nExecuting methods:\n";

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
	    //	    cout << "v.nelem(): " << v.nelem() << endl;
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
	    const ArrayOfIndex& v(mdd.Output());
	    for (Index s=0; s<v.nelem(); ++s) workspace.set(v[s]);
	  }

	  { // Flag the generic output variables as occupied:
	    const ArrayOfIndex& v(mrr.Output());
	    for (Index s=0; s<v.nelem(); ++s) workspace.set(v[s]);
	  }

	  // Call the getaway function:
	  getaways[mrr.Id()]
	    ( workspace, mrr );

	}
      catch (runtime_error x)
	{
	  out0 << "Run-time error in method: " << mdd.Name() << '\n'
	       << x.what() << '\n';
	  exit(1);
	}
      catch (logic_error x)
	{
	  out0 << "Logic error in method: " << mdd.Name() << '\n'
	       << x.what() << '\n';
	  exit(1);
	}
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
