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

  \param workspace The workspace to act on.
*/
void Agenda::execute(WorkSpace& workspace)
{
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

	  // Call the getaway function:
	  getaways[mrr.Id()]
	    ( workspace, mrr );

	  { // Flag the specific output variables as occupied:
	    const ArrayOfIndex& v(mdd.Output());
	    for (Index s=0; s<v.nelem(); ++s) workspace.set(v[s]);
	  }

	  { // Flag the generic output variables as occupied:
	    const ArrayOfIndex& v(mrr.Output());
	    for (Index s=0; s<v.nelem(); ++s) workspace.set(v[s]);
	  }

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

//! Remove contents and set size to 0.
void Agenda::resize(Index n)
{
  mml.resize(n);
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
