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
  \file   agenda.h
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Mar 14 08:49:33 2002
  
  \brief  Declarations for agendas.
*/

#ifndef agenda_h
#define agenda_h

#include "token.h"

// Declare existance of class WorkSpace
class WorkSpace;

// ... and MRecord
class MRecord;

//! The Agenda class.
/*! An agenda is a list of workspace methods (including keyword data)
  to be executed. There are workspace variables of class agenda that
  can contain a list of methods to execute for a particular purpose,
  for example to compute the lineshape in an absorption
  calculation. 
*/
class Agenda {
public:
  void push_back(MRecord n);
  void execute() const;
  void resize(Index n);
  Index nelem() const;
  Agenda& operator=(const Agenda& x);
  bool is_input(Index var) const;
  bool is_output(Index var) const;
private:
  Array<MRecord> mml;	/*!< The actual list of methods to execute. */
};

/** Method runtime data. In contrast to MdRecord, an object of this
    class contains the runtime information for one method: The method
    id and the keyword parameter values. This is all that the engine
    needs to execute the stack of methods.

    An MRecord includes a member magenda, which can contain an entire
    agenda, i.e., a list of other MRecords. 

    @author Stefan Buehler */
class MRecord {
public:
  MRecord(){ /* Nothing to do here. */ }
  MRecord(const Index id,
	  const Array<TokVal>& values,
	  const ArrayOfIndex& output,
	  const ArrayOfIndex& input,
	  const Agenda&       tasks)
    : mid(id),
      mvalues( values ),
      moutput( output ),
      minput(  input  ),
      mtasks( tasks )
  { 
    // Initialization of arrays from other array now works correctly.
  }
  Index                Id()       const { return mid;     }
  const Array<TokVal>& Values()   const { return mvalues; }
  const ArrayOfIndex&  Output()   const { return moutput; }
  const ArrayOfIndex&  Input()    const { return minput;  }
  const Agenda&        Tasks()    const { return mtasks;  }

private:
  /** Method id. */
  Index mid;
  /** List of parameter values (see methods.h for definition of
      TokVal). */
  Array<TokVal> mvalues;
  /** Output workspace variables (for generic methods). */
  ArrayOfIndex moutput;
  /** Input workspace variables (for generic methods). */
  ArrayOfIndex minput;
  /** An agenda, which can be given in the controlfile instead of
      keywords. */
  Agenda mtasks;
};


#endif
