/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \file   agenda_class.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Thu Mar 14 08:49:33 2002
  
  \brief  Declarations for agendas.
*/

#ifndef agenda_class_h
#define agenda_class_h

#include <set>
#include "token.h"
#include "messages.h"

class MRecord;

class Workspace;

//! The Agenda class.
/*! An agenda is a list of workspace methods (including keyword data)
  to be executed. There are workspace variables of class agenda that
  can contain a list of methods to execute for a particular purpose,
  for example to compute the lineshape in an absorption
  calculation. 
*/
class Agenda {
public:

  Agenda() : mname(), mml(), moutput_push(), moutput_dup(), main_agenda(false),
             mchecked(false)
  { /* Nothing to do here */ }

  /*! 
    Copies an agenda.
  */
  Agenda(const Agenda& x) : mname(x.mname),
                            mml(x.mml),
                            moutput_push(x.moutput_push),
                            moutput_dup(x.moutput_dup),
                            main_agenda(x.main_agenda),
                            mchecked(x.mchecked)
  { /* Nothing to do here */ }


  void append(const String& methodname, const TokVal& keywordvalue);
  void check(Workspace& ws, const Verbosity& verbosity);
  void push_back(const MRecord& n);
  void execute(Workspace& ws) const;
  inline void resize(Index n);
  inline Index nelem() const;
  inline Agenda& operator=(const Agenda& x);
  const Array<MRecord>& Methods() const { return mml; }
  bool has_method(const String& methodname) const;
  void set_methods(const Array<MRecord>& ml) { mml = ml; mchecked = false; }
  void set_outputs_to_push_and_dup(const Verbosity& verbosity);
  bool is_input(Workspace& ws, Index var) const;
  bool is_output(Index var) const;
  void set_name(const String& nname);
  String name() const;
  const ArrayOfIndex& get_output2push() const { return moutput_push; }
  const ArrayOfIndex& get_output2dup() const { return moutput_dup; }
  void print(ostream& os, const String& indent) const;
  void set_main_agenda() { main_agenda = true; mchecked = true; }
  bool is_main_agenda() const { return main_agenda; }
  bool checked() const { return mchecked; }

private:
  String         mname; /*!< Agenda name. */
  Array<MRecord> mml;   /*!< The actual list of methods to execute. */

  ArrayOfIndex moutput_push;

  ArrayOfIndex moutput_dup;

  //! Is set to true if this is the main agenda.
  bool main_agenda;

  /** Flag indicating that the agenda was checked for consistency */
  bool mchecked;
};

// Documentation with implementation.
ostream& operator<<(ostream& os, const Agenda& a);


/** Method runtime data. In contrast to MdRecord, an object of this
    class contains the runtime information for one method: The method
    id. This is all that the engine needs to execute the stack of methods.

    An MRecord includes a member magenda, which can contain an entire
    agenda, i.e., a list of other MRecords. 

    @author Stefan Buehler */
class MRecord {
public:
  MRecord() : mid(-1),
              moutput(),
              minput(),
              msetvalue(),
              mtasks(),
              minternal(false)
  { /* Nothing to do here. */ }

  MRecord(const MRecord& x) : mid(x.mid),
                              moutput(x.moutput),
                              minput(x.minput),
                              msetvalue(x.msetvalue),
                              mtasks(x.mtasks),
                              minternal(x.minternal)
  { /* Nothing to do here */ }

  MRecord(const Index         id,
          const ArrayOfIndex& output,
          const ArrayOfIndex& input,
          const TokVal&       setvalue,
          const Agenda&       tasks,
                bool          internal = false) : mid(id),
                                          moutput(output),
                                          minput(input),
                                          msetvalue(setvalue),
                                          mtasks(tasks),
                                          minternal(internal)
  { /* Nothing to do here */ }

  Index               Id()       const { return mid;       }
  const ArrayOfIndex& Out()      const { return moutput;   }
  const ArrayOfIndex& In()       const { return minput;    }
  const TokVal&       SetValue() const { return msetvalue; }
  const Agenda&       Tasks()    const { return mtasks;    }

  //! Indicates the origin of this method.
  /*!
    Returns true if this method originates from a controlfile and false
    if it was added by the engine.
    E.g., for Create and Delete that are added for variables that are
    added internally to handle literals in the controlfile this flag
    will be set to true.
   */
  bool isInternal() const { return minternal; }

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
  MRecord& operator=(const MRecord& x)
  {
    mid = x.mid;

    msetvalue = x.msetvalue;

    moutput.resize(x.moutput.nelem());
    moutput = x.moutput;

    minput.resize(x.minput.nelem());
    minput = x.minput;

    mtasks.resize(x.mtasks.nelem());
    mtasks = x.mtasks;

    return *this;
  }

  //! Get list of generic input only WSVs.
  /*!
    This function returns an array with the indexes of WSVs which are
    only input variables but not output.

    \param[out] ginonly Index array of input only WSVs.

    \author Oliver Lemke
    \date   2008-02-27
  */
  void ginput_only(ArrayOfIndex& ginonly) const
  {
    ginonly = minput;      // Input
    for (ArrayOfIndex::const_iterator j = moutput.begin(); j < moutput.end(); ++j)
      for (ArrayOfIndex::iterator k = ginonly.begin(); k < ginonly.end(); ++k)
        if (*j == *k)
          {
            //              erase_vector_element(vi,k);
            k = ginonly.erase(k) - 1;
            // We need the -1 here, otherwise due to the
            // following increment we would miss the element
            // behind the erased one, which is now at the
            // position of the erased one.
          }
  }

  // Output operator:
  void print(ostream& os, const String& indent) const;

private:
  /** Method id. */
  Index mid;
  /** Output workspace variables. */
  ArrayOfIndex moutput;
  /** Input workspace variables. */
  ArrayOfIndex minput;
  /** Keyword value for Set methods. */
  TokVal msetvalue;
  /** An agenda, which can be given in the controlfile instead of
      keywords. */
  Agenda mtasks;
  /** Flag if this method is called internally by the engine */
  bool minternal;
};

//! Resize the method list.
/*!
  Resizes the agenda's method list to n elements
 */
inline void Agenda::resize(Index n)
{
  mml.resize(n);
}

//! Return the number of agenda elements.
/*!  
  This is needed, so that we can find out the correct size for
  resize, befor we do a copy.

  \return Number of agenda elements.
*/
inline Index Agenda::nelem() const
{
  return mml.nelem();
}


//! Append a new method to end of list.
/*! 
  This is used by the parser to fill up the agenda.

  \param n New method to add.
*/
inline void Agenda::push_back(const MRecord& n)
{
  mml.push_back(n);
  mchecked = false;
}


//! Assignment operator.
/*! 
  Copies an agenda.
*/
inline Agenda& Agenda::operator=(const Agenda& x)
{
  mml = x.mml;
  mname = x.mname;
  moutput_push = x.moutput_push;
  moutput_dup = x.moutput_dup;
  mchecked = x.mchecked;
  return *this;
}

// Documentation is with implementation.
ostream& operator<<(ostream& os, const MRecord& a);

/** An array of Agenda. */
typedef Array<Agenda> ArrayOfAgenda;

#endif

