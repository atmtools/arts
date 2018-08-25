/* Copyright (C) 2017 Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   interactive_workspace.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2017-08-21

  \brief This file contains the declaration and partly the implementation
         of the InteractiveWorkspace class.
*/

#ifndef INTERACTIVE_WORKSPACE_INCLUDED
#define INTERACTIVE_WORKSPACE_INCLUDED

#include "workspace_ng.h"
#include "agenda_class.h"

class InteractiveWorkspace : private Workspace {
public:

    InteractiveWorkspace(const Index verbosity= 1,
                         const Index agenda_verbosity = 0);

    using Workspace::is_initialized;
    using Workspace::operator[];

    static void initialize();

    const char * execute_agenda(const Agenda *a);


    const char * execute_workspace_method(long id,
                                          const ArrayOfIndex &output,
                                          const ArrayOfIndex &input);

    void set_agenda_variable(Index id, const Agenda &src);
    void set_index_variable(Index id, const Index &src);
    void set_numeric_variable(Index id, const Numeric &src);
    void set_string_variable(Index id, const char *src);
    void set_array_of_string_variable(Index id, size_t n, const char * const *src);
    void set_array_of_index_variable(Index id, size_t n, const Index *src);
    void set_vector_variable(Index id, size_t n, const Numeric *src);
    void set_matrix_variable(Index id, size_t m, size_t n, const Numeric *src);
    void set_tensor3_variable(Index id, size_t l, size_t m, size_t n, const Numeric *src);
    void set_tensor4_variable(Index id, size_t k, size_t l, size_t m, size_t n,
                              const Numeric *src);
    void set_tensor5_variable(Index id, size_t k, size_t l, size_t m, size_t n,
                              size_t o, const Numeric *src);
    void set_tensor6_variable(Index id, size_t k, size_t l, size_t m, size_t n,
                              size_t o, size_t p, const Numeric *src);
    void set_tensor7_variable(Index id, size_t k, size_t l, size_t m, size_t n,
                              size_t o, size_t p, size_t q, const Numeric *src);
    void set_sparse_variable(Index id, Index m, Index n, Index nnz,
                             const Numeric *src,
                             const int *inner_ptr,
                             const int *outer_ptr);
    void resize();

    //! Initialize workspace variable.
    /*!
      If unitialized, initializes the workspace variable. If variable
      exists it is reinitialized to be empty.

      \param id Workspace variable index of the variable to (re)initialize.
    */
    void initialize_variable(Index id);

    //! Push a stack for a new variable to the workspace.
    /*!
    Registers a new variable with and adds a new stack to the given workspace.
    If name is nullptr, a unique name is created.

    \param group_id Index of the group of the variable to add to the workspace.
    \param name Char pointer to the name of the variable.
    \return Index of the newly pushed stack
    */
    Index add_variable(Index group_id, const char *name);

    //! Remove a variable stack from the workspace.
    /*!
    Remove the stack given by index id from the workspace. This is used by
    the C API to control the workspace size.

    \param i Index of the stack to erase.
    \param group_id Index of the group of the variable.
    */
    void erase_variable(Index i, Index group_id);

    //! Remove a variable stack from the workspace.
    /*!
    Remove the stack given by index id from the workspace. This is used by
    the C API to control the workspace size.

    \param i Index of the stack to erase.
    */
    void swap(Index i, Index j);

private:

    static size_t     n_anonymous_variables_;
};

#endif // INTERACTIVE_WORKSPACE_INCLUDED
