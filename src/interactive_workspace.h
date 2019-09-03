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
/**
 * @file   interactive_workspace.h
 * @author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
 * @date   2017-08-21
 *
 * @brief This file contains the declaration and partly the implementation
 *        of the InteractiveWorkspace class.
 */
#ifndef _ARTS_INTERACTIVE_WORKSPACE_H_
#define _ARTS_INTERACTIVE_WORKSPACE_H_

#include "agenda_class.h"
#include "workspace_ng.h"

extern void (*getaways[])(Workspace &, const MRecord &);

class InteractiveWorkspace;

/** External callbacks.
 *
 * The Callback class holds references to external callback functions, that
 * should be executed within an ARTS agenda. It is essentially a container
 * that holds a function pointer and provides an execute method to
 * execute the callback stored in the container.
 *
 * Callbacks are assumed to be provided in the form of function pointers that
 * take as single argument the pointer to the interactive workspace within
 * which the callback is executed.
 */
class Callback {
 public:
  /** Create callback
   *
   * @param[in] f Function pointer to the function implementing the
   * callback.
   */
  Callback(void (*f)(InteractiveWorkspace *))
      : callback_(f){
            // Nothing to do here.
        };

  Callback(const Callback &) = default;
  Callback(Callback &&) = default;
  Callback &operator=(const Callback &) = default;
  Callback &operator=(Callback &&) = default;
  ~Callback() = default;

  /** Execute callback
   * 
   * @param[in, out]: The workspace on which to provides as argument to
   * the callback
   */
  virtual void execute(InteractiveWorkspace &ws) {(*callback_)(&ws);};

 private:
  void (*callback_)(InteractiveWorkspace *);
};

/** Interactive ARTS workspace
 *
 * The InteractiveWorkspace class extends the ARTS Workspace class with features
 * that are required to allow an interactive simulation session.
 */
class InteractiveWorkspace : private Workspace {
 public:
  InteractiveWorkspace(const Index verbosity = 1,
                       const Index agenda_verbosity = 0);

  using Workspace::is_initialized;
  using Workspace::operator[];

  /** Workspace intialization
   *
   * This calls the required functions to get the workspace up and running,
   * esentially executing what in called in the ARTS main, when ARTS is
   * used in script mode.
   */
  static void initialize();

  /** Execute agenda
   *
   * Executes the given agenda on this workspace.
   *
   * @param[in] a Pointer to the agenda to execute.
   * @return 0 if execution was successful, otherwise c-string pointer
   * to the generated error message.
   */
  const char *execute_agenda(const Agenda *a);

  /** Execute workspace method
   *
   * Executes the WSM with given id on this workspace.
   *
   * @param[in] id Workspace id of the workspace method to execute.
   * @param[in] output Array containing the indices of the output variables
   * in which to store the results.
   * @param[in] input Array containing the indices of the input variables
   * to the function calls.
   * @return 0 if execution was successful, otherwise c-string pointer
   * to the generated error message.
   */
  const char *execute_workspace_method(long id,
                                       const ArrayOfIndex &output,
                                       const ArrayOfIndex &input);

  /** Set agenda variable into workspace.
   *
   * Copies the given agenda object into the Agenda variable with
   * the given id in this workspace.
   * 
   * Note: This triggers checking of the agenda which can result
   * in an exception being thrown from the call.
   *
   * \param[in] id Workspace id of the agenda to set
   * \param[in] agenda  The agenda to copy into the workspace
   */
  void set_agenda_variable(Index id, const Agenda &src);
  /** Deep-copy of Index variable into workspace.
   *
   * \param[in] id Workspace id of the Index variable to set
   * \param[in] src The index value to copy into the workspace.
   */
  void set_index_variable(Index id, const Index &src);
  /** Deep-copy of ArrayOfIndex variable into workspace.
   *
   * \param[in] id Workspace id of the Index variable to set
   * \param[in] src Pointer to the c-style index array to copy.
   */
  void set_array_of_index_variable(Index id, size_t n, const Index *src);
  /** Deep-copy of Numeric variable into workspace.
   *
   * \param[in] id Workspace id of the Numeric variable to set
   * \param[in] src The index value to copy into the workspace.
   */
  void set_numeric_variable(Index id, const Numeric &src);
  /** Deep-copy of String variable into workspace.
   *
   * \param[in] id Workspace id of the String variable to set
   * \param[in] src Pointer to the null-terminated c-tring to copy
   * into the workspace.
   */
  void set_string_variable(Index id, const char *src);
  /** Deep-copy of ArrayOfString variable into workspace.
   *
   * \param[in] id Workspace id of the ArrayOfString variable to set
   * \param[in] src Pointer to the null-terminated c-string to copy
   * into the workspace.
   */
  void set_array_of_string_variable(Index id, size_t n, const char *const *src);
  /** Deep-copy of Vector variable into workspace.
   *
   * \param[in] id Workspace id of the ArrayOfString variable to set
   * \param[in] src Pointer to the c-array of Numeric containing the
   * vector elements.
   */
  void set_vector_variable(Index id, size_t n, const Numeric *src);
  /** Deep-copy of Matrix variable into workspace.
   *
   * \param[in] id Workspace id of the Vector variable to set
   * \param[in] src Pointer to the two-dimensional  c-array of Numeric containing
   * the matrix elements in column-major order.
   */
  void set_matrix_variable(Index id, size_t m, size_t n, const Numeric *src);
  /** Deep-copy of Tensor3 variable into workspace.
   *
   * \param[in] id Workspace id of the Tensor3 variable to set
   * \param[in] src Pointer to the three-dimensional  c-array of Numeric containing
   * the tensor elements.
   */
  void set_tensor3_variable(
      Index id, size_t l, size_t m, size_t n, const Numeric *src);
  /** Deep-copy of Tensor4 variable into workspace.
   *
   * \param[in] id Workspace id of the Tensor4 variable to set
   * \param[in] src Pointer to the four-dimensional  c-array of Numeric containing
   * the tensor elements.
   */
  void set_tensor4_variable(
      Index id, size_t k, size_t l, size_t m, size_t n, const Numeric *src);
  /** Deep-copy of Tensor5 variable into workspace.
   *
   * \param[in] id Workspace id of the Tensor5 variable to set
   * \param[in] src Pointer to the five-dimensional c-array of Numeric containing
   * the tensor elements.
   */
  void set_tensor5_variable(Index id,
                            size_t k,
                            size_t l,
                            size_t m,
                            size_t n,
                            size_t o,
                            const Numeric *src);
  /** Deep-copy of Tensor6 variable into workspace.
   *
   * \param[in] id Workspace id of the Tensor5 variable to set
   * \param[in] src Pointer to the six-dimensional  c-array of Numeric containing
   * the tensor elements.
   */
  void set_tensor6_variable(Index id,
                            size_t k,
                            size_t l,
                            size_t m,
                            size_t n,
                            size_t o,
                            size_t p,
                            const Numeric *src);
  /** Deep-copy of Tensor5 variable into workspace.
   *
   * \param[in] id Workspace id of the Tensor4 variable to set
   * \param[in] src Pointer to the four-dimensional  c-array of Numeric containing
   * the tensor elements.
   */
  void set_tensor7_variable(Index id,
                            size_t k,
                            size_t l,
                            size_t m,
                            size_t n,
                            size_t o,
                            size_t p,
                            size_t q,
                            const Numeric *src);
  /** Deep-copy of Sparse matrix into workspace.
   *
   * Copies a sparse matrix in coordinate format into the workspace.
   *
   * \param[in] id Workspace id of the Sparse variable to set
   * \param[in] m Number of rows of the matrix
   * \param[in] n Number of column of the matrix
   * \param[in] nnz Number of non-zeros of the matrix
   * \param[in] src Pointer to the one-dimensional c-array containing the matrix
   * elements.
   * \param[in] inner_ptr Pointer to the c-array of Index containing the
   * row indices of the elements in src.
   * \param[in] outer_ptr Pointer to the c-array of Index containing the
   * column indices of the elements in src.
   */
  void set_sparse_variable(Index id,
                           Index m,
                           Index n,
                           Index nnz,
                           const Numeric *src,
                           const int *row_indices,
                           const int *column_indices);

  void execute_callback(Index callback_id) {
    callbacks_[callback_id]->execute(*this);
  }

  /** Add callback to workspace
   *
   * So that a callback can be executed from within an agenda it
   * must be stored within the workspace. This is done with this
   * method. The returned index can then be inserted as a WSM call
   * into a agenda, which will then trigger execution of the workspace
   * when ARTS executes the given agenda.
   *
   * @param[in] cb The callback to add to the workspace.
   * @return The WSM index representing the newly added callback.
   */
  static Index add_callback(Callback *cb) {
    Index id = callbacks_.size();
    callbacks_.push_back(cb);
    return id;
  }

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
  
  /** Resize workspace stack
   *
   * This method resizes the array of stacks that manages the
   * runtime memory of ARTS to match the number of WSV variables.
   */
  void resize();

  static size_t n_anonymous_variables_;
  static std::vector<Callback *> callbacks_;
};

#endif  // _ARTS_INTERACTIVE_WORKSPACE_H_
