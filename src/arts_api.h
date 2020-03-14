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

  \brief This file contains all declarations of the ARTS C API.
*/

#ifndef _ARTS_ARTS_API_H_
#define _ARTS_ARTS_API_H_

#define DLL_PUBLIC __attribute__((visibility("default")))

#include "covariance_matrix.h"
#include "interactive_workspace.h"

extern "C" {

////////////////////////////////////////////////////////////////////////////
// Structs Definitions
////////////////////////////////////////////////////////////////////////////

/** Representation of ARTS WSVs */
struct VariableStruct {
  /** Pointer to the c_str of the name of the variable. */
  const char *name;
  /** Pointer to the c_str of the description of the variable.  */
  const char *description;
  /** The Index value representing the group this variable belongs to. */
  long group;
};

/** Representation of ARTS values
 *
 * This struct is used to transfer values of variables from inside an
 * ARTS workspace to an outside application.
 */
struct VariableValueStruct {
  /** Data pointer
   *
   * Pointer to the data of the variable or NULL, if the variable
   * is uninitialized.
   */
  const void *ptr;
  /** Initialization status
   *
   * Bool value indicating whether the variable is initialized.
   */
  bool initialized;
  /** Dimensions of array data
   *
   * Array holding the dimensions of the variables. This is required information
   * for transferring vector, matrix or tensor data. For an n-dimensional
   * tensor (Vector, Matrix, Tensor3, ...) the first n entries will be set to
   * the number of elements the variable contains in the this dimension. The
   * number of columns is stored in the last position.
   */
  long dimensions[7];
  /** Additional array data for sparse matrices
   *
   * This field is only used to return sparse matrices. In this case inner pointer
   * will point to the array of column indices.
   */
  const int *inner_ptr;
  /** Additional array data for sparse matrices
   *
   * This field is only used to return sparse matrices. In this case outer pointer
   * will point to the array of column indices.
   */
  const int *outer_ptr;
};

/** Representation of workspace methods
 *
 * This struct is used to return descriptions of a workspace method.
 */
struct MethodStruct {
  /** The WSMs index in md_data.  */
  long id;
  /** Method name
   *
   * Pointer to the c string containing the name of the method as defined
   * in methods.cc
   */
  const char *name;
  /** Method description
   *
   * Pointer to the c string containing the description of the workspace method
   * as defined in methods.cc
   */
  const char *description;
  /** Number of non-generic output variables 
   *
   * I.e. length of array pointed to by out.
   * */
  unsigned long n_out;
  /** Output variables
   *
   * Pointer to ArrayOfIndex  holding the indices of the output WSVs.
   */
  const long *out;
  /** Number of generic output variables
   *
   * I.e. length of array pointed to be out.
   */
  unsigned long n_g_out;
  /** Generic output types
   *
   * Pointer to the ArrayOfIndex holding the type indices of the generic output
   * variables.
   */
  const long *g_out_types;
  /** Number of non-generic input variables 
   *
   * I.e. length of array pointed to by in.
   * */
  unsigned long n_in;
  /** Input variables
   *
   * Pointer to ArrayOfIndex  holding the indices of the input WSVs.
   */
  const long *in;
  /** Number of generic input variables
   *
   * Holds the length of the array pointed to by g_in_types.
   */
  unsigned long n_g_in;
  /** Generic input types
   *
   * Pointer to the ArrayOfIndex holding the type indices of the generic input
   * variables.
   */
  const long *g_in_types;
};

/** ARTS version */
struct VersionStruct {
  /** Major version number of ARTS */
  long major;
  /** Minor version number of this ARTS major version. */
  long minor;
  /** Revision number */
  long revision;
};


/** Covariance matrix block
 *
 * This struct is used to return the value of a block of a covariance
 * matrix through the API
 */
struct CovarianceMatrixBlockStruct {
  /** Quantity indices
   *
   * The indices of the retrieval quantities that this block corresponds
   * to.
   */
  long indices[2];

  /** Start indices of block
   *
   * Row and column indices of the upper- and left-most elements
   * of this block w.r.t. to the full covariance matrix.
   */
  long position[2];

  /** Block size
   * 
   * Number of rows and columns of this block.
   */
  long dimensions[2];

  /** Element pointer
   *
   * Pointer to the data of the matrix the block consists of. If the
   * block holds a sparse matrix, this pointer will point to the
   * element vector. If the block holds a dense matrix, this pointer
   * will point to the 2D data array in row-major order.
   */
  const void *ptr;

  /** Number of non-zero elements
   *
   * Contains the number of non-zero elements in the block if
   * it is represented by a sparse matrix. 0 otherwise.
   */
  long nnz;

  /** Inner pointer for sparse matrices
   *
   * If the block contains a sparse matrix, this pointer
   * will point to the array of row indices. Otherwise
   * it will be 0.
   */
  const int *inner_ptr;
  /** Outer pointer for sparse matrices
   * 
   * If the block contains a sparse matrix, this pointer
   * will point to the array of row indices. Otherwise
   * it will be 0.
   */
  const int *outer_ptr;
};

////////////////////////////////////////////////////////////////////////////
// Setup and Finalization.
////////////////////////////////////////////////////////////////////////////

/** Add include path.
 *
 * Pushes a given path to the back of the include path.
 *
 * @param[in] path Pointer to cstring containing the path to
 * add to the include path.
 */
DLL_PUBLIC
void include_path_push(const char *path);

/** Remove last include path.
 * 
 * Remove the most recently added include path. Will result
 * in an exception if the include path is empty.
 */
DLL_PUBLIC
void include_path_pop();

/** Add data path.
 *
 * Adds a given path to the ARTS data path.
 *
 * @param[in] path Pointer to cstring containing the path to
 * add to the data path.
 */
DLL_PUBLIC
void data_path_push(const char *path);

/** Remove last data path.
 *
 * Remove the most recently added data path. Will result
 * in an exception if the data path is empty.
 */
DLL_PUBLIC
void data_path_pop();

/** Initalize ARTS runtime.
 *
 * This function must be called before any other function to initialize the
 * ARTS C API otherwise bad things will probably happen.
 */
DLL_PUBLIC
void initialize();

/** Finalize ARTS runtime.
 *
 * Deletes the error buffer.
 */
DLL_PUBLIC
void finalize();

/** Get most recent error.
 *
 * @return Pointer to the c string holding the most recent error message.
 */
DLL_PUBLIC
const char *get_error();

/** Set the ARTS basename.
 *
 * The basename is a global variable that is used in some cases
 * by some WSMs for example to store or read data.
 */
DLL_PUBLIC
void set_basename(const char *name);

////////////////////////////////////////////////////////////////////////////
// Parsing and executing agendas.
////////////////////////////////////////////////////////////////////////////

/** Parse Controlfile
 *
 * Parses a controlfile into an Agenda and returns a handle to
 * to this agenda, which can be executed on a workspace using
 * execute_agenda(...).
 *
 * Note that parsing the same controlfile twice fails if this controlfile
 * creates new WSVs.
 *
 * @param{in] filename The path to the controlfile relative to the ARTS search path.
 * @return Pointer to the agenda holding the parsed controlfile or NULL if
 *   an error occurred during parsing.
 */
DLL_PUBLIC
Agenda *parse_agenda(const char *filename);

/** Create Agenda
 *
 * Creates an empyt agenda.
 *
 *  @param[in] name Name of the agenda
 *  @return Pointer to the newly created agenda.
 */
DLL_PUBLIC
Agenda *create_agenda(const char *name);

/** Add method to agenda
 *
 * Add a method to the given agenda.
 *  @param a Pointer to agenda object.
 *  @param id Index of the method to add to the agenda.
 *  @param n_args_out Number of non-generic and generic output arguments.
 *  @param args_out Indices of the output variables.
 *  @param n_args_in Number of non-generic and generic input arguments.
 *  @param args_in Indices of the input variables.
 */
DLL_PUBLIC
void agenda_add_method(Agenda *a,
                       const Index id,
                       unsigned long n_args_out,
                       const long *args_out,
                       unsigned long n_args_in,
                       const long *args_in);

/** Insert callback into agenda.
 *
 * Inserts callback to the given function object into the agenda.
 */
DLL_PUBLIC
void agenda_insert_callback(Agenda *a, void (*f)(InteractiveWorkspace *));

/** Insert a set method into an agenda.
 *
 * This inserts a set method into the given agenda which sets the a workspace
 * variable given by its ID to to the value it currently has in the given
 * workspace.
 *
 *  @param Pointer to the workspace
 *  @param Pointer to the agenda
 *  @param Pointer to the callback function object
 *  @param Language type of the callback.
 */
DLL_PUBLIC
void agenda_insert_set(InteractiveWorkspace *ws,
                       Agenda *a,
                       long id,
                       long group_id);

/** Append agendas.
 *
 * Appends the methods in agenda src to agenda dst. This can be
 * used to emulate include statements in agenda definitions.
 *
 * @param dst The agenda to append the methods to
 * @param src The agenda whose methods to append to src
 */
DLL_PUBLIC
void agenda_append(Agenda *dst, const Agenda *src);

/** Clear Agenda
 *
 *  Resets a given agenda to be emtpy.
 *  @param a Pointer to agenda object.
 */
DLL_PUBLIC
void agenda_clear(Agenda *a);

/** Execute Agenda
 *
 * Executes the given agenda on a given workspace.
 *
 *  @param workspace Pointer of the InteractiveWorkspace object to execute the agenda on.
 *  @param a Pointer to the agenda to execute.
 *  @return NULL if execution was successfull, otherwise pointer to the c_str holding the
 *  error message from ARTS.
 */
DLL_PUBLIC
const char *execute_agenda(InteractiveWorkspace *workspace, const Agenda *a);

/** Destroy Agenda
 * 
 * @param a Pointer to the agenda.
 */
DLL_PUBLIC
void destroy_agenda(Agenda *a);

////////////////////////////////////////////////////////////////////////////
// Creating Workspaces
////////////////////////////////////////////////////////////////////////////

/** Create new workspace.
 *
 * @return Pointer to a newly created InteractiveWorkspace object.
 */
DLL_PUBLIC
InteractiveWorkspace *create_workspace(const Index verbosity = 1,
                                       const Index agenda_verbosity = 0);

/** Destroy given workspace. */
DLL_PUBLIC
void destroy_workspace(InteractiveWorkspace *workspace);

////////////////////////////////////////////////////////////////////////////
// Accessing WSV Group Information
////////////////////////////////////////////////////////////////////////////

/** Return number of WSV groups. */
DLL_PUBLIC
unsigned long get_number_of_groups();

/** Get pointer to name of given group. */
DLL_PUBLIC
const char *get_group_name(int i);

////////////////////////////////////////////////////////////////////////////
// Accessing and Executing WSMs
////////////////////////////////////////////////////////////////////////////

/** Return number of WSMs.*/
DLL_PUBLIC
unsigned long get_number_of_methods();

/** Return MethodStruct describing method with index i. */
DLL_PUBLIC
MethodStruct get_method(Index i);

/** Get name of generic input argument.
 *	
 * @param[in] i The index of the WSM
 * @param[in] The index of the generic input arguments
 * @return Pointer to c_str holding the name of the input argument as defined in
 * methods.cc.
 */
DLL_PUBLIC
const char *get_method_g_in(Index i, Index j);

/** Get default value of generic input argument.
 *	
 * @param[in] i The index of the WSM
 * @param[in] j The index of the generic input argument
 * @return Pointer to c_str holding the default value as defined in methods.cc.
 */
DLL_PUBLIC
const char *get_method_g_in_default(Index i, Index j);

/**
 * Get string defining missing default parameter.
 */
DLL_PUBLIC
const char *get_g_in_nodef() {
    return NODEF;
}

/** Get name value of generic output argument.
 *
 * @param[in] i The index of the WSM
 * @param[in] j The index of the generic output argument
 * @return Pointer to c_str holding the name of the output argument as
 *         defined in methods.cc.
 */
DLL_PUBLIC
const char *get_method_g_out(Index i, Index j);

/** Execute workspace method.
 *
 * Execute workspace method with given index. Input arguments are passed using pointers
 * to two arrays holding the indices of the input arguments and the output arguments,
 * respectively. The length of the input (output) arguments must of course match the
 * number of non-generic and generic input (output) arguments.
 *
 * \param i The index of the WSM
 * \param args_out Pointer to the array holding the indices of the output WSVs.
 * \param args_in Pointer to the array holding the indices of the input WSVs.
 * \return NULL if execution of workspace method was successful, pointer to c_str
 *         hodling the ARTS error message otherwise.
 */
DLL_PUBLIC
const char *execute_workspace_method(InteractiveWorkspace *workspace,
                                     long id,
                                     unsigned long n_args_out,
                                     const long *args_out,
                                     unsigned long n_args_in,
                                     const long *args_in);

/** Print method documentation.
 *
 * This prints the documentation of the method as it is found for example
 * in the HTML browser to the stream buffer and returns a pointer to
 * it.
 *
 * @param id The id of the method.
 * @return Pointer to c string holding the documentation of the method
 */
DLL_PUBLIC
const char *method_print_doc(long id);

////////////////////////////////////////////////////////////////////////////
// Accessing and Manipulating WSVs
////////////////////////////////////////////////////////////////////////////
/** Lookup workspace variable by name.
 *
 * @param s Pointer to c string holding the name of the variable.
 * @return Index of the variable or -1 if no variable with the given name
 *         exists.
 */
DLL_PUBLIC
long lookup_workspace_variable(const char *s);

/** Number of defined WSVs. */
DLL_PUBLIC
unsigned long get_number_of_variables();

/** Get WSV by index.
 *
 * @param The index of the WSV
 * @return VariableStruct representing the WSV
 */
DLL_PUBLIC
VariableStruct get_variable(Index i);

/** Get value WSV in given workspace.
 *
 * @param workspace Pointer to a InteractiveWorkspace object.
 * @param id Index of the workspace variable.
 * @param group_id Index of the group the variable belongs to.
 * @return VariableValueStruct providing access to the variable.
 */
DLL_PUBLIC
VariableValueStruct get_variable_value(InteractiveWorkspace *workspace,
                                       Index id,
                                       Index group_id);

/** Return block of covariance matrix.
 *
 *
 * @param m Pointer to covariance matrix from which to return the
 *        block.
 * @param block_index Index of the block to return.
 * @param inverse Whether the block should be returned from the list
 *        of inverse blocks or not.
 *
 * @return Returns the block identified by the given pointer
 *         represented by CovarianceMatrixBlockStruct.
 */
DLL_PUBLIC
CovarianceMatrixBlockStruct get_covariance_matrix_block(CovarianceMatrix *m,
                                                        long block_index,
                                                        bool inverse);

/** Sets the value of a WSV in a given workspace.
 *
 * This method can be used to set the value of an ARTS WSV from external data. Currently
 * supported data types are:
 *
 * - Agenda
 * - ArrayOfIndex
 * - ArrayOfString
 * - Index
 * - Matrix
 * - Numeric
 * - Sparse
 * - String
 * - Tensor3, ..., Tensor7
 * - Vector
 *
 * For Index and Numeric the data pointer of the VariableValueStruct is assumed to
 * point to the value to set the WSV to. For variables of type String the data pointer
 * is interpreted as pointer to a c_str. For ArrayOfString the data pointer is assumed
 * to point to an array of pointers to c_str. For Vectors and Matrix data is assumed
 * to point an array of double with c-style memory layout and size given by the
 * corresponding values in the dimension field of VariableValueStruct.
 *
 * @param workspace Pointer to a InteractiveWorkspace object.
 * @param id Index of the workspace variable.
 * @param group_id Index of the group the variable belongs to.
 * @return Poiter to null-terminated string containing the error message if setting of
 * variable fails.
 */
DLL_PUBLIC
const char *set_variable_value(InteractiveWorkspace *workspace,
                               long id,
                               long group_id,
                               VariableValueStruct value);
/** Add variable of given type to workspace.
 *
 * This adds and initializes a variable in the current workspace and also
 * adds a new entry to the global wsv_data array and updates WsvMap.
 *
 * @param workspace The pointer to the InteractiveWorkspace object where to create
 * the object.
 * @param group_id The index of the group of the variable
 * @param name Pointer to a null-terminated string containing the name of the variable.
 *             If nullptr, a name will be auto-generated.
 * @return The index which indentifies the variable in the given workspace.
 */
DLL_PUBLIC
long add_variable(InteractiveWorkspace *workspace,
                  long group_id,
                  const char *name);

/** Erase variable from workspace.
 *
 * This variable removes a variable from the workspace and the global wsv_data
 * and WsvMap. This is likely to invalidate indices of existing variables so it
 * should be used with care.
 *
 * @param workspace Pointer to the InteractiveWorkspace object owning the variable
 * @param id The index of the variable
 * @param group_id The index of the group of the variable
 */
DLL_PUBLIC
void erase_variable(InteractiveWorkspace *workspace, long id, long group_id);

/** Get ARTS Version
 *
 * This function returns the ARTS version number as a double precision floating
 * point number.
 *
 * @return The ARTS version.
 */
DLL_PUBLIC
VersionStruct get_version();
}

#endif // _ARTS_ARTS_API_
