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

#ifndef  ARTS_API_INCLUDED
#define ARTS_API_INCLUDED

#define DLL_PUBLIC __attribute__ ((visibility ("default")))

#include "interactive_workspace.h"

extern "C" {

    ////////////////////////////////////////////////////////////////////////////
    // Structs Definitions
    ////////////////////////////////////////////////////////////////////////////

    /** @struct VariableStruct
     * Struct to describe a symbolic WSV
     */
    struct VariableStruct {
        /** @var VariableStruct::name
         * Pointer to the c_str of the name of the variable.
         */
        const char* name;
        /** @var VariableStruct::description
         * Pointer to the c_str of the description of the variable.
         */
        const char* description;
        /** @var VariableStruct::group
         * The Index value representing the group this variable belongs to.
         */
        long group;
    };

    /** @struct VariableValueStruct
     * This struct is used to return the value of an instantiation of a symbolic
     * WSV in a given workspace.
     */
    struct VariableValueStruct {
        /** @var VariableValueStruct::ptr
         * Pointer to the data of the variable or NULL, if the variable
         * is uninitialized.
         */
        const void* ptr;
        /** @var VariableValueStruct::initialized
         * Boolean indicating whether the variable is initialized.
         */
        bool        initialized;
        /** @var VariableValueStruct::dimensions
         * Index array holding the dimensions of the variable. For an n-dimensional
         * tensor (Vector, Matrix, Tensor3, ...) the first n entries will be set to
         * the number of elements the variable contains in the this dimension with the
         * number of columns in the last position.
         */
        long        dimensions[6];
        /** @var VariableValueStruct::inner_pointer
         * This field is only used to return sparse matrices. In this case inner pointer
         * will point to the array of column indices.
         */
        const int *inner_ptr;
        /** @var VariableValueStruct::outer_pointer
         * This field is only used to return sparse matrices. In this case outer pointer
         * will point to the array of column indices.
         */
        const int *outer_ptr;

    };

    /** @struct MethodStruct
     * This struct is used to return the a description of a workspace method.
     */
    struct MethodStruct {
        /** @var MethodStruct::id
         * The WSMs index in md_data.
         */
        long id;
        /** @var MethodStruct::name
         * c_str pointer to the name of the method as defined in methods.cc
         */
        const char* name;
        /** @var MethodStruct::description
         * c_str pointer to the description of the method as defined in methods.cc
         */
        const char* description;
        // Output
        /** @var MethodStruct::n_out
         * c_str Number of output variables.
         */
        unsigned long n_out;
        /** @var MethodStruct::out
         * Pointer to the ArrayOfIndex holding the indices of the output WSVs.
         */
        const long * out;
        // Generic Output
        /** @var MethodStruct::n_g_out
         * The number of generic outputs of the methods.
         */
        unsigned long n_g_out;
        // Generic Output
        /** @var MethodStruct::g_out_types
         * Pointer to the ArrayOfIndex holding the types of the generic output arguments.
         */
        const long * g_out_types;
        // Input
        /** @var MethodStruct::n_in
         * Number of the methods non-generic input arguments.
         */
        unsigned long n_in;
        /** @var MethodStruct::in
         * Pointer to the ArrayOfIndex holding the indices of the input arguments.
         */
        const long * in;
        // Generic Input
        /** @var MethodStruct::n_g_in
         * The number of generic input arguments of the method.
         */
        unsigned long n_g_in;
        /** @var MethodStruct::g_in
         * Pointer to the ArrayOfIndex hodling the types of the generic input arguments.
         */
        const long * g_in_types;
    };

    /** @struct VersionStruct
     * This struct is used to return the version of ARTS
     */
    struct VersionStruct {
        /** @var VersionStruct::major
         * long The major version number of ARTS
         */
        long major;
        /** @var VersionStruct::release
         * long The major release number of this ARTS major version.
         */
        long minor;
        /** @var VersionStruct::
         * long The minor release number of this ARTS major version.
         */
        long revision;
        /** @var MethodStruct::description
         * long The revision number of this ARTS major version.
         */
    };

    ////////////////////////////////////////////////////////////////////////////
    // Setup and Finalization.
    ////////////////////////////////////////////////////////////////////////////

    //! Set ARTS runtime parameters
    /**
     * Add paths to the include- and datapath of the ARTS runtime.
     *
     * \param includepath Pointer to cstring to be added to include path. Use NULL
     *        if nothing should be added to the current path.
     * \param datapath Pointer to cstring to be added to data path. Use NULL
     *        if nothing should be added to the current path.
     */
    DLL_PUBLIC
    void set_parameters(const char *includepath, const char *datapath);

    //! Initalize ARTS runtime.
    /**
     * This function must be called before any other function to initialize the
     * ARTS C API otherwise bad things will probably happen.
     */
    DLL_PUBLIC
    void initialize();

    //! Finalize ARTS runtime.
    /**
     * Deletes the error buffer.
     */
    DLL_PUBLIC
    void finalize();

    //! Get most recent error.
    /**
     * \return c_str pointing to the c_str holding the most recent error message.
     */
    DLL_PUBLIC
    const char * get_error();

    ////////////////////////////////////////////////////////////////////////////
    // Parsing and executing agendas.
    ////////////////////////////////////////////////////////////////////////////

    //! Parse Controlfile
    /** Parses a controlfile into an Agenda and returns a handle to
     *  to this agenda, which can be executed on a workspace using
     *  execute_agenda(...).
     *
     *  Note that parsing the same controlfile twice files if this controlfile
     *  creates new WSVs.
     *
     *  \param filename The path to the controlfile relative to the ARTS search path.
     *  \return Pointer to the agenda holding the parsed controlfile or NULL if
     *   an error occurred during parsing.
     */
    DLL_PUBLIC
    Agenda * parse_agenda(const char *filename);

    //! Create Agenda
    /** Creates an empyt agenda.
     *
     *  \param name Name of the agenda
     *  \return Pointer to the newly created agenda.
     */
    DLL_PUBLIC
    Agenda * create_agenda(const char *name);

    //! Add method to agenda
    /**
     * Add a method to the given agenda.
     *  \param a Pointer to agenda object.
     *  \param id Index of the method to add to the agenda.
     *  \param n_args_out Number of non-generic and generic output arguments.
     *  \param args_out Indices of the output variables.
     *  \param n_args_in Number of non-generic and generic input arguments.
     *  \param args_in Indices of the input variables.
     */
    DLL_PUBLIC
    void agenda_add_method(Agenda *a, const Index id,
                           unsigned long n_args_out,
                           const long * args_out,
                           unsigned long n_args_in,
                           const long * args_in);

    DLL_PUBLIC
    //! Insert a set method into an agenda.
    /**
     * This inserts a set method into the given agenda which sets the a workspace variable
     * given by its ID to to the value it currently has in the given workspace.
     *
     *  \param Pointer to the workspace
     *  \param Pointer to the agenda
     *  \param The workspace variables which should be set to the value it currently has in
     *  the given workspace.
     *  \param The group_id of the workspace variable.
     */
    void agenda_insert_set(InteractiveWorkspace *ws,
                           Agenda * a,
                           long id,
                           long group_id);

    //! Clear Agenda
    /**
     *  Resets a given agenda to be emtpy.
     *  \param a Pointer to agenda object.
     */
    DLL_PUBLIC
    void agenda_clear(Agenda *a);

    //! Execute Agenda
    /** Executes the given agenda on a given workspace.
     *
     *  \param workspace Pointer of the InteractiveWorkspace object to execute the agenda on.
     *  \param a Pointer to the agenda to execute.
     *  \return NULL if execution was successfull, otherwise pointer to the c_str holding the
     *  error message from ARTS.
     */
    DLL_PUBLIC
    const char * execute_agenda(InteractiveWorkspace *workspace, const Agenda *a);

    //! Destroy Agenda
    /**
     * \param a Pointer to the agenda.
     */
    DLL_PUBLIC
    void destroy_agenda(Agenda *a);

    ////////////////////////////////////////////////////////////////////////////
    // Creating Workspaces
    ////////////////////////////////////////////////////////////////////////////

    //! Create new workspace.
    /**
     * \return Pointer to a newly created InteractiveWorkspace object.
     */
    DLL_PUBLIC
    InteractiveWorkspace* create_workspace();

    //! Destroy given workspace.
    DLL_PUBLIC
    void destroy_workspace(InteractiveWorkspace* workspace);

    ////////////////////////////////////////////////////////////////////////////
    // Accessing WSV Group Information
    ////////////////////////////////////////////////////////////////////////////
    //! Return number of WSV groups.
    DLL_PUBLIC
    unsigned long get_number_of_groups();

    //! Get pointer to name of given group.
    DLL_PUBLIC
    const char * get_group_name(int i);

    ////////////////////////////////////////////////////////////////////////////
    // Accessing and Executing WSMs
    ////////////////////////////////////////////////////////////////////////////
    //! Return number of WSMs.
    DLL_PUBLIC
    unsigned long get_number_of_methods();

    //! Return MethodStruct describing method with index i.
    DLL_PUBLIC
    MethodStruct get_method(Index i);

    //! Get name of generic input argument.
    /**
     * \param i The index of the WSM
     * \param j The index of the generic input arguments
     * \return Pointer to c_str holding the name of the input argument as defined in
               methods.cc.
     */
    DLL_PUBLIC
    const char * get_method_g_in(Index i, Index j);

    //! Get default value of generic input argument.
    /**
     * \param i The index of the WSM
     * \param j The index of the generic input argument
     * \return Pointer to c_str holding the default value as defined in methods.cc.
     */
    DLL_PUBLIC
    const char * get_method_g_in_default(Index i, Index j);

    //! Get name value of generic output argument.
    /**
     * \param i The index of the WSM
     * \param j The index of the generic output argument
     * \return Pointer to c_str holding the name of the output argument as
     *         defined in methods.cc.
     */
    DLL_PUBLIC
    const char * get_method_g_out(Index i, Index j);

    //! Execute workspace method.
    /**
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
    const char * execute_workspace_method(InteractiveWorkspace *workspace,
                                          long id,
                                          unsigned long n_args_out,
                                          const long * args_out,
                                          unsigned long n_args_in,
                                          const long * args_in);

    ////////////////////////////////////////////////////////////////////////////
    // Accessing and Manipulating WSVs
    ////////////////////////////////////////////////////////////////////////////
    //! Lookup workspace variable by name.
    /**
     * \param s Pointer to the c_str holding the name of the variable.
     * \return Index of the variable or -1 if no variable with the given name
     *         exists.
     */
    DLL_PUBLIC
    long lookup_workspace_variable(const char *s);

    //! Get number of WSVs.
    DLL_PUBLIC
    unsigned long get_number_of_variables();

    //! Get WSV by index.
    DLL_PUBLIC
    VariableStruct get_variable(Index i);

    //! Get value WSV in given workspace.
    /**
     * \param workspace Pointer to a InteractiveWorkspace object.
     * \param id Index of the workspace variable.
     * \param group_id Index of the group the variable belongs to.
     * \return VariableValueStruct providing access to the variable.
     */
    DLL_PUBLIC
    VariableValueStruct get_variable_value(InteractiveWorkspace *workspace,
                                           Index id,
                                           Index group_id);
    //! Sets the value of a WSV in a given workspace.
    /**
     * This method can be used to set the value of an ARTS WSV from external data. Currently
     * supported data types are:
     *
     * - Index
     * - Numeric
     * - String
     * - ArrayOfString
     * - Vector
     * - Matrix
     *
     * For Index and Numeric the data pointer of the VariableValueStruct is assumed to
     * point to the value to set the WSV to. For variables of type String the data pointer
     * is interpreted as pointer to a c_str. For ArrayOfString the data pointer is assumed
     * to point to an array of pointers to c_str. For Vectors and Matrix data is assumed
     * to point an array of double with c-style memory layout and size given by the
     * corresponding values in the dimension field of VariableValueStruct.
     *
     * \param workspace Pointer to a InteractiveWorkspace object.
     * \param id Index of the workspace variable.
     * \param group_id Index of the group the variable belongs to.
     * \return Poiter to null-terminated string containing the error message if setting of
     * variable fails.
     */
    DLL_PUBLIC
    const char * set_variable_value(InteractiveWorkspace *workspace,
                                    long id,
                                    long group_id,
                                    VariableValueStruct value);
    //! Add variable of given type to workspace.
    /**
     * This adds and initializes a variable in the current workspace and also
     * adds a new entry to the global wsv_data array and updates WsvMap.
     *
     * \param workspace The pointer to the InteractiveWorkspace object where to create
     * the object.
     * \param group_id The index of the group of the variable
     * \param name Pointer to a null-terminated string containing the name of the variable.
     *             If nullptr, a name will be auto-generated.
     * \return The index which indentifies the variable in the given workspace.
     */
    DLL_PUBLIC
    long add_variable(InteractiveWorkspace *workspace, long group_id, const char *name);

    //! Erase variable from workspace.
    /**
     * This variable removes a variable from the workspace and the global wsv_data
     * and WsvMap. This is likely to invalidate indices of existing variables so it
     * should be used with care.
     *
     * \param workspace Pointer to the InteractiveWorkspace object owning the variable
     * \param id The index of the variable
     * \param group_id The index of the group of the variable
     */
    DLL_PUBLIC
    void erase_variable(InteractiveWorkspace *workspace, long id, long group_id);

    //! Get ARTS Version
    /**
     * This function returns the ARTS version number as a double precision floating
     * point number.
     *
     * \return The ARTS version.
     */
    DLL_PUBLIC
    VersionStruct get_version();
}

#endif
