#include "agenda_class.h"
#include "arts.h"
#include "arts_api.h"
#include "auto_md.h"
#include "global_data.h"
#include "interactive_workspace.h"
#include "parser.h"
#include "parameters.h"
#include "workspace_ng.h"

using global_data::wsv_group_names;
using global_data::md_data;
extern Parameters parameters;
extern Verbosity verbosity_at_launch;
extern void (*getaways[])(Workspace&, const MRecord&);
Index get_wsv_id(const char*);

using global_data::MdMap;

extern std::string *error_buffer;

//extern "C" {

////////////////////////////////////////////////////////////////////////////
// Internal Helper Functions
////////////////////////////////////////////////////////////////////////////

void copy_output_and_input(ArrayOfIndex &output,
                           ArrayOfIndex &input,
                           unsigned long n_args_out,
                           const long *args_out,
                           unsigned long n_args_in,
                           const long *args_in)
{
    output.reserve(n_args_out);
    for (size_t i = 0; i < n_args_out; ++i) {
        output.push_back(args_out[i]);
    }
    input.reserve(n_args_in);
    for (size_t i = 0; i < n_args_in; ++i) {
        input.push_back(args_in[i]);
    }
}

////////////////////////////////////////////////////////////////////////////
// Setup and Finalization.
////////////////////////////////////////////////////////////////////////////
void set_parameters(const char *includepath,
                    const char *datapath)
{
    if (includepath) {
        parameters.includepath.push_back(includepath);
    }
    if (datapath) {
        parameters.datapath.push_back(datapath);
    }
}

void initialize()
{
    error_buffer = new std::string;
    InteractiveWorkspace::initialize();
}

void finalize()
{
    delete error_buffer;
}

const char * get_error()
{
    return error_buffer->c_str();
}

////////////////////////////////////////////////////////////////////////////
// Parsing and executing agendas.
////////////////////////////////////////////////////////////////////////////
Agenda * parse_agenda(const char *filename)
{
    Agenda *a = new Agenda;
    ArtsParser parser = ArtsParser(*a, filename, verbosity_at_launch);

    try {
        parser.parse_tasklist();
        a->set_name(filename);
        a->set_main_agenda();
    } catch(const std::runtime_error &e) {
        *error_buffer = std::string(e.what());
        return nullptr;
    }
    return a;
}

Agenda * create_agenda(const char *name)
{
    Agenda * ptr = new Agenda;
    ptr->set_name(name);
    return ptr;
}

void agenda_insert_set(InteractiveWorkspace *workspace,
                       Agenda *a,
                       long id,
                       long group_id)
{
    TokVal t;
    // Index
    if (wsv_group_names[group_id] == "Index") {
        t = TokVal(*reinterpret_cast<Index*>(workspace->operator[](id)));
    }
    // ArrayOfIndex
    if (wsv_group_names[group_id] == "ArrayOfIndex") {
        t = TokVal(*reinterpret_cast<ArrayOfIndex*>(workspace->operator[](id)));
    }
    // String
    if (wsv_group_names[group_id] == "String") {
        t = TokVal(*reinterpret_cast<String*>(workspace->operator[](id)));
    }
    // ArrayOfString
    if (wsv_group_names[group_id] == "ArrayOfString") {
        t = TokVal(*reinterpret_cast<ArrayOfString*>(workspace->operator[](id)));
    }
    // ArrayOfIndex
    if (wsv_group_names[group_id] == "ArrayOfIndex") {
        t = TokVal(*reinterpret_cast<ArrayOfString*>(workspace->operator[](id)));
    }
    // Vector
    if (wsv_group_names[group_id] == "Vector") {
        t = TokVal(*reinterpret_cast<Vector*>(workspace->operator[](id)));
    }
    // Matrix
    if (wsv_group_names[group_id] == "Matrix") {
        t = TokVal(*reinterpret_cast<Matrix*>(workspace->operator[](id)));
    }

    std::stringstream s;
    s << wsv_group_names[group_id] << "Set";

    Index m_id = MdMap.at(String(s.str()));
    ArrayOfIndex output(1), input(0);
    output[0] = id;
    MRecord mr(m_id, output, input, t, Agenda{});
    a->push_back(mr);
}

void agenda_add_method(Agenda * a,
                       const long id,
                       unsigned long n_output_args,
                       const long *output_args,
                       unsigned long n_input_args,
                       const long *input_args)
{
    ArrayOfIndex output, input;
    Agenda aa{};
    TokVal t{};
    copy_output_and_input(output, input, n_output_args, output_args, n_input_args, input_args);
    MRecord mr(id, output, input, t, Agenda{});
    a->push_back(mr);
}

void agenda_clear(Agenda * a)
{
    a->operator=(Agenda());
}

const char * execute_agenda(InteractiveWorkspace *workspace, const Agenda *a)
{
    return workspace->execute_agenda(a);
}

void destroy_agenda(Agenda *a)
{
    delete a;
}

////////////////////////////////////////////////////////////////////////////
// Creating Workspaces
////////////////////////////////////////////////////////////////////////////
InteractiveWorkspace* create_workspace()
{
    return new InteractiveWorkspace();
}

void destroy_workspace(InteractiveWorkspace* workspace)
{
    delete workspace;
}

////////////////////////////////////////////////////////////////////////////
// Accessing WSV Group Information
////////////////////////////////////////////////////////////////////////////

unsigned long get_number_of_groups() {
    return wsv_group_names.size();
}

const char * get_group_name(int i)
{
    return wsv_group_names[i].c_str();
}

////////////////////////////////////////////////////////////////////////////
// Accessing and Executing WSMs
////////////////////////////////////////////////////////////////////////////
unsigned long get_number_of_methods()
{
    return md_data.size();
}

MethodStruct get_method(Index i)
{
    MethodStruct m{
        i, md_data[i].Name().c_str(), md_data[i].Description().c_str(),
        // Output
        md_data[i].Out().size(), nullptr,
        // Generic Output
        md_data[i].GOut().size(), nullptr,
        // Input
        md_data[i].In().size(), nullptr,
        // Generic Input
        md_data[i].GIn().size(), nullptr
    };
    if (m.n_out > 0) {
        m.out = &md_data[i].Out()[0];
    }
    if (m.n_g_in > 0) {
        m.g_in_types = &md_data[i].GInType()[0];
    }
    if (m.n_in > 0) {
        m.in = &md_data[i].In()[0];
    }
    if (m.n_g_out > 0) {
        m.g_out_types = &md_data[i].GOutType()[0];
    }
    return m;
}

const char * get_method_in(Index i, Index j)
{
    return md_data[i].GIn()[j].c_str();
}

const char * get_method_g_in(Index i, Index j)
{
    return md_data[i].GIn()[j].c_str();
}

const char * get_method_g_in_default(Index i, Index j)
{
    return md_data[i].GInDefault()[j].c_str();
}

const char * get_method_out(Index i, Index j)
{
    return md_data[i].GIn()[j].c_str();
}

const char * get_method_g_out(Index i, Index j)
{
    return md_data[i].GOut()[j].c_str();
}

const char * execute_workspace_method(InteractiveWorkspace *workspace,
                                        long id,
                                        unsigned long n_args_out,
                                        const long * args_out,
                                        unsigned long n_args_in,
                                        const long * args_in)
{
    ArrayOfIndex output, input;
    copy_output_and_input(output, input, n_args_out, args_out, n_args_in, args_in);
    return workspace->execute_workspace_method(id, output, input);
}

////////////////////////////////////////////////////////////////////////////
// Accessing and Manipulating WSVs
////////////////////////////////////////////////////////////////////////////
long lookup_workspace_variable(const char *s)
{
    auto it = Workspace::WsvMap.find(s);
    if (it == Workspace::WsvMap.end()) {
        return -1;
    }
    return it->second;
}

unsigned long get_number_of_variables()
{
    return Workspace::wsv_data.size();
}

VariableStruct get_variable(Index i)
{
    const WsvRecord& r = Workspace::wsv_data[i];
    return VariableStruct {
        r.Name().c_str(),
        r.Description().c_str(),
        r.Group()
    };
}

VariableValueStruct get_variable_value(InteractiveWorkspace *workspace, Index id, Index group_id)
{
    VariableValueStruct value{nullptr,
                              workspace->is_initialized(id),
                              {0, 0, 0, 0, 0, 0},
                              nullptr,
                              nullptr};
    // Index
    if (wsv_group_names[group_id] == "Index") {
        if (value.initialized) {
            value.ptr = workspace->operator[](id);
        }
    }
    // ArrayOfIndex
    else if (wsv_group_names[group_id] == "ArrayOfIndex") {
        if (value.initialized) {
            ArrayOfIndex *a = reinterpret_cast<ArrayOfIndex*>(workspace->operator[](id));
            value.dimensions[0] = a->size();
            value.ptr = &a->operator[](0);
        }
    }
    // String
    else if (wsv_group_names[group_id] == "String") {
        if (value.initialized) {
            value.ptr = reinterpret_cast<String*>(workspace->operator[](id))->c_str();
        }
    }
    // Numeric
    else if (wsv_group_names[group_id] == "Numeric") {
        if (value.initialized) {
            value.ptr = reinterpret_cast<Numeric*>(workspace->operator[](id));
        }
    }
    // Vector
    else if (wsv_group_names[group_id] == "Vector") {
        if (value.initialized) {
            Vector *v = reinterpret_cast<Vector*>(workspace->operator[](id));
            value.dimensions[0] = v->nelem();
            if (!v->empty()) {
                value.ptr = v->get_c_array();
            }
        }
    }
    // Matrix
    else if (wsv_group_names[group_id] == "Matrix") {
        if (value.initialized) {
            Matrix *m = reinterpret_cast<Matrix*>(workspace->operator[](id));
            value.dimensions[0]   = m->nrows();
            value.dimensions[1]   = m->ncols();
            if (!m->empty()) {
                value.ptr = m->get_c_array();
            }
        }
    }
    // Tensor3
    else if (wsv_group_names[group_id] == "Tensor3") {
        if (value.initialized) {
            Tensor3* t = reinterpret_cast<Tensor3*>(workspace->operator[](id));
            value.dimensions[0]   = t->npages();
            value.dimensions[1]   = t->nrows();
            value.dimensions[2]   = t->ncols();
            if (!t->empty()) {
                value.ptr = t->get_c_array();
            }
        }
    }
    // Tensor4
    else if (wsv_group_names[group_id] == "Tensor4") {
        if (value.initialized) {
            Tensor4* t = reinterpret_cast<Tensor4*>(workspace->operator[](id));
            value.dimensions[0]   = t->nbooks();
            value.dimensions[1]   = t->npages();
            value.dimensions[2]   = t->nrows();
            value.dimensions[3]   = t->ncols();
            if (!t->empty()) {
                value.ptr = t->get_c_array();
            }
        }
    }
    // Tensor5
    else if (wsv_group_names[group_id] == "Tensor5") {
        if (value.initialized) {
            Tensor5 *t = reinterpret_cast<Tensor5*>(workspace->operator[](id));
            value.dimensions[0]   = t->nshelves();
            value.dimensions[1]   = t->nbooks();
            value.dimensions[2]   = t->npages();
            value.dimensions[3]   = t->nrows();
            value.dimensions[4]   = t->ncols();
            if (!t->empty()) {
                value.ptr = t->get_c_array();
            }
        }
    }
    // Tensor6
    else if (wsv_group_names[group_id] == "Tensor6") {
        if (value.initialized) {
            Tensor6 *t = reinterpret_cast<Tensor6*>(workspace->operator[](id));
            value.ptr = t->get_c_array();
            value.dimensions[0] = t->nvitrines();
            value.dimensions[1] = t->nshelves();
            value.dimensions[2] = t->nbooks();
            value.dimensions[3] = t->npages();
            value.dimensions[4] = t->nrows();
            value.dimensions[5] = t->ncols();
        }
    }
    // Sparse
    else if (wsv_group_names[group_id] == "Sparse") {
        if (value.initialized) {
            Sparse *s = reinterpret_cast<Sparse*>(workspace->operator[](id));
            value.dimensions[0] = s->nrows();
            value.dimensions[1] = s->ncols();
            value.dimensions[2] = s->nnz();
            value.ptr = s->get_element_pointer();
            value.inner_ptr = s->get_column_index_pointer();
            value.outer_ptr = s->get_row_start_pointer();
        }
    } else {
        if (value.initialized) {
            value.ptr = workspace->operator[](id);
        }
    }
    return value;
}

void set_variable_value(InteractiveWorkspace *workspace,
                        long id,
                        long group_id,
                        VariableValueStruct value)
{
    // Agenda
    if (wsv_group_names[group_id] == "Agenda") {
        const Agenda *ptr = reinterpret_cast<const Agenda *>(value.ptr);
        workspace->set_agenda_variable(id, *ptr);
    }
    // Index
    if (wsv_group_names[group_id] == "Index") {
        const Index *ptr = reinterpret_cast<const Index *>(value.ptr);
        workspace->set_index_variable(id, *ptr);
    }
    // Numeric
    if (wsv_group_names[group_id] == "Numeric") {
        const Numeric *ptr = reinterpret_cast<const Numeric *>(value.ptr);
        workspace->set_numeric_variable(id, *ptr);
    }
    // String
    else if (wsv_group_names[group_id] == "String") {
        const char * ptr = reinterpret_cast<const char *>(value.ptr);
        workspace->set_string_variable(id, ptr);
    }
    // Array of String
    else if (wsv_group_names[group_id] == "ArrayOfString") {
        const char * const * ptr = reinterpret_cast<const char * const *>(value.ptr);
        workspace->set_array_of_string_variable(id, value.dimensions[0], ptr);
    }
    // Array of Index
    else if (wsv_group_names[group_id] == "ArrayOfIndex") {
        const Index *ptr = reinterpret_cast<const Index *>(value.ptr);
        workspace->set_array_of_index_variable(id, value.dimensions[0], ptr);
    }
    // Vector
    else if (wsv_group_names[group_id] == "Vector") {
        const Numeric * ptr = reinterpret_cast<const Numeric *>(value.ptr);
        workspace->set_vector_variable(id, value.dimensions[0], ptr);
    }
    // Matrix
    else if (wsv_group_names[group_id] == "Matrix") {
        const Numeric * ptr = reinterpret_cast<const Numeric *>(value.ptr);
        workspace->set_matrix_variable(id, value.dimensions[0], value.dimensions[1], ptr);
    }
}

long add_variable(InteractiveWorkspace *workspace, long group_id, const char *name)
{
    return workspace->add_variable(group_id, name);
}

DLL_PUBLIC
void erase_variable(InteractiveWorkspace *workspace, long id, long group_id)
{
    return workspace->erase_variable(id, group_id);
}
//}
