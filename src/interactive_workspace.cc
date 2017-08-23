#include "interactive_workspace.h"
#include "global_data.h"
#include "auto_workspace.h"

using global_data::md_data;

extern Verbosity verbosity_at_launch;
extern void (*getaways[])(Workspace&, const MRecord&);
extern WorkspaceMemoryHandler wsmh;

Index get_wsv_id(const char*);

std::string *error_buffer;

size_t InteractiveWorkspace::n_anonymous_variables_ = 0;

InteractiveWorkspace::InteractiveWorkspace() : Workspace(), verbosity_(1, 1, 1)
{
    Workspace::initialize();
    verbosity_at_launch = verbosity_;
}

void InteractiveWorkspace::initialize() {
    define_wsv_group_names();
    Workspace::define_wsv_data();
    Workspace::define_wsv_map();
    define_md_data_raw();
    expand_md_data_raw_to_md_data();
    define_md_map();
    define_md_raw_map();
    define_agenda_data();
    define_agenda_map();
    assert( check_agenda_data() );
    define_species_data();
    define_species_map();
    define_lineshape_data();
    define_lineshape_norm_data();
}

const char * InteractiveWorkspace::execute_agenda(const Agenda *a)
{
    resize();
    try {
        a->execute(*this);
    } catch(const std::runtime_error &e) {
        *error_buffer = e.what();
        return error_buffer->c_str();
    }
    return nullptr;
}

const char * InteractiveWorkspace::execute_workspace_method(long id,
                                                            const long * args_out,
                                                            const long * args_in)
{
    // Make sure verbosity is set.
    Index wsv_id_verbosity = get_wsv_id("verbosity");
    Verbosity& verbosity = *((Verbosity*)this->operator[](wsv_id_verbosity));
    verbosity.set_main_agenda(true);

    CREATE_OUTS;

    const MdRecord &m = md_data[id];
    ArrayOfIndex output{};
    ArrayOfIndex input{};
    TokVal t{};
    Agenda a{};

    if (m.SetMethod()) {
        swap(args_out[0], args_in[0]);
        return nullptr;
    }

    size_t arg_index = 0;
    for (size_t i = 0; i < m.Out().size(); ++i) {
        output.push_back(args_out[arg_index]);
        ++arg_index;
    }
    try {
        for (size_t i = 0; i < m.GOut().size(); ++i) {
            output.push_back(args_out[arg_index]);
            ++arg_index;
        }
    } catch (...) {}

    arg_index = 0;
    for (size_t i = 0; i < m.In().size(); ++i) {
        input.push_back(args_in[arg_index]);
        ++arg_index;
    }
    for (size_t i = 0; i < m.GIn().size(); ++i) {
        input.push_back(args_in[arg_index]);
        ++arg_index;
    }

    try {
        MRecord mr(id, output, input, t, a);
        getaways[id](*this, mr);
    } catch (const std::runtime_error &e) {
        *error_buffer = e.what();
        return error_buffer->c_str();
    }
    return nullptr;
}

void InteractiveWorkspace::set_index_variable(Index id, const Index &src)
{
    *reinterpret_cast<Index*>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_numeric_variable(Index id, const Numeric &src)
{
    *reinterpret_cast<Numeric*>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_string_variable(Index id, const char *src)
{
    *reinterpret_cast<String*>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_array_of_string_variable(Index id,
                                                        size_t n,
                                                        const char * const *src)
{
    ArrayOfString *dst = reinterpret_cast<ArrayOfString*>(this->operator[](id));
    dst->resize(n);
    for (size_t i = 0; i < n; ++i) {
        dst->operator[](i) = String(src[i]);
    }
}

void InteractiveWorkspace::set_array_of_index_variable(Index id,
                                                       size_t n,
                                                       const Index *src)
{
    ArrayOfIndex *dst = reinterpret_cast<ArrayOfIndex*>(this->operator[](id));
    dst->resize(n);
    for (size_t i = 0; i < n; ++i) {
        dst->operator[](i) = src[i];
    }
}

void InteractiveWorkspace::set_vector_variable(Index id,
                                               size_t n,
                                               const Numeric *src)
{
    Vector *dst = reinterpret_cast<Vector*>(this->operator[](id));
    dst->resize(n);
    for (size_t i = 0; i < n; ++i) {
        dst->operator[](i) = src[i];
    }
}

void InteractiveWorkspace::set_matrix_variable(Index id,
                                               size_t m,
                                               size_t n,
                                               const Numeric *src)
{
    Matrix *dst = reinterpret_cast<Matrix*>(this->operator[](id));
    dst->resize(m, n);
    for (size_t i = 0; i < n * m; ++i) {
        dst->get_c_array()[i] = src[i];
    }
}

void InteractiveWorkspace::resize()
{
    Array<stack<WsvStruct *>> ws_new(wsv_data.nelem());
    std::copy(ws.begin(), ws.end(), ws_new.begin());
    //std::copy(ws.end() - n_anonymous_variables_, ws.end(), ws_new.begin() + wsv_data.nelem());
    std::swap(ws, ws_new);
}

Index InteractiveWorkspace::add_variable(Index group_id)
{
    Index id = static_cast<Index>(ws.size());

    ws.push_back(stack<WsvStruct *>());
    push(ws.size()-1, nullptr);
    ws.back().top()->wsv = wsmh.allocate (group_id);
    ws.back().top()->auto_allocated = true;
    ws.back().top()->initialized = true;

    std::stringstream stream;
    stream << "anonymous_variable_" << n_anonymous_variables_;
    String s(stream.str());
    wsv_data.push_back(WsvRecord(s.c_str(), "Created by C API.", group_id));
    WsvMap[s] = id;

    ++n_anonymous_variables_;
    return id;
}

void InteractiveWorkspace::erase_variable(Index i, Index group_id)
{
    WsvStruct *wsvs;
    while (ws[i].size())
    {
        wsvs = ws[i].top();
        if (wsvs->auto_allocated && wsvs->wsv)
        {
            wsmh.deallocate (group_id, wsvs->wsv);
        }
        delete (wsvs);
        ws[i].pop ();
    }
    ws.erase(ws.begin() + i);

    WsvMap.erase(wsv_data[i].Name());
    wsv_data.erase(wsv_data.begin() + i);
    --n_anonymous_variables_;
}

void InteractiveWorkspace::swap(Index i, Index j)
{
    if (is_initialized(i) && is_initialized(j)) {
        std::swap(ws[i], ws[j]);
    } else if (!is_initialized(i) && is_initialized(j)) {
        ws[i].push(ws[j].top());
        ws[j].pop();
    } else if (is_initialized(i) && !is_initialized(j)) {
        ws[j].push(ws[i].top());
        ws[i].pop();
    }
}
