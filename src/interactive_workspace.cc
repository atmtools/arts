#include "interactive_workspace.h"
#include "agenda_class.h"
#include "agenda_record.h"
#include "auto_workspace.h"

extern Verbosity verbosity_at_launch;
extern WorkspaceMemoryHandler wsmh;
extern std::string string_buffer;

namespace global_data {
extern map<String, Index> AgendaMap;
extern Array<MdRecord> md_data;
extern map<String, Index> WsvGroupMap;
}
using global_data::md_data;

Index get_wsv_id(const char *);

size_t InteractiveWorkspace::n_anonymous_variables_ = 0;
std::vector<Callback *> InteractiveWorkspace::callbacks_{};

void callback_getaway(Workspace &ws, const MRecord &mr) {
  InteractiveWorkspace &iws = *reinterpret_cast<InteractiveWorkspace *>(&ws);
  iws.execute_callback(static_cast<Index>(mr.SetValue()));
}

MdRecord callback_mr = MdRecord("APICallback",
                                "",
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                ArrayOfString(),
                                false,
                                false,
                                false,
                                false,
                                false);

InteractiveWorkspace::InteractiveWorkspace(const Index verbosity,
                                           const Index agenda_verbosity)
    : Workspace() {
  Workspace::initialize();
  verbosity_at_launch.set_screen_verbosity(verbosity);
  verbosity_at_launch.set_agenda_verbosity(agenda_verbosity);
  // No report file is used for the C interface
  verbosity_at_launch.set_file_verbosity(0);

  // Set names of agendas.
  using global_data::AgendaMap;
  using global_data::WsvGroupMap;
  const Index WsvAgendaGroupIndex = WsvGroupMap.find("Agenda")->second;
  for (const auto & agenda_iterator : AgendaMap) {
      const auto & wsv_iterator = Workspace::WsvMap.find(agenda_iterator.first);
      if (!(wsv_iterator == Workspace::WsvMap.end())
          && (wsv_data[wsv_iterator->second].Group() == WsvAgendaGroupIndex)) {
          Agenda &agenda = *reinterpret_cast<Agenda *>(this->operator[](wsv_iterator->second));
          agenda.set_name(agenda_iterator.first);
      }
  }
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
  assert(check_agenda_data());
  define_species_data();
  define_species_map();

  // Add getaway for callbacks.
  size_t n_methods = md_data.size();
  getaways[n_methods] = &callback_getaway;
  md_data.push_back(callback_mr);
}

const char *InteractiveWorkspace::execute_agenda(const Agenda *a) {
  // Need to check size of stack as agenda definitions may have
  // added variables.
  if (wsv_data.size() != ws.size()) {
      resize();
  }
  try {
    a->execute(*this);
  } catch (const std::exception &e) {
    string_buffer = e.what();
    return string_buffer.c_str();
  }
  return nullptr;
}

const char *InteractiveWorkspace::execute_workspace_method(
    long id, const ArrayOfIndex &output, const ArrayOfIndex &input) {
  const MdRecord &m = md_data[id];

  // Need to check size of stack as agenda definitions may have
  // added variables.
  if (wsv_data.size() != ws.size()) {
    resize();
  }

  // Check if all input variables are initialized.
  for (Index i : input) {
    if (!is_initialized(i)) {
      string_buffer = "Method " + m.Name() + " needs input " +
                      wsv_data[i].Name() + " but it is uninitialized.";
      return string_buffer.c_str();
    }
  }

  // Make sure verbosity is set.
  Index wsv_id_verbosity = get_wsv_id("verbosity");
  Verbosity &verbosity = *((Verbosity *)this->operator[](wsv_id_verbosity));
  verbosity.set_main_agenda(true);

  CREATE_OUTS;

  if (m.SetMethod()) {
    swap(output[0], input[0]);
    return nullptr;
  }

  TokVal t{};
  Agenda a{};
  try {
    MRecord mr(id, output, input, t, a);
    if (mr.isInternal()) {
      out3 << "- " + m.Name() + "\n";
    } else {
      out1 << "- " + m.Name() + "\n";
    }
    getaways[id](*this, mr);
  } catch (const std::exception &e) {
    string_buffer = e.what();
    return string_buffer.c_str();
  }
  return nullptr;
}

void InteractiveWorkspace::set_agenda_variable(Index id, const Agenda &src) {
  using global_data::AgendaMap;
  Agenda &dst = *reinterpret_cast<Agenda *>(this->operator[](id));
  String old_name = dst.name();
  dst = src;
  if (old_name != "") {
    dst.set_name(old_name);
  }
  dst.check(*this, verbosity_at_launch);
}

void InteractiveWorkspace::set_index_variable(Index id, const Index &src) {
  *reinterpret_cast<Index *>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_numeric_variable(Index id, const Numeric &src) {
  *reinterpret_cast<Numeric *>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_string_variable(Index id, const char *src) {
  *reinterpret_cast<String *>(this->operator[](id)) = src;
}

void InteractiveWorkspace::set_array_of_string_variable(
    Index id, size_t n, const char *const *src) {
  ArrayOfString *dst = reinterpret_cast<ArrayOfString *>(this->operator[](id));
  dst->resize(n);
  for (size_t i = 0; i < n; ++i) {
    dst->operator[](i) = String(src[i]);
  }
}

void InteractiveWorkspace::set_array_of_index_variable(Index id,
                                                       size_t n,
                                                       const Index *src) {
  ArrayOfIndex *dst = reinterpret_cast<ArrayOfIndex *>(this->operator[](id));
  dst->resize(n);
  for (size_t i = 0; i < n; ++i) {
    dst->operator[](i) = src[i];
  }
}

void InteractiveWorkspace::set_vector_variable(Index id,
                                               size_t n,
                                               const Numeric *src) {
  Vector *dst = reinterpret_cast<Vector *>(this->operator[](id));
  dst->resize(n);
  for (size_t i = 0; i < n; ++i) {
    dst->operator[](i) = src[i];
  }
}

void InteractiveWorkspace::set_matrix_variable(Index id,
                                               size_t m,
                                               size_t n,
                                               const Numeric *src) {
  Matrix *dst = reinterpret_cast<Matrix *>(this->operator[](id));
  dst->resize(m, n);
  for (size_t i = 0; i < n * m; ++i) {
    dst->get_c_array()[i] = src[i];
  }
}

void InteractiveWorkspace::set_tensor3_variable(
    Index id, size_t l, size_t m, size_t n, const Numeric *src) {
  Tensor3 *dst = reinterpret_cast<Tensor3 *>(this->operator[](id));
  dst->resize(l, m, n);
  for (size_t i = 0; i < l * n * m; ++i) {
    dst->get_c_array()[i] = src[i];
  }
}

void InteractiveWorkspace::set_tensor4_variable(
    Index id, size_t k, size_t l, size_t m, size_t n, const Numeric *src) {
  Tensor4 *dst = reinterpret_cast<Tensor4 *>(this->operator[](id));
  dst->resize(k, l, m, n);
  for (size_t i = 0; i < k * l * m * n; ++i) {
    dst->get_c_array()[i] = src[i];
  }
}

void InteractiveWorkspace::set_tensor5_variable(Index id,
                                                size_t k,
                                                size_t l,
                                                size_t m,
                                                size_t n,
                                                size_t o,
                                                const Numeric *src) {
  Tensor5 *dst = reinterpret_cast<Tensor5 *>(this->operator[](id));
  dst->resize(k, l, m, n, o);
  for (size_t i = 0; i < k * l * m * n * o; ++i) {
    dst->get_c_array()[i] = src[i];
  }
}

void InteractiveWorkspace::set_tensor6_variable(Index id,
                                                size_t k,
                                                size_t l,
                                                size_t m,
                                                size_t n,
                                                size_t o,
                                                size_t p,
                                                const Numeric *src) {
  Tensor6 *dst = reinterpret_cast<Tensor6 *>(this->operator[](id));
  dst->resize(k, l, m, n, o, p);
  for (size_t i = 0; i < k * l * m * n * o * p; ++i) {
    dst->get_c_array()[i] = src[i];
  }
}

void InteractiveWorkspace::set_tensor7_variable(Index id,
                                                size_t k,
                                                size_t l,
                                                size_t m,
                                                size_t n,
                                                size_t o,
                                                size_t p,
                                                size_t q,
                                                const Numeric *src) {
  Tensor7 *dst = reinterpret_cast<Tensor7 *>(this->operator[](id));
  dst->resize(k, l, m, n, o, p, q);
  for (size_t i = 0; i < k * l * m * n * o * p * q; ++i) {
    dst->get_c_array()[i] = src[i];
  }
}

void InteractiveWorkspace::set_sparse_variable(Index id,
                                               Index m,
                                               Index n,
                                               Index nnz,
                                               const Numeric *src,
                                               const int *row_inds,
                                               const int *col_inds) {
  Sparse *dst = reinterpret_cast<Sparse *>(this->operator[](id));
  *dst = Sparse(m, n);

  Vector elements(nnz);
  ArrayOfIndex row_indices(nnz), column_indices(nnz);

  for (size_t i = 0; i < (size_t)nnz; ++i) {
    elements[i] = src[i];
    row_indices[i] = row_inds[i];
    column_indices[i] = col_inds[i];
  }

  dst->insert_elements(nnz, row_indices, column_indices, elements);
}

void InteractiveWorkspace::resize() {
  Array<stack<WsvStruct *>> ws_new(wsv_data.nelem());
  std::copy(ws.begin(), ws.end(), ws_new.begin());
  std::swap(ws, ws_new);
}

void InteractiveWorkspace::initialize_variable(Index i) {
  this->operator[](i);
  this->pop_free(i);
  this->operator[](i);
}

Index InteractiveWorkspace::add_variable(Index group_id, const char *name) {

  if (wsv_data.size() != ws.size()) {
    resize();
  }

  Index id = static_cast<Index>(ws.size());

  ws.push_back(stack<WsvStruct *>());
  push(ws.size() - 1, nullptr);
  ws.back().top()->wsv = wsmh.allocate(group_id);
  ws.back().top()->auto_allocated = true;
  ws.back().top()->initialized = true;

  String s;
  if (name) {
    s = String(name);
  } else {
    std::stringstream stream;
    stream << "anonymous_variable_" << n_anonymous_variables_;
    s = stream.str();
  }

  wsv_data.push_back(WsvRecord(s.c_str(), "Created by C API.", group_id));
  WsvMap[s] = id;

  ++n_anonymous_variables_;
  return id;
}

void InteractiveWorkspace::erase_variable(Index i, Index group_id) {
  WsvStruct *wsvs;
  while (ws[i].size()) {
    wsvs = ws[i].top();
    if (wsvs->auto_allocated && wsvs->wsv) {
      wsmh.deallocate(group_id, wsvs->wsv);
    }
    delete (wsvs);
    ws[i].pop();
  }
  ws.erase(ws.begin() + i);

  WsvMap.erase(wsv_data[i].Name());
  wsv_data.erase(wsv_data.begin() + i);
  --n_anonymous_variables_;
}

void InteractiveWorkspace::swap(Index i, Index j) {
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
