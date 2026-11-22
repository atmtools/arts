#include "workspace_agenda_class.h"

#include <auto_wsa.h>
#include <auto_wsm.h>
#include <compare.h>
#include <debug.h>
#include <enumsWorkspaceInitialization.h>

#include <algorithm>
#include <exception>
#include <ranges>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "workspace_class.h"
#include "workspace_method_class.h"

void Agenda::add(const Method& method) {
  checked = false;
  method.add_defaults_to_agenda(*this);
  methods.push_back(method);
}

void Agenda::append_methods(const Agenda& other) { methods.append_range(other.methods); }

void Agenda::prepend_methods(const Agenda& other) { methods.insert_range(methods.begin(), other.methods); }

void Agenda::insert_methods(const Agenda& other, const std::string_view before) {
  methods.insert_range(stdr::find_if(methods, Cmp::eq(before), &Method::get_name), other.methods);
}

namespace {
InAndOut agendas_ins_and_outs(const std::vector<Method>& methods) {
  InAndOut out{};
  auto& [ins_first, outs_first, in_then_out] = out;

  for (const Method& method : methods) {
    const auto& ins  = method.get_ins();
    const auto& outs = method.get_outs();

    stdr::copy_if(ins, std::back_inserter(ins_first), [&](const std::string& i) {
      const auto cmp = Cmp::eq(i);
      return not stdr::any_of(ins_first, cmp) and not stdr::any_of(outs_first, cmp) and
             not stdr::any_of(in_then_out, cmp);
    });

    stdr::copy_if(outs, std::back_inserter(in_then_out), [&](const std::string& i) {
      const auto cmp = Cmp::eq(i);
      return stdr::any_of(ins_first, cmp) and not stdr::any_of(in_then_out, cmp);
    });

    stdr::copy_if(outs, std::back_inserter(outs_first), [&](const std::string& i) {
      const auto cmp = Cmp::eq(i);
      return not stdr::any_of(outs_first, cmp) and not stdr::any_of(ins_first, cmp) and
             not stdr::any_of(in_then_out, cmp);
    });
  }

  auto sort_and_erase_copies = [](std::vector<std::string>& vec) {
    stdr::sort(vec);
    auto [s, e] = stdr::unique(vec);
    vec.erase(s, e);
  };
  sort_and_erase_copies(ins_first);
  sort_and_erase_copies(outs_first);
  sort_and_erase_copies(in_then_out);

  return out;
}

struct BadOrRecoverable {
  std::vector<std::string> required, can_be_fixed;
};

BadOrRecoverable missing_inputs(const std::vector<std::string>& must_in,
                                const std::vector<std::string>& ins_first,
                                const std::vector<std::string>& outs_first) {
  BadOrRecoverable out{};

  auto& [required, can_be_fixed] = out;

  for (const std::string& i : must_in) {
    if (stdr::binary_search(outs_first, i)) required.push_back(i);
    if (not stdr::binary_search(ins_first, i)) can_be_fixed.push_back(i);
  }

  return out;
}

BadOrRecoverable missing_outputs(const std::vector<std::string>& must_out,
                                 const std::vector<std::string>& ins_first,
                                 const std::vector<std::string>& outs_first,
                                 const std::vector<std::string>& in_then_out) {
  BadOrRecoverable out{};

  auto& [required, can_be_fixed] = out;

  for (const std::string& o : must_out) {
    if (stdr::binary_search(ins_first, o) and not stdr::binary_search(in_then_out, o)) required.push_back(o);
    if (not stdr::binary_search(outs_first, o) and not stdr::binary_search(in_then_out, o)) can_be_fixed.push_back(o);
  }
  return out;
}
}  // namespace

InAndOut Agenda::get_required_inputs() const { return agendas_ins_and_outs(methods); }

void Agenda::finalize(bool fix) try {
  const static auto& wsa = workspace_agendas();

  auto                            ag_ptr   = wsa.find(name);
  const std::vector<std::string>  empty    = {};
  const std::vector<std::string>& must_out = ag_ptr == wsa.end() ? empty : ag_ptr->second.output;
  const std::vector<std::string>& must_in  = ag_ptr == wsa.end() ? empty : ag_ptr->second.input;

  auto [ins_first, outs_first, in_then_out] = agendas_ins_and_outs(methods);
  const auto [rin, fin]                     = missing_inputs(must_in, ins_first, outs_first);
  const auto [rut, fut]                     = missing_outputs(must_out, ins_first, outs_first, in_then_out);

  std::string error{};
  if (not rin.empty()) error += arts::format("Required input(s) {:B,} are instead output(s)\n", rin);
  if (not rut.empty()) error += arts::format("Required output(s) {:B,} are instead input(s)\n", rut);
  if (not fix) {
    if (not fin.empty()) error += arts::format("Must use {:B,} - set fix to ignore this error\n", fin);
    if (not fut.empty()) error += arts::format("Must set {:B,} - set fix to ignore this error\n", fut);
  }
  if (not error.empty()) throw std::runtime_error(error);

  using T1 = std::vector<std::string>;
  using T2 = std::unordered_map<std::string, std::string>;
  for (const std::string& i : fin) methods.emplace_back("Ignore", T1{i}, T2{});
  for (const std::string& o : fut) methods.emplace_back("Touch", T1{o}, T2{});

  std::erase_if(ins_first, [&must_in](const auto& str) { return stdr::any_of(must_in, Cmp::eq(str)); });

  std::erase_if(ins_first, [&in_then_out](const auto& str) { return stdr::binary_search(in_then_out, str); });

  std::erase_if(in_then_out, [&must_out, &must_in](const auto& str) {
    return stdr::any_of(must_out, Cmp::eq(str)) and stdr::any_of(must_in, Cmp::eq(str));
  });

  copy  = in_then_out;
  share = ins_first;

  checked = true;
} catch (std::exception& e) {
  throw std::runtime_error(arts::format(R"(Cannot finalize agenda "{}"

{}
)",
                                        name,
                                        e.what()));
}

void Agenda::set_name(const std::string& v, bool finalize_fix) {
  name = v;
  finalize(finalize_fix);
}

void Agenda::share_workspace(Workspace& out, const Workspace& in) const try {
  for (auto& str : share) {
    if (not out.contains(str)) out.set(str, in.share(str));
  }

  for (auto& str : copy) {
    if (not out.contains(str)) out.set(str, in.share(str));
  }
} catch (std::exception& e) {
  throw std::runtime_error(arts::format(
      R"(
Cannot share workspace for "{}"

Workspace contains:
{:s}

{})",
      name,
      in,
      e.what()));
}

namespace {
void agenda_add_inner_logic(Workspace& out, const Workspace& in, WorkspaceAgendaBoolHandler handle) {
startover:
  for (auto& var : out) {
    if (var.second.holds<Agenda>()) {
      if (not handle.has(var.first)) {
        auto& ag = var.second.get<Agenda>();
        handle.set(var.first);
        ag.share_workspace(out, in);
        goto startover;
      }
    }
  }
}
}  // namespace

void Agenda::copy_workspace(Workspace& out, const Workspace& in) const try {
  for (auto& str : share) out.set(str, in.share(str));

  //! If copy and share are the same, copy will overwrite share (keep them unique!)
  //! Also if named-agenda call has set output variable, copy will take a copy of variable
  for (auto& str : copy) {
    if (out.contains(str)) {
      out.overwrite(str, out.copy(str));
    } else {
      out.set(str, in.copy(str));
    }
  }

  WorkspaceAgendaBoolHandler handle;
  handle.set(name);
  agenda_add_inner_logic(out, in, handle);
} catch (std::exception& e) {
  throw std::runtime_error(arts::format(
      R"(
Cannot copy workspace for "{}"

Workspace contains:
{:s}

{})",
      name,
      in,
      e.what()));
}

void Agenda::copy_only_workspace(Workspace& out, const Workspace& in) const try {
  for (auto& str : copy) {
    if (out.contains(str)) {
      out.overwrite(str, out.copy(str));
    } else {
      out.set(str, in.copy(str));
    }
  }
} catch (std::exception& e) {
  throw std::runtime_error(arts::format(
      R"(
Cannot copy only workspace in Agenda "{}"

Workspace contains:
{:s}

{})",
      name,
      in,
      e.what()));
}

void Agenda::execute(Workspace& ws) const {
  for (auto& method : methods) {
    try {
      method(ws);
    } catch (std::exception& e) {
      throw std::runtime_error(arts::format(R"(Error executing method "{}":
{})",
                                            method.get_name(),
                                            e.what()));
    }
  }
}

bool Agenda::has_method(const std::string& method) const {
  for (auto& m : methods) {
    if (m.get_name() == method) return true;
  }
  return false;
}

Agenda::Agenda(std::string                     n,
               const std::vector<Method>&      m,
               const std::vector<std::string>& s,
               const std::vector<std::string>& c,
               bool                            check)
    : name(std::move(n)), methods(m), share(s), copy(c), checked(check) {}

std::string Agenda::sphinx_list(const std::string_view prep) const {
  std::string out{};

  for (auto& item : share) { out += arts::format("{}Shares the global *{}*\n", prep, item); }

  for (auto& item : copy) { out += arts::format("{}Copies the global *{}*\n", prep, item); }

  for (auto& method : methods) { out += arts::format("{}{}\n", prep, method.sphinx_list_item()); }

  return out;
}

void Agenda::change_default(const std::string_view name, Wsv value) {
  assert(not name.empty());

  for (Method& method : methods) {
    const auto& method_name = method.get_name();
    if (method_name.starts_with(internal_prefix) and method_name.ends_with(name)) {
      method.change_default(std::move(value));
      return;
    }
  }

  arts::user_error(arts::format(R"(Agenda "{}" does not have a default input "{}")", this->name, name));
}

void xml_io_stream<Agenda>::write(std::ostream& os, const Agenda& x, bofstream* pbofs, std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.get_name(), pbofs);
  xml_write_to_stream(os, x.get_methods(), pbofs);
  xml_write_to_stream(os, x.get_share(), pbofs);
  xml_write_to_stream(os, x.get_copy(), pbofs);
  xml_write_to_stream(os, x.is_checked(), pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Agenda>::read(std::istream& is, Agenda& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  std::string              name{};
  std::vector<Method>      methods;
  std::vector<std::string> share{};
  std::vector<std::string> copy{};
  bool                     checked{false};
  xml_read_from_stream(is, name, pbifs);
  xml_read_from_stream(is, methods, pbifs);
  xml_read_from_stream(is, share, pbifs);
  xml_read_from_stream(is, copy, pbifs);
  xml_read_from_stream(is, checked, pbifs);
  x = Agenda{name, methods, share, copy, checked};

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
