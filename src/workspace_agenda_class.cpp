#include "workspace_agenda_class.h"

#include <auto_wsa.h>
#include <auto_wsm.h>
#include <compare.h>
#include <debug.h>

#include <algorithm>
#include <exception>
#include <ranges>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "enumsWorkspaceInitialization.h"
#include "time_report.h"
#include "workspace_class.h"
#include "workspace_method_class.h"

void Agenda::add(const Method& method) {
  checked = false;
  method.add_defaults_to_agenda(*this);
  methods.push_back(method);
}

void Agenda::finalize(bool fix) try {
  const static auto& wsa = workspace_agendas();

  auto                            ag_ptr   = wsa.find(name);
  const std::vector<std::string>  empty    = {};
  const std::vector<std::string>& must_out = ag_ptr == wsa.end() ? empty : ag_ptr->second.output;
  const std::vector<std::string>& must_in  = ag_ptr == wsa.end() ? empty : ag_ptr->second.input;

  std::vector<std::string> ins_first;
  std::vector<std::string> outs_first;
  std::vector<std::string> in_then_out;

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

  for (const std::string& i : must_in) {
    if (stdr::binary_search(outs_first, i)) {
      throw std::runtime_error(std::format(R"(Agenda "{}" first uses "{}" as an input but it is an output)", name, i));
    }

    if (not stdr::binary_search(ins_first, i)) {
      if (fix) {
        methods.emplace_back("Ignore", std::vector<std::string>{i}, std::unordered_map<std::string, std::string>{});
      } else {
        throw std::runtime_error(std::format(R"(Agenda "{}" does not use "{}")", name, i));
      }
    }
  }

  for (const std::string& o : must_out) {
    if (stdr::binary_search(ins_first, o) and not stdr::binary_search(in_then_out, o)) {
      throw std::runtime_error(std::format(
          R"(Agenda "{}" uses "{}" only as an input but it is Agenda in-out

Agenda user input:      {:B,}
Agenda user output:     {:B,}
Agenda user in-out:     {:B,}

Agenda required input:  {:B,}
Agenda required output: {:B,}
)",
          name,
          o,
          ins_first,
          outs_first,
          in_then_out,
          must_in,
          must_out));
    }

    if (not stdr::binary_search(outs_first, o) and not stdr::binary_search(in_then_out, o)) {
      if (fix) {
        methods.emplace_back("Touch", std::vector<std::string>{o}, std::unordered_map<std::string, std::string>{});
      } else {
        throw std::runtime_error(std::format(R"(Agenda "{}" does not set "{}")", name, o));
      }
    }
  }

  std::erase_if(ins_first, [&must_in](const auto& str) { return stdr::any_of(must_in, Cmp::eq(str)); });

  std::erase_if(ins_first, [&in_then_out](const auto& str) { return stdr::binary_search(in_then_out, str); });

  std::erase_if(in_then_out, [&must_out, &must_in](const auto& str) {
    return stdr::any_of(must_out, Cmp::eq(str)) and stdr::any_of(must_in, Cmp::eq(str));
  });

  copy  = in_then_out;
  share = ins_first;

  checked = true;
} catch (std::exception& e) {
  throw std::runtime_error(std::format(R"(Cannot finalize agenda "{}"

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
  throw std::runtime_error(std::format(
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
  throw std::runtime_error(std::format(
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
  throw std::runtime_error(std::format(
      R"(
Cannot copy only workspace in Agenda "{}"

Workspace contains:
{:s}

{})",
      name,
      in,
      e.what()));
}

void Agenda::execute(Workspace& ws) const try {
  for (auto& method : methods) method(ws);
} catch (std::exception& e) {
  throw std::runtime_error(std::format(R"(Cannot execute "{}"

{})",
                                       name,
                                       e.what()));
}

std::vector<Agenda> Agenda::par_tasks(Workspace& ws) const try {
  ARTS_TIME_REPORT

  auto filter = stdv::filter([](const std::string& s) {
    return not s.starts_with(named_input_prefix) and not s.starts_with(internal_prefix);
  });
  auto concat = stdv::join | filter;

  stdr::for_each(methods | stdv::transform(&Method::get_outs) | concat,
                 [&ws](const std::string& o) { ws.init_if_new(o); });

  std::unordered_set<std::string> all_outputs{};

  std::vector<Agenda>             tasks{};
  std::vector<Method>             task_methods{};
  std::string                     task_name{};
  std::unordered_set<std::string> task_output{};
  std::unordered_set<std::string> task_share{};
  std::vector<Method>             possible_task_methods{};

  auto flush_batch = [&]() {
    if (task_methods.empty()) return;
    tasks.emplace_back(task_name,
                       task_methods,
                       std::vector<std::string>(std::from_range, task_share),
                       std::vector<std::string>{},
                       false);
    task_methods.clear();
    task_share.clear();
    task_name = {};

    all_outputs.insert(task_output.begin(), task_output.end());
    task_output.clear();
  };

  const auto insert_share = [&](const Method& method) {
    for (const auto& in : method.get_ins() | filter) task_share.insert(in);
    for (const auto& out : method.get_outs() | filter) task_share.insert(out);
  };

  const auto insert_out = [&](const Method& method) {
    for (const auto& out : method.get_outs() | filter) task_output.insert(out);
  };

  for (auto& method : methods) {
    if (method.get_setval().has_value()) {
      if (method.is_callback()) {
        task_name = method.get_name();
        task_methods.insert_range(task_methods.end(), possible_task_methods);
        task_methods.push_back(method);
        insert_share(method);
        insert_out(method);
        flush_batch();
      } else {
        possible_task_methods.push_back(method);
      }

      continue;
    }

    if (stdr::any_of(method.get_ins() | filter,
                     [&all_outputs](const std::string& s) { return all_outputs.contains(s); })) {
      throw std::runtime_error(std::format(
          R"(Cannot execute in parallel: method "{}" has a dependency (input) on the output of a previously defined parallel task.)",
          method.get_name()));
    }

    if (stdr::any_of(method.get_outs() | filter,
                     [&all_outputs](const std::string& s) { return all_outputs.contains(s); })) {
      throw std::runtime_error(std::format(
          R"(Cannot execute in parallel: method "{}" has a dependency (output) on the output of a previously defined parallel task.)",
          method.get_name()));
    }

    const bool overlap_current =
        stdr::any_of(method.get_ins() | filter,
                     [&task_output](const std::string& s) { return task_output.contains(s); }) or
        stdr::any_of(method.get_outs() | filter,
                     [&task_output](const std::string& s) { return task_output.contains(s); });

    if (not overlap_current) flush_batch();

    task_methods.insert_range(task_methods.end(), possible_task_methods);
    possible_task_methods.clear();
    task_methods.push_back(method);
    insert_share(method);
    insert_out(method);
    task_name += std::format("{}{}", task_name.empty() ? ""sv : "\n"sv, method.get_name());
  }
  flush_batch();

  return tasks;
} catch (std::exception& e) {
  throw std::runtime_error(std::format(R"(Cannot perform parallelization:

{}

{})",
                                       name,
                                       e.what()));
}

void Agenda::par_execute(Workspace& ws) const try {
  ARTS_TIME_REPORT

  const auto tasks = par_tasks(ws);

  std::string error_message;

#pragma omp parallel for schedule(dynamic)
  for (auto& task : tasks) {
    try {
      Workspace local_ws{WorkspaceInitialization::Empty};
      task.share_workspace(local_ws, ws);
      task.execute(local_ws);
    } catch (std::exception& e) {
#pragma omp critical
      if (error_message.empty()) {
        error_message = std::format(R"(Failed to execute method "{}":

{})",
                                    task.get_name(),
                                    e.what());
      }
    }
  }

  if (not error_message.empty()) throw std::runtime_error(error_message);
} catch (std::exception& e) {
  throw std::runtime_error(std::format(R"(Cannot perform parallel execution of:

{}

{})",
                                       name,
                                       e.what()));
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

  for (auto& item : share) { out += std::format("{}Shares the global *{}*\n", prep, item); }

  for (auto& item : copy) { out += std::format("{}Copies the global *{}*\n", prep, item); }

  for (auto& method : methods) { out += std::format("{}{}\n", prep, method.sphinx_list_item()); }

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

  ARTS_USER_ERROR(std::format(R"(Agenda "{}" does not have a default input "{}")", this->name, name));
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
