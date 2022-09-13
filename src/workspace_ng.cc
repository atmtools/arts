/* Copyright (C) 2004-2012 Oliver Lemke <olemke@core-dump.info>

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

/** This file implements Workspace related functionality.
 *
 * @file   workspace_ng.cc
 * @author Oliver Lemke <olemke@core-dump.info>
 * @date   2004-11-05
 */

#include "workspace_ng.h"
#include "debug.h"
#include "workspace_global_data.h"
#include <memory>

#include "agenda_class.h"
#include "global_data.h"
#include "workspace_memory_handler.h"
#include "wsv_aux.h"

namespace global_data {
WorkspaceMemoryHandler workspace_memory_handler{};
}  // namespace global_data

using global_data::workspace_memory_handler;

Index Workspace::add_wsv(const WsvRecord &wsv) {
  wsv_data_ptr->push_back(wsv);
  WsvMap_ptr->operator[](wsv.Name()) = wsv_data_ptr->nelem() - 1;
  ws.emplace_back();
  return wsv_data_ptr->nelem() - 1;
}

void Workspace::set_empty(Index i) {
  if (ws[i].size()) {
    ws[i].pop();
    emplace(i);
  }
}

void Workspace::duplicate(Index i) {
  static const auto agenda_index = global_data::WsvGroupMap.at("Agenda");

  WorkspaceVariableStruct wsvs;

  if (ws[i].size()) {
    wsvs.wsv = workspace_memory_handler.duplicate((*wsv_data_ptr)[i].Group(),
                                                  ws[i].top().wsv);
    wsvs.initialized = true;
  } else {
    if ((*wsv_data_ptr)[i].Group() == agenda_index) {
      wsvs.wsv = std::make_shared<Agenda>(*this);
    } else {
      wsvs.wsv = workspace_memory_handler.allocate((*wsv_data_ptr)[i].Group());
    }
    wsvs.initialized = false;
  }
  ws[i].push(std::move(wsvs));
}

Workspace::Workspace(const Workspace &workspace)
    : ws(workspace.ws.nelem()),
      wsv_data_ptr(workspace.wsv_data_ptr),
      WsvMap_ptr(workspace.WsvMap_ptr),
      original_workspace(workspace.original_workspace) {
  for (Index i = 0; i < workspace.ws.nelem(); i++) {
    if (workspace.ws[i].size() && workspace.ws[i].top().wsv) {
      WorkspaceVariableStruct wsvs;
      wsvs.wsv = workspace.ws[i].top().wsv;
      wsvs.initialized = workspace.ws[i].top().initialized;
      ws[i].push(std::move(wsvs));
    }
  }
}

void Workspace::pop(Index i) {
  ws[i].pop(); }

void Workspace::swap(Workspace &other) noexcept {
  ws.swap(other.ws);
  wsv_data_ptr.swap(other.wsv_data_ptr);
  WsvMap_ptr.swap(other.WsvMap_ptr);
  std::swap(original_workspace, other.original_workspace);
}

bool Workspace::is_initialized(Index i) const {
  return ws[i].size() and ws[i].top().initialized;
}

Index Workspace::depth(Index i) const { return static_cast<Index>(ws[i].size()); }

void Workspace::emplace(Index i) {
  static const auto agenda_index = global_data::WsvGroupMap.at("Agenda");

  if ((*wsv_data_ptr)[i].Group() == agenda_index) {
    ws[i].emplace(
        WorkspaceVariableStruct{std::make_shared<Agenda>(*this), false});
  } else {
    ws[i].emplace(WorkspaceVariableStruct{
        workspace_memory_handler.allocate((*wsv_data_ptr)[i].Group()), false});
  }
}

std::shared_ptr<void> Workspace::operator[](Index i) {
  if (ws[i].size() == 0) emplace(i);
  ws[i].top().initialized = true;
  return ws[i].top().wsv;
}

Workspace::Workspace()
    : ws(global_data::wsv_data.nelem()),
      wsv_data_ptr(std::make_shared<Array<WsvRecord>>(global_data::wsv_data)),
      WsvMap_ptr(std::make_shared<map<String, Index>>(global_data::WsvMap)),
      original_workspace(this) {
  ARTS_ASSERT(wsv_data_ptr -> size() == WsvMap_ptr->size())

  for (Index i = 0; i < (*wsv_data_ptr).nelem(); i++) {
    if ((*wsv_data_ptr)[i].has_defaults()) {
      push_move(i, (*wsv_data_ptr)[i].get_copy());
    }
  }
}

std::shared_ptr<Workspace> Workspace::deepcopy() {
  std::shared_ptr<Workspace> out{new Workspace{}};
  out->wsv_data_ptr = std::shared_ptr<Workspace::wsv_data_type>(
      new Workspace::wsv_data_type{*wsv_data_ptr});
  out->WsvMap_ptr = std::shared_ptr<Workspace::WsvMap_type>(
      new Workspace::WsvMap_type{*WsvMap_ptr});
  out->ws.resize(nelem());

  for (Index i = 0; i < out->nelem(); i++) {
    auto &wsv_data = out->wsv_data_ptr->operator[](i);

    if (depth(i) > 0) {
      // Set the WSV by copying the top value
      out->ws[i].emplace(WorkspaceVariableStruct{
          workspace_memory_handler.duplicate(
              wsv_data_ptr->operator[](i).Group(), ws[i].top().wsv),
          is_initialized(i)});

      // Copy the agenda to the new workspace
      if (wsv_data.Group() == WorkspaceGroupIndexValue<Agenda>) {
        Agenda *ag = static_cast<Agenda *>(out->operator[](i).get());
        *ag = ag->deepcopy_if(*out);
      } else if (wsv_data.Group() == WorkspaceGroupIndexValue<ArrayOfAgenda>) {
        for (auto &a :
             *static_cast<ArrayOfAgenda *>(out->operator[](i).get())) {
          a = a.deepcopy_if(*out);
        }
      }
    }

    // If we have any default agenda types, we must copy them to the new workspace as well
    if (wsv_data.has_defaults()) {
      if (wsv_data.Group() == WorkspaceGroupIndexValue<Agenda>) {
        wsv_data.update_default_value(
            Agenda(wsv_data.default_value()).deepcopy_if(*out));
      }
      if (wsv_data.Group() == WorkspaceGroupIndexValue<ArrayOfAgenda>) {
        wsv_data.update_default_value(
            deepcopy_if(*out, wsv_data.default_value()));
      }
    }
  }

  return out;
}
