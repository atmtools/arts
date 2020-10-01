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

#include "workspace_memory_handler.h"
#include "wsv_aux.h"

namespace global_data {
WorkspaceMemoryHandler workspace_memory_handler{};
}
using global_data::workspace_memory_handler;

Array<WsvRecord> Workspace::wsv_data;

map<String, Index> Workspace::WsvMap;

Workspace::Workspace()
    : ws(0)
#ifndef NDEBUG
      ,
      context("")
#endif
{
}

void Workspace::define_wsv_map() {
  for (Index i = 0; i < wsv_data.nelem(); ++i) {
    WsvMap[wsv_data[i].Name()] = i;
  }
}

Index Workspace::add_wsv(const WsvRecord &wsv) {
  Workspace::wsv_data.push_back(wsv);
  Workspace::define_wsv_map();
  return wsv_data.nelem() - 1;
}

void Workspace::del(Index i) {
  WsvStruct *wsvs = ws[i].top();

  if (wsvs && wsvs->wsv) {
    workspace_memory_handler.deallocate(wsv_data[i].Group(), wsvs->wsv);
    wsvs->wsv = NULL;
    wsvs->auto_allocated = false;
    wsvs->initialized = false;
  }
}

void Workspace::duplicate(Index i) {
  WsvStruct *wsvs = new WsvStruct;

  wsvs->auto_allocated = true;
  if (ws[i].size() && ws[i].top()->wsv) {
    wsvs->wsv = workspace_memory_handler.duplicate(wsv_data[i].Group(),
                                                   ws[i].top()->wsv);
    wsvs->initialized = true;
  } else {
    wsvs->wsv = NULL;
    wsvs->initialized = false;
  }
  ws[i].push(wsvs);
}

Workspace::Workspace(const Workspace &workspace) : ws(workspace.ws.nelem()) {
#ifndef NDEBUG
  context = workspace.context;
#endif
  for (Index i = 0; i < workspace.ws.nelem(); i++) {
    WsvStruct *wsvs = new WsvStruct;
    wsvs->auto_allocated = false;
    if (workspace.ws[i].size() && workspace.ws[i].top()->wsv) {
      wsvs->wsv = workspace.ws[i].top()->wsv;
      wsvs->initialized = workspace.ws[i].top()->initialized;
    } else {
      wsvs->wsv = NULL;
      wsvs->initialized = false;
    }
    ws[i].push(wsvs);
  }
}

Workspace::~Workspace() {
#ifndef NDEBUG
#pragma omp critical(ws_destruct)
  if (context != "") cout << "WS destruct: " << context << endl;
#endif
  for (int i = 0; i < ws.nelem(); i++) {
    WsvStruct *wsvs;

    while (ws[i].size()) {
      wsvs = ws[i].top();
      if (wsvs->auto_allocated && wsvs->wsv) {
        workspace_memory_handler.deallocate(wsv_data[i].Group(), wsvs->wsv);
      }
      delete (wsvs);
      ws[i].pop();
    }
  }
  ws.empty();
}

void *Workspace::pop(Index i) {
  WsvStruct *wsvs = ws[i].top();
  void *vp = NULL;
  if (wsvs) {
    vp = wsvs->wsv;
    delete wsvs;
    ws[i].pop();
  }
  return vp;
}

void Workspace::pop_free(Index i) {
  WsvStruct *wsvs = ws[i].top();

  if (wsvs) {
    if (wsvs->wsv)
      workspace_memory_handler.deallocate(wsv_data[i].Group(), wsvs->wsv);

    delete wsvs;
    ws[i].pop();
  }
}

void Workspace::push(Index i, void *wsv) {
  WsvStruct *wsvs = new WsvStruct;
  wsvs->auto_allocated = false;
  wsvs->initialized = true;
  wsvs->wsv = wsv;
  ws[i].push(wsvs);
}

void Workspace::push_uninitialized(Index i, void *wsv) {
  WsvStruct *wsvs = new WsvStruct;
  wsvs->auto_allocated = false;
  wsvs->initialized = false;
  wsvs->wsv = wsv;
  ws[i].push(wsvs);
}

void *Workspace::operator[](Index i) {
  if (!ws[i].size()) push(i, NULL);

  if (!ws[i].top()->wsv) {
    ws[i].top()->auto_allocated = true;
    ws[i].top()->wsv = workspace_memory_handler.allocate(wsv_data[i].Group());
  }

  ws[i].top()->initialized = true;

  return (ws[i].top()->wsv);
}
