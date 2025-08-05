#include "callback.h"

#include <workspace.h>

#include "workspace_class.h"

void CallbackOperator::operator()(Workspace& ws_in) const try {
  ARTS_USER_ERROR_IF(
      not callback.f, "No callback function set for operator:\n{}", *this);

  Workspace ws(WorkspaceInitialization::Empty);

  for (auto& n : inputs) ws.set(n, ws_in.share(n));
  for (auto& n : outputs) {
    if (ws_in.contains(n)) {
      ws.set(n, ws_in.share(n));
    } else {
      ws.init(n);
      ws_in.set(n, ws.share(n));
    }
  }
  callback(ws);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error in callback operator:\n{}\n{}", *this, e.what()));
}

void xml_io_stream<CallbackOperator>::write(std::ostream& os,
                                            const CallbackOperator& x,
                                            bofstream* pbofs,
                                            std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.callback, pbofs);
  xml_write_to_stream(os, x.inputs, pbofs);
  xml_write_to_stream(os, x.outputs, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<CallbackOperator>::read(std::istream& is,
                                           CallbackOperator& x,
                                           bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.callback, pbifs);
  xml_read_from_stream(is, x.inputs, pbifs);
  xml_read_from_stream(is, x.outputs, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
