#include "workspace_groups.h"

#include <iostream>
#include <stdexcept>

WorkspaceGroupRecord::WorkspaceGroupRecord(std::string&& f, std::string&& d) try
    : file(std::move(f)), desc(std::move(d)) {
  if (file.empty()) throw std::runtime_error("File name must not be empty");
  if (desc.empty()) throw std::runtime_error("Description must not be empty");
  if (desc.back() != '\n') desc.push_back('\n');
} catch (std::exception& e) {
  std::cerr << R"--(
      ERROR

        CANNOT CREATE WORKSPACE GROUP.  THIS IS A BUG.

      ERROR

Inputs:
  file:)--" << f
            << "\ndesc: " << d << "\nExtra:\n"
            << e.what() << '\n';
  std::exit(EXIT_FAILURE);
}

std::unordered_map<std::string, WorkspaceGroupRecord>
internal_workspace_groups() {
  std::unordered_map<std::string, WorkspaceGroupRecord> wsv_groups;

  wsv_groups["Any"] = WorkspaceGroupRecord {
    "supergeneric.h",
    "Test any"
  };

  wsv_groups["Index"] = WorkspaceGroupRecord {
    "matpack.h",
    "Test index"
  };

   wsv_groups["Numeric"] = WorkspaceGroupRecord {
    "matpack.h",
    "Test numeric"
  };

   wsv_groups["Vector"] = WorkspaceGroupRecord {
    "matpack.h",
    "Test vector"
  };

  return wsv_groups;
}
