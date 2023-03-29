#include <iostream>
#include <map>

#include "global_data.h"
#include "methods.h"
#include "workspace.h"

struct Types {
  Index type;
  String val;
};

bool is_default(const String& x) { return x == NODEF; }

int main() {
  define_wsv_groups();
  define_wsv_data();
  define_wsv_map();
  define_md_data_raw();

  std::ofstream("default_gins.h") << R"--(#pragma once

#include "workspace.h"

Index create_workspace_gin_default_internal(Workspace& ws, const std::string_view method, const std::string_view gin);
)--";

  std::map<String, Types> has;

  for (auto& method : global_data::md_data_raw) {
    for (std::size_t i = 0; i < method.GIn().size(); i++) {
      if (auto& def = method.GInDefault()[i]; not is_default(def)) {
        auto group_ind = method.GInType()[i];
        auto& group = global_data::wsv_groups[group_ind].name;
        auto& vals =
            has[var_string("::", method.Name(), "::", method.GIn()[i])];
        vals.type = group_ind;
        auto& val = vals.val;

        if (group == "String") {
          val = std::string("\"") + def + std::string("\"");
        } else if (group == "Numeric") {
          if ("NaN" == def or "nan" == def) {
            val = "std::numeric_limits<Numeric>::quiet_NaN()";
          } else if ("Inf" == def or "inf" == def) {
            val = "std::numeric_limits<Numeric>::infinity()";
          } else if ("-Inf" == def or "-inf" == def) {
            val = "-std::numeric_limits<Numeric>::infinity()";
          } else {
            val = def;
          }
        } else if (def == "[]") {
          val = var_string(group, "{}");
        } else if (group == "ArrayOfIndex") {
          if (def[0] == '[')
            val = "std::initializer_list<Index>" + def;
          else
            val = def;
        } else if (group == "Vector") {
          if (def[0] == '[')
            val = "std::initializer_list<Numeric>" + def;
          else
            val = def;
        } else {
          val = def;
        }

        if (val == "") val = var_string(group, "{}");

        for (auto& x : val) {
          if (x == '[')
            x = '{';
          else if (x == ']')
            x = '}';
        }
      }
    }
  }

  std::ofstream os("default_gins.cc");

  os << R"--(#include "auto_md.h"

template <typename T, Index group>
Index get_and_set_wsv_gin_pos(Workspace& ws, Index pos, T&& data) {
  static std::size_t anon=0;

  if (pos < 0) {
    pos = ws.add_wsv(WsvRecord(var_string("::defgin", anon++).c_str(), "do not modify", group));
  }

  *static_cast<T *>(ws[pos].get()) = std::forward<T>(data);
  
  return pos;
}

Index create_workspace_gin_default_internal(Workspace& ws, const std::string_view method, const std::string_view gin) {
  const String key{var_string("::", method, "::", gin)};
  const static std::map<String, Index> gins {
)--";

  {
    Index counter = 0;
    for (auto& [key, items] : has) {
      os << "    {\"" << key << "\", " << counter++ << "},\n";
    }
  }

  os << R"--(  };

  auto ptr = ws.WsvMap_ptr -> find(key);
  Index pos = ptr == ws.WsvMap_ptr -> end() ? -1 : ptr -> second;

  switch (gins.at(key)) {
)--";
  {
    Index counter = 0;
    for (auto& [key, items] : has) {
      os << "    case " << counter++ << ": return get_and_set_wsv_gin_pos<"
         << global_data::wsv_groups[items.type] << ", " << items.type
         << ">(ws, pos, " << items.val << ");\n";
    }
  }
  os << R"--(  }

  return pos;
}
)--";
}