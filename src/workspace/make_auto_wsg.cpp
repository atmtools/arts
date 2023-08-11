#include <algorithm>
#include <iostream>
#include <ranges>
#include <vector>

#include "workspace_groups.h"

const static auto data = internal_workspace_groups();

std::vector<std::pair<std::string, std::vector<std::string>>> files() {
  std::vector<std::pair<std::string, std::vector<std::string>>> files;

  for (const auto& [group_name, group] : data) {
    auto ptr = std::find_if(files.begin(),
                            files.end(),
                            [f = group.file](auto& p) { return f == p.first; });
    if (ptr == files.end()) {
      files.emplace_back(group.file, std::vector<std::string>{group_name});
    } else {
      ptr->second.push_back(group_name);
    }
  }

  std::sort(files.begin(), files.end(), [](auto& a, auto& b) {
    return a.first < b.first;
  });
  for (auto& [file, groups] : files) {
    std::sort(groups.begin(), groups.end());
  }

  return files;
}

std::vector<std::string> groups() {
  std::vector<std::string> groups;
  groups.reserve(data.size());

  for (const auto& [group_name, group] : data) {
    groups.push_back(group_name);
  }

  std::sort(groups.begin(), groups.end());

  return groups;
}

void header(std::ostream& os) {
  os << R"--(#pragma once

//! auto-generated by make_auto_wsg.cpp

#include <memory>
#include <variant>

)--";

  for (const auto& [file, groups] : files()) {
    os << "// ";
    for (const auto& group : groups) {
      os << group << ", ";
    }
    os << "\n#include <" << file << ">\n\n";
  }

  os << "template <typename T>\nconcept WorkspaceGroup =\n     std::is_same_v<T, "
     << groups()[0] << ">";
  for (auto& group : std::ranges::drop_view{groups(), 1}) {
    os << "\n  || std::is_same_v<T, " << group << ">";
  }
  os << ";\n\n";

  os << R"(template <typename T> struct WorkspaceGroupInfo {
  static constexpr std::string_view name = "<Unknown>";
  static constexpr std::string_view file = "<Unknown>";
  static constexpr std::string_view desc = "<Unknown>";
};

)";
  for (auto& group : groups()) {
    os << "template <> struct WorkspaceGroupInfo<" << group << "> {\n";
    os << "  static constexpr std::string_view name = \"" << group << "\";\n";
    os << "  static constexpr std::string_view file = \"" << data.at(group).file
       << "\";\n";
    os << "  static constexpr std::string_view desc = R\"--("
       << data.at(group).desc << ")--\";\n";
    os << "};\n\n";
  }

  os << "using WsvValue = std::variant<\n  std::shared_ptr<" << groups()[0]
     << ">";
  for (auto& group : std::ranges::drop_view{groups(), 1}) {
    os << ",\n  std::shared_ptr<" << group << ">";
  }
  os << ">;\n\n";

  os << "[[nodiscard]] std::string_view name_wsg(const WsvValue&);\n";
  os << "[[nodiscard]] bool valid_wsg(const std::string_view&);\n\n";

  os << R"--(
struct Wsv {
  WsvValue value;

  //! Move value into the workspace variable
  template <WorkspaceGroup T> Wsv(T&& x) : value(std::make_shared<T>(std::move(x))) {}

  //! Copy value into the workspace variable
  template <WorkspaceGroup T> Wsv(const T& x) : value(std::make_shared<T>(x)) {}

  //! Borrow value as workspace variable
  template <WorkspaceGroup T> Wsv(T* x) : value(std::shared_ptr<T>(x, [](void*){})) {}

  //! Share value as workspace variable
  template <WorkspaceGroup T> Wsv(std::shared_ptr<T>&& x) : value(std::move(x)) {}

  //! Must declare destructor to avoid incomplete type error
  Wsv() = default;
  Wsv(const Wsv&) = default;
  Wsv(Wsv&&) = default;
  Wsv& operator=(const Wsv&) = default;
  Wsv& operator=(Wsv&&) = default;

  [[nodiscard]] const std::string_view type_name() const;

  [[nodiscard]] constexpr std::size_t index() const {return value.index();}
  [[nodiscard]] constexpr bool holds_same(const Wsv& other) const {return index() == other.index();}

  template <WorkspaceGroup T>
  [[nodiscard]] constexpr bool holds() const {
    return std::holds_alternative<std::shared_ptr<T>>(value);
  }

  [[nodiscard]] Wsv copy() const;

  template <WorkspaceGroup T>
  [[nodiscard]] std::shared_ptr<T> share_unsafe() const {
    return *std::get_if<std::shared_ptr<T>>(&value);
  }

  template <WorkspaceGroup T>
  [[nodiscard]] std::shared_ptr<T> share() const {
    if (not holds<T>()) {
      throw std::runtime_error(var_string("Cannot use workspace variable of workspace group ",
                                          std::quoted(type_name()),
                                          " as ",
                                          std::quoted(WorkspaceGroupInfo<T>::name)));
    }

    return share_unsafe<T>();
  }

  template <WorkspaceGroup T>
  [[nodiscard]] T& get_unsafe() const {
    return *share_unsafe<T>();
  }

  template <WorkspaceGroup T>
  [[nodiscard]] T& get() const {
    return *share<T>();
  }

  static Wsv from_named_type(const std::string& type);
};
)--";
}

void implementation(std::ostream& os) {
  os << R"--(//! auto-generated by make_auto_wsg.cpp

#include <typeinfo>

#include "auto_wsg.h"

#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

)--";

  os << "std::string_view name_wsg(const WsvValue& x) {\n  const static std::unordered_map<std::size_t, std::string> val {\n";
  for (auto& group : groups()) {
    os << "    {Wsv(" << group << "{}).index(), \"" << group << "\"}, \n";
  }
  os << R"--(  };

  return val.at(x.index());
}

)--";

  os << "bool valid_wsg(const std::string_view& x) {\n  constexpr static std::array val {\n";
  for (auto& group : groups()) {
    os << "    \"" << group << "\", \n";
  }

  os << R"--(  };
  return std::binary_search(val.begin(), val.end(), x);
}

Wsv Wsv::copy() const {
  return std::visit([](const auto& val) -> Wsv {
    return Wsv{*val};
  }, value);
}

const std::string_view Wsv::type_name() const { return name_wsg(value); }

Wsv Wsv::from_named_type(const std::string& type) {
)--";

  for (auto& group : groups()) {
    os << "  if (type == \"" << group << "\") return " << group << "{};\n";
  }

os << R"--( throw std::runtime_error(var_string("Unknown workspace group ", std::quoted(type)));
}
)--";
}

int main() {
  header(std::cout);
  implementation(std::cerr);
}
