#pragma once

#include <array.h>

#include <ostream>
#include <vector>

class Method;
class Workspace;
class Wsv;

class Agenda {
  std::string name;
  std::vector<Method> methods;
  std::vector<std::string> output{};
  std::vector<std::string> pure_input{};
  bool checked{false};

public:
Agenda() = default;
Agenda(std::string name);

void add(const Method& method);

void finalize();
void copy_workspace(Workspace& out, const Workspace& in) const;
[[nodiscard]] Workspace copy_workspace(const Workspace& in) const;
void execute(Workspace& ws) const;
[[nodiscard]] bool is_checked() const { return checked; }
[[nodiscard]] const std::string& get_name() const {return name;}
[[nodiscard]] bool has_method(const std::string& method) const;

friend std::ostream& operator<<(std::ostream& os, const Agenda& a);
};

using ArrayOfAgenda = Array<Agenda>;
