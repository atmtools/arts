#pragma once

#include <array.h>

#include <ostream>
#include <vector>

class Method;
class Workspace;
struct Wsv;

class Agenda {
  std::string name;
  std::vector<Method> methods;
  std::vector<std::string> share{};
  std::vector<std::string> copy{};
  bool checked{false};

public:
  Agenda() = default;
  Agenda(std::string name);

  void add(const Method& method);

  //! Must be called before named agendas, will deal with input and output variables for copy_workspace
  void finalize();

  //! Copies the required workspace variables from the agenda
  void copy_workspace(Workspace& out, const Workspace& in, bool share_only=false) const;

  //! Copies the required workspace variables from the agenda
  [[nodiscard]] Workspace copy_workspace(const Workspace& in) const;

  //! Executes the agenda without checks on the current workspace
  void execute(Workspace& ws) const;

  [[nodiscard]] bool is_checked() const { return checked; }

  [[nodiscard]] const std::string& get_name() const {return name;}

  void set_name(const std::string& v) {name = v;}

  [[nodiscard]] bool has_method(const std::string& method) const;

  friend std::ostream& operator<<(std::ostream& os, const Agenda& a);
};

using ArrayOfAgenda = Array<Agenda>;
