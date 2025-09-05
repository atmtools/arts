#pragma once

#include <array.h>
#include <xml.h>

#include <vector>

class Method;
struct Workspace;
class Wsv;

class Agenda {
  std::string name{};
  std::vector<Method> methods;
  std::vector<std::string> share{};
  std::vector<std::string> copy{};
  bool checked{false};

 public:
  Agenda(std::string name = "not-a-name");

  Agenda(std::string name,
         const std::vector<Method>& methods,
         const std::vector<std::string>& share,
         const std::vector<std::string>& copy,
         bool checked);

  void add(const Method& method);

  //! Must be called before named agendas, will deal with input and output variables for copy_workspace
  void finalize(bool fix = false);

  //! Copies the required workspace variables from the agenda
  void copy_workspace(Workspace& out,
                      const Workspace& in,
                      bool share_only = false) const;

  //! Copies the required workspace variables from the agenda
  [[nodiscard]] Workspace copy_workspace(const Workspace& in) const;

  //! Executes the agenda without checks on the current workspace
  void execute(Workspace& ws) const;

  [[nodiscard]] bool is_checked() const { return checked; }

  [[nodiscard]] const std::string& get_name() const { return name; }

  void set_name(const std::string& v) { name = v; }

  [[nodiscard]] bool has_method(const std::string& method) const;

  [[nodiscard]] const std::vector<Method>& get_methods() const {
    return methods;
  }
  [[nodiscard]] const std::vector<std::string>& get_share() const {
    return share;
  }
  [[nodiscard]] const std::vector<std::string>& get_copy() const {
    return copy;
  }

  [[nodiscard]] std::string sphinx_list(
      const std::string_view prep = "- ") const;
};

template <>
struct xml_io_stream<Agenda> {
  static constexpr std::string_view type_name = "Agenda"sv;

  static void write(std::ostream& os,
                    const Agenda& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, Agenda& x, bifstream* pbifs = nullptr);
};
