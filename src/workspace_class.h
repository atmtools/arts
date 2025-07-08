#pragma once

#include <wsv_value_wrapper.h>
#include <format_tags.h>

#include <memory>
#include <unordered_map>

enum class WorkspaceInitialization : bool { FromGlobalDefaults, Empty };

struct Workspace {
  std::unordered_map<std::string, Wsv> wsv;

  Workspace(WorkspaceInitialization how_to_initialize =
                WorkspaceInitialization::FromGlobalDefaults);

  Workspace(std::unordered_map<std::string, Wsv>);

  //! Returns a shared pointer to the workspace variable with the given name.
  [[nodiscard]] const Wsv& share(const std::string& name) const;

  //! Returns a copy of the workspace variable with the given name.
  [[nodiscard]] Wsv copy(const std::string& name) const;

  //! Sets the workspace variable with the given name to the given value.
  void set(const std::string& name, const Wsv& data);

  //! Copy the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, const T& data) {
    set(name, Wsv{data});
  }

  //! Move the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, T&& data) {
    set(name, Wsv{std::move(data)});
  }

  //! Borrows the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void set(const std::string& name, T* data) {
    set(name, Wsv{data});
  }

  //! Overwrites a variable with another variable of the same type
  void overwrite(const std::string& name, const Wsv& data);

  //! Copy the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, const T& data) {
    overwrite(name, Wsv{data});
  }

  //! Move the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, T&& data) {
    overwrite(name, Wsv{std::move(data)});
  }

  //! Borrows the workspace variable with the given name to the given value.
  template <WorkspaceGroup T>
  void overwrite(const std::string& name, T* data) {
    overwrite(name, Wsv{data});
  }

  //! Returns a type directly based on the name of the workspace variable.
  template <WorkspaceGroup T>
  [[nodiscard]] T& get(const std::string& name) const;

  //! Returns a type directly based on the name of the workspace variable, creating it in-place if it is not there
  template <WorkspaceGroup T>
  [[nodiscard]] std::shared_ptr<T> share_or(const std::string& name);

  //! Returns a type directly based on the name of the workspace variable, creating it in-place if it is not there
  template <WorkspaceGroup T>
  [[nodiscard]] T& get_or(const std::string& name);

  //! Checks if the workspace variable with the given name exists.
  [[nodiscard]] bool contains(const std::string& name) const;

  //! As contains, but also checks if the variable is a workspace variable
  [[nodiscard]] bool wsv_and_contains(const std::string& name) const;

  friend std::ostream& operator<<(std::ostream& os, const Workspace& ws);

  [[nodiscard]] auto begin() { return wsv.begin(); }

  [[nodiscard]] auto end() { return wsv.end(); }

  [[nodiscard]] auto begin() const { return wsv.begin(); }

  [[nodiscard]] auto end() const { return wsv.end(); }

  [[nodiscard]] auto size() const { return wsv.size(); }

  void init(const std::string& name);

  [[nodiscard]] Workspace deepcopy() const;
};

template <>
struct std::formatter<Workspace> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Workspace& v, FmtContext& ctx) const {
    const std::string_view sep   = tags.sep(true);
    const std::string_view quote = tags.quote();

    tags.add_if_bracket(ctx, '{');
    for (auto& [name, wsv] : v) {
      tags.format(ctx, quote, name, quote, ": "sv, wsv, sep);
    }
    tags.add_if_bracket(ctx, '}');

    return ctx.out();
  }
};

template <>
struct xml_io_stream_name<Workspace> {
  static constexpr std::string_view name = "Workspace";
};

template <>
struct xml_io_stream_aggregate<Workspace> {
  static constexpr bool value = true;
};
