#pragma once

#include <auto_wsg.h>

struct WsvValueWrapper;

class Wsv {
  WsvValueWrapper* data{nullptr};

 public:
  Wsv();
  Wsv(const Wsv&);
  Wsv(Wsv&&) noexcept;
  Wsv& operator=(const Wsv&);
  Wsv& operator=(Wsv&&) noexcept;
  ~Wsv();

  template <WorkspaceGroup T>
  Wsv(T&&) noexcept;
  template <WorkspaceGroup T>
  Wsv(const T&);
  template <WorkspaceGroup T>
  Wsv(T*) noexcept;
  template <WorkspaceGroup T>
  Wsv(std::shared_ptr<T>&&) noexcept;

  [[nodiscard]] Wsv copy() const;

  [[nodiscard]] std::string_view type_name() const;

  [[nodiscard]] std::size_t index() const;

  [[nodiscard]] bool holds_same(const Wsv& other) const;

  template <WorkspaceGroup T>
  [[nodiscard]] bool holds() const;

  template <WorkspaceGroup T>
  [[nodiscard]] const std::shared_ptr<T>& share_unsafe() const;

  template <WorkspaceGroup T>
  [[nodiscard]] const std::shared_ptr<T>& share() const;

  template <WorkspaceGroup T>
  [[nodiscard]] T& get_unsafe() const;

  template <WorkspaceGroup T>
  [[nodiscard]] T& get() const;

  [[nodiscard]] static Wsv from_named_type(const std::string& type);

  [[nodiscard]] const WsvValue& value() const;

  [[nodiscard]] WsvValue& value();
};
