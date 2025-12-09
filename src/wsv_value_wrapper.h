#pragma once

//! auto-generated file, edits do not persist

#include <auto_wsg.h>

class Wsv {
  std::shared_ptr<void> data;  // A shared_ptr custom variant
  std::size_t index;           // The index of the variant

 public:
  std::size_t value_index() const noexcept { return index; }

  //! Default constructor, creates an Any
  Wsv() : Wsv{std::shared_ptr<Any>(new Any{})} {}

  Wsv(Wsv&&) noexcept            = default;
  Wsv& operator=(Wsv&&) noexcept = default;

  //! We must be able to make a Wsv from any of the allowed workspace groups
  template <WorkspaceGroup T>
  Wsv(std::shared_ptr<T>&& x) noexcept
      : data{std::move(x)}, index{WorkspaceGroupInfo<T>::index} {}

  //! We must be able to take over the data of a workspace group
  template <WorkspaceGroup T>
  Wsv(T&& x) noexcept : Wsv(std::shared_ptr<T>(new T{std::move(x)})) {}

  //! We must have a sneaky share constructor (note: this is unsafe)
  template <WorkspaceGroup T>
  Wsv(T* x) noexcept : Wsv(std::shared_ptr<T>(x, [](void*) {})) {}

  //! We must be able to copy workspace groups into the variant
  template <WorkspaceGroup T>
  Wsv(const T& x) : Wsv(std::shared_ptr<T>(new T{x})) {}

  Wsv& operator=(const Wsv& x) { return *this = Wsv(x); }

  [[nodiscard]] bool holds_same(const Wsv& x) const noexcept {
    return index == x.index;
  }

  template <WorkspaceGroup T>
  [[nodiscard]] bool holds() const noexcept {
    return index == WorkspaceGroupInfo<T>::index;
  }

  //! Unsafely share into the variant (this may segfault)
  template <WorkspaceGroup T>
  [[nodiscard]] std::shared_ptr<T> share_unsafe() const {
    assert(holds<T>());
    return std::static_pointer_cast<T>(data);
  }

  //! Safely share into the variant (this may throw)
  template <WorkspaceGroup T>
  [[nodiscard]] std::shared_ptr<T> share() const {
    if (not holds<T>())
      throw std::runtime_error(std::format(
          R"(Casting "{}" to "{}")", WorkspaceGroupInfo<T>::name, type_name()));
    return share_unsafe<T>();
  }

  //! Unsafely reference into the variant (this may segfault)
  template <WorkspaceGroup T>
  [[nodiscard]] T& get_unsafe() const {
    assert(holds<T>());
    return *share_unsafe<T>();
  }

  //! Safely reference into the variant (this may throw)
  template <WorkspaceGroup T>
  [[nodiscard]] T& get() const {
    if (not holds<T>())
      throw std::runtime_error(std::format(
          R"(Casting "{}" to "{}")", WorkspaceGroupInfo<T>::name, type_name()));
    return get_unsafe<T>();
  }

  Wsv(const Wsv&) = default;

  //! Get a Wsv that holds the named type
  static Wsv from_named_type(const std::string_view);

  //! Get a Wsv that holds the named type
  static Wsv from_index(Index);

  //! Get the type name of the variant
  [[nodiscard]] std::string_view type_name() const;

  //! Get the type name of the variant
  [[nodiscard]] std::string vformat(std::string_view) const;

  [[nodiscard]] Wsv shared() const { return *this; }
  [[nodiscard]] Wsv copied() const;

  //! Write to an existing stream (only the type is written, nothing about Wsv)
  std::ostream& write_to_stream(std::ostream&,
                                bofstream*,
                                const std::string&) const;

  //! Write to a file (only the type is written, nothing about Wsv)
  std::string write_to_file(const std::string& fn,
                            FileType ftype,
                            bool clobber) const;

  //! Read from an existing stream (only the type is read, nothing about Wsv)
  std::istream& read_from_stream(std::istream&, bifstream*);

  //! Read from a file (only the type is read, nothing about Wsv)
  std::string read_from_file(const std::string& fn);
};
