#pragma once

#include <auto_wsg.h>

class Wsv {
  template <WorkspaceGroup T>
  static constexpr std::size_t ind = WorkspaceGroupInfo<T>::index;

  template <WorkspaceGroup T>
  static constexpr std::string_view name = WorkspaceGroupInfo<T>::name;

  std::shared_ptr<void> data{std::shared_ptr<Any>{new Any{}}};
  std::size_t index{ind<Any>};

 public:
  [[nodiscard]] std::size_t value_index() const noexcept { return index; }

  Wsv()                          = default;
  Wsv(const Wsv&)                = default;
  Wsv(Wsv&&) noexcept            = default;
  Wsv& operator=(const Wsv&)     = default;
  Wsv& operator=(Wsv&&) noexcept = default;

  //! Share
  template <WorkspaceGroup T>
  Wsv(std::shared_ptr<T>&& x) noexcept : data{std::move(x)}, index{ind<T>} {}

  //! Take over
  template <WorkspaceGroup T>
  Wsv(T&& x) noexcept : Wsv(std::shared_ptr<T>(new T{std::forward<T>(x)})) {}

  //! Share (for agenda sharing)
  template <WorkspaceGroup T>
  Wsv(T* x) noexcept : Wsv(std::shared_ptr<T>(x, [](void*) {})) {}

  //! Copy
  template <WorkspaceGroup T>
  Wsv(const T& x) : Wsv(std::shared_ptr<T>(new T{x})) {}

  [[nodiscard]] bool holds_same(const Wsv& x) const noexcept {
    return index == x.index;
  }

  template <WorkspaceGroup T>
  [[nodiscard]] bool holds() const noexcept {
    return index == ind<T>;
  }

  //! Unsafely share into the variant
  template <WorkspaceGroup T>
  [[nodiscard]] std::shared_ptr<T> share_unsafe() const {
    assert(holds<T>());
    return std::static_pointer_cast<T>(data);
  }

  //! Safely share into the variant
  template <WorkspaceGroup T>
  [[nodiscard]] std::shared_ptr<T> share() const {
    if (not holds<T>())
      throw std::runtime_error(
          std::format(R"(Casting "{}" to "{}")", name<T>, type_name()));
    return share_unsafe<T>();
  }

  //! Unsafely reference into the variant
  template <WorkspaceGroup T>
  [[nodiscard]] T& get_unsafe() const {
    assert(holds<T>());
    return *share_unsafe<T>();
  }

  //! Safely reference into the variant
  template <WorkspaceGroup T>
  [[nodiscard]] T& get() const {
    if (not holds<T>())
      throw std::runtime_error(
          std::format(R"(Casting "{}" to "{}")", name<T>, type_name()));
    return get_unsafe<T>();
  }

  //! Get a Wsv that holds the named type
  static Wsv from_named_type(const std::string_view);

  //! Get a Wsv that holds the indexed type
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
