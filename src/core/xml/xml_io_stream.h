#pragma once

#include <bifstream.h>
#include <bofstream.h>
#include <configtypes.h>

#include <istream>
#include <ostream>
#include <string_view>

using namespace std::literals;

/** Basic XML I/O interface.
 * 
 * @tparam T The type to store
 * 
 * The default methods are all deleted.  You need to specialize this for your type.
 * 
 * The first two are required.  The second two allows some optimization for binary-layout compatible types
 */
template <typename T>
struct xml_io_stream {
  constexpr static std::string_view type_name = "<unknown>"sv;

  // Core IO
  static void write(std::ostream&,
                    const T&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv)                = delete;
  static void read(std::istream&, T&, bifstream* = nullptr) = delete;

  // Binary/streaming IO (optional)
  static void get(const T* const, bofstream*, Size) = delete;
  static void put(T*, bifstream*, Size)             = delete;
};

//! Test that the type can be written via XML-IO
template <typename T>
concept xml_io_writable = requires(T a) {
  xml_io_stream<T>::write(std::declval<std::ostream&>(), a, nullptr, "");
};

//! Test that the type can be read via XML-IO
template <typename T>
concept xml_io_readable = requires(T a) {
  xml_io_stream<T>::read(std::declval<std::istream&>(), a, nullptr);
};

//! Test that the type can be read via XML-IO
template <typename T>
concept xml_io_nameable = requires(xml_io_stream<T> a) {
  { a.type_name } -> std::convertible_to<std::string_view>;
};

//! Test that the type is ready
template <typename T>
concept arts_xml_ioable =
    xml_io_writable<T> and xml_io_readable<T> and xml_io_nameable<T>;

//! Test that the type is ready for binary data also
template <typename T>
concept xml_io_binary = arts_xml_ioable<T> and requires(const T* const a) {
  xml_io_stream<T>::put(a, nullptr, 1);
} and requires(T* a) { xml_io_stream<T>::get(a, nullptr, 1); };

template <typename T>
concept bitshift_readable = requires(T a, std::istream& is, std::ostream& os) {
  is >> a;
  std::println(os, "{:IO}", a);
};
