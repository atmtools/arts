#pragma once

#include <bifstream.h>
#include <bofstream.h>

#include <istream>
#include <ostream>
#include <span>
#include <string_view>

using namespace std::literals;

//! Name overload for all the automatic methods
template <typename T>
struct xml_io_stream_name {
  static constexpr std::string_view name = "<unknown>"sv;
};

template <typename T>
constexpr std::string_view xml_io_stream_name_v = xml_io_stream_name<T>::name;

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
  constexpr static std::string_view type_name = xml_io_stream_name_v<T>;

  // Core IO
  static void write(std::ostream&,
                    const T&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv)                  = delete;
  static void read(std::istream&, T&, bifstream* = nullptr)   = delete;
  static void extend(std::istream&, T&, bifstream* = nullptr) = delete;

  // Binary streaming IO (optional)
  static void get(std::span<const T>, bofstream*) = delete;
  static void put(std::span<T>, bifstream*)       = delete;

  // Text streaming (optional)
  static void parse(std::span<T>, std::istream&) = delete;
};

//! Test that the type can be written via XML-IO
template <typename T>
concept xml_io_writable = requires(T a) {
  xml_io_stream<T>::write(std::declval<std::ostream&>(), a, nullptr, "");
  xml_io_stream<T>::write(std::declval<std::ostream&>(), a);
};

//! Test that the type can be read via XML-IO
template <typename T>
concept xml_io_readable = requires(T a) {
  xml_io_stream<T>::read(std::declval<std::istream&>(), a, nullptr);
  xml_io_stream<T>::read(std::declval<std::istream&>(), a);
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

//! Test that the type can be extended via XML-IO
template <typename T>
concept arts_xml_extendable = arts_xml_ioable<T> and requires(T a) {
  xml_io_stream<T>::extend(std::declval<std::istream&>(), a, nullptr);
  xml_io_stream<T>::extend(std::declval<std::istream&>(), a);
};

//! Test that the type is ready for binary data also
template <typename T>
concept xml_io_binary =
    arts_xml_ioable<T> and requires(std::span<const T, 1> a) {
      xml_io_stream<T>::put(a, nullptr);
    } and requires(std::span<T> a) { xml_io_stream<T>::get(a, nullptr); };

//! Test if parsing contiguous is possible
// For most types this means doing the reverse of std::println(os, "{:IO}", x).
// This is not an explicit demand since, e.g., Numeric cannot, being builtin
template <typename T>
concept xml_io_parseable = requires(std::span<T> b, std::istream& is) {
  xml_io_stream<T>::parse(b, is);
};
