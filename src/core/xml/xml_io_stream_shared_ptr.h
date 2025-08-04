#pragma once

#include <memory>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename T>
struct xml_io_stream_name<std::shared_ptr<T>> {
  static constexpr std::string_view name = "Shared"sv;
};

template <typename T>
  requires(arts_xml_ioable<std::remove_const_t<T>> and
           std::is_default_constructible_v<T>)
struct xml_io_stream<std::shared_ptr<T>> {
  using mutT = std::remove_const_t<T>;

  static constexpr std::string_view type_name =
      xml_io_stream_name_v<std::shared_ptr<T>>;

  static void write(std::ostream& os,
                    const std::shared_ptr<T>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    XMLTag tag{type_name,
               "name",
               name,
               "type",
               xml_io_stream<T>::type_name,
               "notnull",
               Index{bool{x}}};
    tag.write_to_stream(os);

    if (x) xml_io_stream<mutT>::write(os, *x, pbofs, name);

    tag.write_to_end_stream(os);
  }

  static void read(std::istream& is,
                   std::shared_ptr<T>& x,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);
    tag.check_attribute("type", xml_io_stream<mutT>::type_name);

    Index notnull;
    tag.get_attribute_value("notnull", notnull);

    if (notnull) {
      mutT v{};
      xml_io_stream<mutT>::read(is, v, pbifs);
      x = std::make_shared<T>(std::move(v));
    } else {
      x = nullptr;
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::format("Error reading {}<{}>:\n{}",
                                         type_name,
                                         xml_io_stream<T>::type_name,
                                         e.what()));
  }
};
