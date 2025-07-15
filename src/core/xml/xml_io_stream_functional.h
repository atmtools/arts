#pragma once

#include <concepts>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <variant>

#include "debug.h"
#include "xml_io_base.h"
#include "xml_io_stream.h"
#include "xml_io_stream_core.h"
#include "xml_io_stream_variant.h"

template <typename R, typename... Ts>
struct xml_io_stream_name<std::function<R(Ts...)>> {
  static constexpr std::string_view name = "Function"sv;
};

template <typename T>
struct xml_io_stream_functional {
  using func_t    = void (*)();
  using structs_t = void;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <typename R, typename... Ts>
struct xml_io_stream_functional<std::function<R(Ts...)>> {
  using func_t    = R (*)(Ts...);
  using structs_t = void;  // Overload to std::variant to use
  static constexpr std::array<func_t*, 0> funcs{};
};

template <typename R, typename... Ts>
  requires(
      std::same_as<
          typename xml_io_stream_functional<std::function<R(Ts...)>>::structs_t,
          void> or
      arts_xml_ioable<typename xml_io_stream_functional<
          std::function<R(Ts...)>>::structs_t>)
struct xml_io_stream<std::function<R(Ts...)>> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<std::function<R(Ts...)>>;

  using func_helper_t = xml_io_stream_functional<std::function<R(Ts...)>>;
  using func_t        = typename func_helper_t::func_t;
  using structs_t     = typename func_helper_t::structs_t;

  template <Size I = 0>
  static structs_t write_struct_helper(const std::function<R(Ts...)>& f
                                       [[maybe_unused]]) {
    if constexpr (I < std::variant_size_v<structs_t>) {
      using T = std::variant_alternative_t<I, structs_t>;

      if (const T* ptr = f.template target<T>(); ptr != nullptr) return *ptr;

      return write_struct_helper<I + 1>(f);
    } else {
      throw std::runtime_error("Unknown function struct");
    }
  }

  static void write(std::ostream& os,
                    const std::function<R(Ts...)>& f,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    const func_t* ptr = f.template target<func_t>();

    std::println(os,
                 R"(<{} name="{}" empty="{}" type="{}">)",
                 type_name,
                 name,
                 Index{not f},
                 ptr ? "array"sv : "struct"sv);

    if (f) {
      if (ptr != nullptr) {
        bool any = false;
        for (Size i = 0; i < func_helper_t::funcs.size(); i++) {
          if (ptr == reinterpret_cast<const func_t*>(func_helper_t::funcs[i])) {
            xml_write_to_stream(os, i, pbofs, "ID");
            any = true;
            break;
          }
        }

        if (not any) throw std::runtime_error("Unknown function pointer");
      } else {
        if constexpr (not std::same_as<void, structs_t>) {
          xml_write_to_stream(os, write_struct_helper(f), pbofs, "struct");
        } else {
          throw std::runtime_error("Cannot write function structs");
        }
      }
    }

    std::println(os, R"(</{}>)", type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Error writing {} (with mangled-name: '{}')\n{}",
                    type_name,
                    f.target_type().name(),
                    e.what()));
  }

  static void read(std::istream& is,
                   std::function<R(Ts...)>& f,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);

    Size empty;
    tag.get_attribute_value("empty", empty);

    if (empty) {
      f = std::function<R(Ts...)>{};
    } else {
      std::string type;
      tag.get_attribute_value("type", type);

      if (type == "array"sv) {
        Size id;
        xml_read_from_stream(is, id, pbifs);
        f = *func_helper_t::funcs.at(id);
      } else if (type == "struct"sv) {
        if constexpr (not std::same_as<void, structs_t>) {
          structs_t x{};
          xml_read_from_stream(is, x, pbifs);
          f = std::visit([](auto& v) -> std::function<R(Ts...)> { return v; },
                         x);
        } else {
          throw std::runtime_error("Cannot read function structs");
        }
      } else {
        throw std::runtime_error(std::format("Unknown type tag: {}", type));
      }
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Error reading {}\n{}", type_name, e.what()));
  }
};
