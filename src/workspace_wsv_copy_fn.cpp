#include "auto_wsg.h"

Wsv Wsv::copy() const {
  return std::visit([](const auto& val) -> Wsv {
    return Wsv{*val};
  }, value);
}
