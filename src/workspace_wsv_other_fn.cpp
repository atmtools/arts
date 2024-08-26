#include "auto_wsg.h"

std::size_t Wsv::index() const { return value.index(); }

bool Wsv::holds_same(const Wsv& other) const {
  return index() == other.index();
}

const std::string_view Wsv::type_name() const { return name_wsg(value); }
