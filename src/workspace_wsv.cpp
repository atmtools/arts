#include "workspace_wsv.h"

#include "workspace_wsvwrapper.h"

Wsv::Wsv() : data{new WsvValueWrapper{}} {}

Wsv::Wsv(const Wsv& x) : data{new WsvValueWrapper{*x.data}} {}

Wsv::Wsv(Wsv&& x) noexcept : data(x.data) { x.data = nullptr; }

Wsv& Wsv::operator=(const Wsv& x) {
  ARTS_ASSERT(data, "Wsv data has been deleted!");
  ARTS_ASSERT(x.data, "Other wsv data has already been deleted!");

  *data = *x.data;

  return *this;
}

Wsv& Wsv::operator=(Wsv&& x) noexcept {
  if (data != x.data) this->~Wsv();

  data = x.data;
  x.data = nullptr;

  return *this;
}

Wsv::~Wsv() {
  delete data;

  data = nullptr;
}

Wsv Wsv::copy() const {
  ARTS_ASSERT(data, "Wsv data has been deleted!");
  return std::visit([](const auto& val) -> Wsv { return Wsv{*val}; }, *data);
}

std::size_t Wsv::index() const {
  ARTS_ASSERT(data, "Wsv data has been deleted!");
  return data->index();
}

bool Wsv::holds_same(const Wsv& x) const {
  ARTS_ASSERT(data, "Wsv data has been deleted!");
  return index() == x.index();
}

const WsvValue& Wsv::value() const {
  ARTS_ASSERT(data, "Wsv data has been deleted!");
  return *data;
}

WsvValue& Wsv::value() {
  ARTS_ASSERT(data, "Wsv data has been deleted!");
  return *data;
}
