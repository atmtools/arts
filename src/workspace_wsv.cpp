#include "workspace_wsv.h"

#include "workspace_wsvwrapper.h"

Wsv::Wsv() : data{new WsvValueWrapper{}} {}

Wsv::Wsv(const Wsv& x) : data{new WsvValueWrapper{*x.data}} {}

Wsv::Wsv(Wsv&& x) noexcept : data(x.data) {x.data=nullptr;}

Wsv& Wsv::operator=(const Wsv& x) {
  *data = *x.data;
  return *this;
}

Wsv& Wsv::operator=(Wsv&& x) noexcept {
  data = x.data;
  x.data = nullptr;
  return *this;
}

Wsv::~Wsv() {
  if (data) data->~WsvValueWrapper();
  data = nullptr;
}

Wsv Wsv::copy() const {
  return std::visit([](const auto& val) -> Wsv { return Wsv{*val}; }, *data);
}

std::size_t Wsv::index() const { return data->index(); }

bool Wsv::holds_same(const Wsv& other) const {
  return index() == other.index();
}

const WsvValue& Wsv::value() const { return *data; }

WsvValue& Wsv::value() { return *data; }
