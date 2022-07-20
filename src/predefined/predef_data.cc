#include "predef_data.h"

#include <iostream>
#include <memory>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include "debug.h"
#include "double_imanip.h"

namespace Absorption::PredefinedModel {

Model::DataHolder construct_empty(PredefinedModelDataKey key) {
  switch (key) {
    case PredefinedModelDataKey::HITRANMTCKDWATER:
      return Hitran::MTCKD::WaterData{};
    case PredefinedModelDataKey::FINAL: { /* do nothing */
    }
  }
  return std::monostate{};
}

Model::Model(const Model& d) {
  for (auto& x : d.data) {
    std::visit([&](auto& v) { data[x.first] = v; }, x.second);
  }
}

Model& Model::operator=(const Model& d) {
  data.clear();

  for (auto& x : d.data) {
    std::visit([&](auto& v) { data[x.first] = v; }, x.second);
  }
  return *this;
}

namespace Hitran::MTCKD {
std::ostream& operator<<(std::ostream& os,
                         const Hitran::MTCKD::WaterData& data) {
  for (auto& x : data.for_absco_ref) os << x << ' ';
  os << '\n';
  for (auto& x : data.self_absco_ref) os << x << ' ';
  os << '\n';
  for (auto& x : data.wavenumbers) os << x << ' ';
  os << '\n';
  for (auto& x : data.self_texp) os << x << ' ';
  return os;
}

std::istream& operator>>(std::istream& is, WaterData& data) {
  for (auto& x : data.for_absco_ref) is >> double_imanip() >> x;
  for (auto& x : data.self_absco_ref) is >> double_imanip() >> x;
  for (auto& x : data.wavenumbers) is >> double_imanip() >> x;
  for (auto& x : data.self_texp) is >> double_imanip() >> x;
  return is;
}

void WaterData::resize(const std::vector<std::size_t>& inds) {
  ARTS_USER_ERROR_IF(inds.size() not_eq 1, "Expects only one size")
  self_absco_ref.resize(inds.front());
  for_absco_ref.resize(inds.front());
  wavenumbers.resize(inds.front());
  self_texp.resize(inds.front());
}
}  // namespace Hitran::MTCKD

void Model::set_all(const PredefinedModelData& d) {
  for (auto& x : d.data) {
    std::visit([&](auto& v) { data[x.first] = v; }, x.second);
  }
}

std::vector<PredefinedModelDataKey> Model::keys() const {
  std::vector<PredefinedModelDataKey> out;
  out.reserve(data.size());
  for (auto& x : data) out.emplace_back(x.first);
  return out;
}

struct SizesInterface {
  template <typename T>
  std::vector<std::size_t> operator()(const T& x) const {
    return x.sizes();
  }

  std::vector<std::size_t> operator()(const std::monostate&) const {
    ARTS_USER_ERROR("The data is corrupt")
  }
};

[[nodiscard]] std::vector<std::size_t> Model::data_size(
    PredefinedModelDataKey key) const {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  return std::visit(SizesInterface{}, data.at(key));
}

struct ReizeInterface {
  const std::vector<std::size_t>& inds;

  template <typename T>
  void operator()(T& x) const {
    x.resize(inds);
  }

  void operator()(std::monostate&) const {
    ARTS_USER_ERROR("The data is corrupt")
  }
};

void Model::resize(const std::vector<std::size_t>& inds,
                   PredefinedModelDataKey key) {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  if (data.find(key) == data.end()) data[key] = construct_empty(key);

  std::visit(ReizeInterface{inds}, data.at(key));
}

struct OutputInterface {
  std::ostream& os;

  template <typename T>
  void operator()(const T& x) {
    os << x;
  }

  void operator()(const std::monostate&) {
    ARTS_USER_ERROR("The data is corrupt")
  }
};

void Model::output_data_to_stream(std::ostream& os,
                                  PredefinedModelDataKey key) const {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  std::visit(OutputInterface{os}, data.at(key));
}

struct InputInterface {
  std::istream& is;

  template <typename T>
  void operator()(T& x) {
    is >> x;
  }

  void operator()(std::monostate&) {
    ARTS_USER_ERROR("The data is corrupt")
  }
};

void Model::set_data_from_stream(std::istream& is, PredefinedModelDataKey key) {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  std::visit(InputInterface{is}, data.at(key));
}

std::ostream& operator<<(std::ostream& os, const Model& data) {
  for (auto& x : data.data) {
    os << x.first << '\n';
    data.output_data_to_stream(os, x.first);
    os << '\n';
  }
  return os;
}
}  // namespace Absorption::PredefinedModel
