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
      return std::make_unique<Hitran::MTCKD::WaterData>();
    case PredefinedModelDataKey::FINAL: { /* do nothing */
    }
  }
  return std::unique_ptr<Hitran::MTCKD::WaterData>(nullptr);
}

Model::Model(const Model& d) {
  for (auto& x : d.data) {
    std::visit(
        [&](auto&& v) {
          using MyType = std::remove_reference_t<decltype(*v)>;
          MyType y{*v};
          data[x.first] = std::make_unique<MyType>(std::move(y));
        },
        x.second);
  }
}

Model& Model::operator=(const Model& d) {
  data.clear();

  for (auto& x : d.data) {
    std::visit(
        [&](auto&& v) {
          using MyType = std::remove_reference_t<decltype(*v)>;
          MyType y{*v};
          data[x.first] = std::make_unique<MyType>(std::move(y));
        },
        x.second);
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
    std::visit(
        [&](auto&& v) {
          using MyType = std::remove_reference_t<decltype(*v)>;
          MyType y{*v};
          data[x.first] = std::make_unique<MyType>(std::move(y));
        },
        x.second);
  }
}

std::vector<PredefinedModelDataKey> Model::keys() const {
  std::vector<PredefinedModelDataKey> out;
  out.reserve(data.size());
  for (auto& x : data) out.emplace_back(x.first);
  return out;
}

[[nodiscard]] std::vector<std::size_t> Model::data_size(
    PredefinedModelDataKey key) const {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  return std::visit([&](auto& v) { return v->sizes(); }, data.at(key));
}

void Model::resize(const std::vector<std::size_t>& inds,
                   PredefinedModelDataKey key) {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  if (data.find(key) == data.end()) data[key] = construct_empty(key);

  std::visit([&](auto& v) { v->resize(inds); }, data.at(key));
}

void Model::output_data_to_stream(std::ostream& os,
                                  PredefinedModelDataKey key) const {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  std::visit([&](auto& v) { os << (*v); }, data.at(key));
}

void Model::set_data_from_stream(std::istream& is, PredefinedModelDataKey key) {
  ARTS_USER_ERROR_IF(not good_enum(key), "Bad key")
  std::visit([&](auto& v) { is >> (*v); }, data.at(key));
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