#include "predef_data.h"

#include <debug.h>
#include <double_imanip.h>

#include <iomanip>
#include <istream>
#include <ostream>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace Absorption::PredefinedModel {
Model::Model()                             = default;
Model::Model(const Model &)                = default;
Model::Model(Model &&) noexcept            = default;
Model &Model::operator=(const Model &)     = default;
Model &Model::operator=(Model &&) noexcept = default;

namespace MT_CKD400 {
std::ostream &operator<<(std::ostream &os, const WaterData &data) {
  os << data.ref_temp << ' ' << data.ref_press << ' ' << data.ref_h2o_vmr
     << '\n';
  for (auto &x : data.for_absco_ref) os << x << ' ';
  os << '\n';
  for (auto &x : data.self_absco_ref) os << x << ' ';
  os << '\n';
  for (auto &x : data.wavenumbers) os << x << ' ';
  os << '\n';
  for (auto &x : data.self_texp) os << x << ' ';
  return os;
}

std::istream &operator>>(std::istream &is, WaterData &data) {
  is >> double_imanip() >> data.ref_temp >> data.ref_press >> data.ref_h2o_vmr;
  for (auto &x : data.for_absco_ref) is >> double_imanip() >> x;
  for (auto &x : data.self_absco_ref) is >> double_imanip() >> x;
  for (auto &x : data.wavenumbers) is >> double_imanip() >> x;
  for (auto &x : data.self_texp) is >> double_imanip() >> x;
  return is;
}

void WaterData::resize(const std::vector<std::size_t> &inds) {
  ARTS_USER_ERROR_IF(inds.size() not_eq 1, "Expects only one size")
  self_absco_ref.resize(inds.front());
  for_absco_ref.resize(inds.front());
  wavenumbers.resize(inds.front());
  self_texp.resize(inds.front());
}

std::vector<std::size_t> WaterData::sizes() const {
  return {static_cast<std::size_t>(self_absco_ref.size())};
};
}  // namespace MT_CKD400

std::ostream &operator<<(std::ostream &os, const ModelName &) { return os; }

std::istream &operator>>(std::istream &is, ModelName &) { return is; }

std::ostream &operator<<(std::ostream &os, const Model &m) {
  std::string_view nline = "";

  for (auto &[tag, data] : m.data) {
    os << std::exchange(nline, "\n") << tag << ':' << '\n';
    std::visit([&os](auto &&arg) { os << arg; }, data);
  }

  return os;
}

std::string_view model_name(const ModelVariant &data) {
  if (std::holds_alternative<ModelName>(data)) {
    return "ModelName";
  }

  if (std::holds_alternative<MT_CKD400::WaterData>(data)) {
    return "MT_CKD400::WaterData";
  }

  throw std::runtime_error("Unspecificed model type, this is a developer bug");
}

ModelVariant model_data(const std::string_view name) {
  if (name == "ModelName") {
    return ModelName{};
  }

  if (name == "MT_CKD400::WaterData") {
    return MT_CKD400::WaterData{};
  }

  throw std::runtime_error(std::format(
      R"(Unknown model name: "{}". Are all models defined?)", name));
}

}  // namespace Absorption::PredefinedModel
