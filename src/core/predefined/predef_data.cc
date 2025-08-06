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
namespace MT_CKD400 {
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

std::string_view model_name(const ModelVariant &data) {
  if (std::holds_alternative<ModelName>(data.data)) {
    return "ModelName";
  }

  if (std::holds_alternative<MT_CKD400::WaterData>(data.data)) {
    return "MT_CKD400::WaterData";
  }

  throw std::runtime_error("Unspecificed model type, this is a developer bug");
}

ModelVariant model_data(const std::string_view name) {
  if (name == "ModelName") {
    return ModelVariant{.data = ModelName{}};
  }

  if (name == "MT_CKD400::WaterData") {
    return ModelVariant{.data = MT_CKD400::WaterData{}};
  }

  throw std::runtime_error(std::format(
      R"(Unknown model name: "{}". Are all models defined?)", name));
}
}  // namespace Absorption::PredefinedModel

void xml_io_stream<Absorption::PredefinedModel::MT_CKD400::WaterData>::write(
    std::ostream &os,
    const Absorption::PredefinedModel::MT_CKD400::WaterData &x,
    bofstream *pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.ref_press, pbofs);
  xml_write_to_stream(os, x.ref_temp, pbofs);
  xml_write_to_stream(os, x.ref_h2o_vmr, pbofs);
  xml_write_to_stream(os, x.self_absco_ref, pbofs);
  xml_write_to_stream(os, x.for_absco_ref, pbofs);
  xml_write_to_stream(os, x.wavenumbers, pbofs);
  xml_write_to_stream(os, x.self_texp, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Absorption::PredefinedModel::MT_CKD400::WaterData>::read(
    std::istream &is,
    Absorption::PredefinedModel::MT_CKD400::WaterData &x,
    bifstream *pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.ref_press, pbifs);
  xml_read_from_stream(is, x.ref_temp, pbifs);
  xml_read_from_stream(is, x.ref_h2o_vmr, pbifs);
  xml_read_from_stream(is, x.self_absco_ref, pbifs);
  xml_read_from_stream(is, x.for_absco_ref, pbifs);
  xml_read_from_stream(is, x.wavenumbers, pbifs);
  xml_read_from_stream(is, x.self_texp, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Absorption::PredefinedModel::ModelName>::write(
    std::ostream &os,
    const Absorption::PredefinedModel::ModelName &,
    bofstream *,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}

void xml_io_stream<Absorption::PredefinedModel::ModelName>::read(
    std::istream &is, Absorption::PredefinedModel::ModelName &, bifstream *) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<PredefinedModelDataVariant>::write(
    std::ostream &os,
    const PredefinedModelDataVariant &x,
    bofstream *pbofs,
    std::string_view name) {
  xml_write_to_stream(os, x.data, pbofs, name);
}

void xml_io_stream<PredefinedModelDataVariant>::read(
    std::istream &is, PredefinedModelDataVariant &x, bifstream *pbifs) {
  xml_read_from_stream(is, x.data, pbifs);
}
