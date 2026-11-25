#pragma once

#include <artstime.h>
#include <enumsFieldComponent.h>
#include <matpack.h>
#include <xml.h>

namespace Atm {
struct MagnitudeField {
  GeodeticField3 magnitude{};
  GeodeticField3 theta{};
  GeodeticField3 phi{};

  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;

  [[nodiscard]] ConstVectorView                        x() const;
  [[nodiscard]] VectorView                             x();
  [[nodiscard]] std::vector<std::pair<Index, Numeric>> w(Numeric alt, Numeric lat, Numeric lon) const;
};

//! The eastward speed [m/s] of a rigidly rotating planet: its angular speed times the distance from its axis
struct RigidPlanetRotation {
  //! The rotation period [s], negative for a retrograde rotation
  Numeric period{};

  //! The ellipsoid of the planet [m]
  Vector2 ell{};

  Numeric operator()(Numeric alt, Numeric lat, Numeric lon) const;
};
}  // namespace Atm

template <> struct xml_io_stream_name<Atm::MagnitudeField> {
  static constexpr std::string_view name = "MagnitudeField";
};

template <> struct xml_io_stream_aggregate<Atm::MagnitudeField> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_name<Atm::RigidPlanetRotation> {
  static constexpr std::string_view name = "RigidPlanetRotation";
};

template <> struct xml_io_stream_aggregate<Atm::RigidPlanetRotation> {
  static constexpr bool value = true;
};
