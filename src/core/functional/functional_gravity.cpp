#include "functional_gravity.h"

#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <geodetic.h>
#include <planet_data.h>

Numeric EllipsoidGravity::operator()(Numeric h,
                                     Numeric lat,
                                     Numeric lon) const {
  const auto [r, lat2, lon2] = geodetic2geocentric({h, lat, lon}, {a, b});

  return GM / Math::pow2(r);
}

EllipsoidGravity EllipsoidGravity::Earth() {
  return {.GM = Body::Earth::GM, .a = Body::Earth::a, .b = Body::Earth::b};
}

EllipsoidGravity EllipsoidGravity::Mars() {
  return {.GM = Body::Mars::GM, .a = Body::Mars::a, .b = Body::Mars::b};
}

EllipsoidGravity EllipsoidGravity::Venus() {
  return {.GM = Body::Venus::GM, .a = Body::Venus::a, .b = Body::Venus::b};
}

EllipsoidGravity EllipsoidGravity::Jupiter() {
  return {
      .GM = Body::Jupiter::GM, .a = Body::Jupiter::a, .b = Body::Jupiter::b};
}

EllipsoidGravity EllipsoidGravity::Moon() {
  return {.GM = Body::Moon::GM, .a = Body::Moon::a, .b = Body::Moon::b};
}

EllipsoidGravity EllipsoidGravity::Mercury() {
  return {
      .GM = Body::Mercury::GM, .a = Body::Mercury::a, .b = Body::Mercury::b};
}

EllipsoidGravity EllipsoidGravity::Saturn() {
  return {.GM = Body::Saturn::GM, .a = Body::Saturn::a, .b = Body::Saturn::b};
}
