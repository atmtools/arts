#include "gravity.h"

#include "arts_constexpr_math.h"
#include "planet_data.h"

Numeric Gravity::operator()(Numeric h, Numeric lat, Numeric lon) const {
  using Conversion::cosd;
  using Conversion::sind;
  using Math::pow2;
  using std::sqrt;

  const Numeric N  = a / sqrt(1 - pow2(e * sind(lat)));
  const Numeric r2 = pow2((N + h) * cosd(lon) * cosd(lat)) +
                     pow2((N + h) * sind(lon) * cosd(lat)) +
                     pow2((N * (1 - pow2(e)) + h) * sind(lat));

  return GM / r2;
}

Gravity Gravity::Earth() {
  return {.GM = Body::Earth::GM,
          .a  = Body::Earth::a,
          .e  = std::sqrt(1 - Math::pow2(Body::Earth::b / Body::Earth::a))};
}

Gravity Gravity::Mars() {
  return {.GM = Body::Mars::GM,
          .a  = Body::Mars::a,
          .e  = std::sqrt(1 - Math::pow2(Body::Mars::b / Body::Mars::a))};
}

Gravity Gravity::Venus() {
  return {.GM = Body::Venus::GM,
          .a  = Body::Venus::a,
          .e  = std::sqrt(1 - Math::pow2(Body::Venus::b / Body::Venus::a))};
}

Gravity Gravity::Jupiter() {
  return {.GM = Body::Jupiter::GM,
          .a  = Body::Jupiter::a,
          .e  = std::sqrt(1 - Math::pow2(Body::Jupiter::b / Body::Jupiter::a))};
}

Gravity Gravity::Moon() {
  return {.GM = Body::Moon::GM,
          .a  = Body::Moon::a,
          .e  = std::sqrt(1 - Math::pow2(Body::Moon::b / Body::Moon::a))};
}

Gravity Gravity::Mercury() {
  return {.GM = Body::Mercury::GM,
          .a  = Body::Mercury::a,
          .e  = std::sqrt(1 - Math::pow2(Body::Mercury::b / Body::Mercury::a))};
}

Gravity Gravity::Saturn() {
  return {.GM = Body::Saturn::GM,
          .a  = Body::Saturn::a,
          .e  = std::sqrt(1 - Math::pow2(Body::Saturn::b / Body::Saturn::a))};
}
