#include "optproperties.h"

std::string_view PTypeToString(const PType ptype) {
  switch (ptype) {
    case PTYPE_GENERAL:     return "general"sv; break;
    case PTYPE_TOTAL_RND:   return "totally_random"sv; break;
    case PTYPE_AZIMUTH_RND: return "azimuthally_random"sv; break;
    default:
      ARTS_USER_ERROR(
          "Internal error: Cannot map PType enum value {} to String.", ptype)
      break;
  }

  return "unknown PType"sv;
}
PType PType2FromString(const std::string_view ptype_string) {
  PType ptype;
  if (ptype_string == "general"sv)
    ptype = PTYPE_GENERAL;
  else if (ptype_string == "macroscopically_isotropic"sv)
    ptype = PTYPE_TOTAL_RND;
  else if (ptype_string == "horizontally_aligned"sv)
    ptype = PTYPE_AZIMUTH_RND;
  else {
    ARTS_USER_ERROR(
        "Unknown ptype: {}"
        "\n"
        "Valid types are: general, macroscopically_isotropic and "
        "horizontally_aligned.",
        ptype_string)
  }

  return ptype;
}

PType PTypeFromString(const std::string_view ptype_string) {
  PType ptype;
  if (ptype_string == "general"sv)
    ptype = PTYPE_GENERAL;
  else if (ptype_string == "totally_random"sv)
    ptype = PTYPE_TOTAL_RND;
  else if (ptype_string == "azimuthally_random"sv)
    ptype = PTYPE_AZIMUTH_RND;
  else {
    ARTS_USER_ERROR(
        "Unknown ptype: {}"
        "\n"
        "Valid types are: general, totally_random and "
        "azimuthally_random.",
        ptype_string)
  }

  return ptype;
}