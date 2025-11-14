#include "workspace_variable_shortnames.h"

namespace {
std::unordered_map<std::string, WsvShortForm>
workspace_variables_shortnames_impl() {
  std::unordered_map<std::string, WsvShortForm> shortnames;

  shortnames["atm"]     = {.desc = "The atmosphere"};
  shortnames["surf"]    = {.desc = "The surface"};
  shortnames["subsurf"] = {.desc = "The subsurface"};
  shortnames["abs"]     = {.desc = "Absorption processes"};
  shortnames["rad"]     = {.desc = "Radiative transfer/radiance"};
  shortnames["propmat"] = {.desc = "The propagation matrix"};
  shortnames["tramat"]  = {.desc = "The transmittance matrix"};
  shortnames["phamat"]  = {.desc = "The phase matrix"};
  shortnames["srcvec"]  = {.desc = "The source vector"};
  shortnames["absvec"]  = {.desc = "The absorption vector"};
  shortnames["scat"]    = {.desc = "Scattering"};
  shortnames["refl"]    = {.desc = "Reflectance"};
  shortnames["jac"]     = {.desc = "The Jacobian matrix"};
  shortnames["covmat"]  = {.desc = "Covariance matrices"};
  shortnames["leg"]     = {.desc = "Legendre polynomials"};
  shortnames["nlte"]    = {.desc = "Non-local thermodynamic equilibrium"};
  shortnames["cia"]     = {.desc = "Collision-induced absorption"};
  shortnames["ecs"]     = {.desc = "Error-corrected sudden"};
  shortnames["xfit"]    = {.desc = "Cross-section fitting"};
  shortnames["meas"]    = {.desc = "Measurements"};
  shortnames["ray"]     = {.desc = "Ray-tracing"};
  shortnames["bkg"]     = {.desc = "Background"};
  shortnames["obs"]     = {.desc = "Observers"};
  shortnames["pos"]     = {.desc = "Position"};
  shortnames["los"]     = {.desc = "Line-of-sight"};
  shortnames["dim"]     = {.desc = "Dimension"};
  shortnames["quad"]    = {.desc = "Quadrature"};

  shortnames["alt"] = {
      .desc  = "Altitude",
      .title = "Coordinates",
  };

  shortnames["lon"] = {
      .desc  = "Longitude [-180, 180) circular.",
      .title = "Coordinates",
  };

  shortnames["lat"] = {
      .desc  = "Latitude [-90, 90].",
      .title = "Coordinates",
  };

  shortnames["freq"] = {
      .desc  = "Frequency",
      .title = "Coordinates",
  };

  shortnames["za"] = {
      .desc  = "Zenith angle [0, 180].",
      .title = "Coordinates",
  };

  shortnames["az"] = {
      .desc  = "Azimuth angle [0, 360) circular.  0° is North, 90° is East.",
      .title = "Coordinates",
  };

  shortnames["grid"] = {
      .desc  = "Sorted one-dimensional grids.",
      .title = "Gridded",
  };

  shortnames["field"] = {
      .desc =
          "Contiguous values defined along one or more dimensions.  May be ungridded.",
      .title = "Gridded",
  };

  shortnames["point"] = {
      .desc  = "Values are extracted at a single coordinate in a field.",
      .title = "Gridded",
  };

  shortnames["profile"] = {
      .desc =
          "Profiles are cuts through a field along a single dimension, e.g., the altitude.",
      .title = "Gridded",
  };

  shortnames["path"] = {
      .desc =
          "Paths contain values sampled along a trajectory through the atmosphere.  Unlike profiles, paths are possibly not aligned with any grid.",
      .title = "Gridded",
  };

  shortnames["spec"] = {
      .desc  = "Defined over a frequency grid.",
      .title = "Gridded",
  };

  shortnames["single"] = {
      .desc  = "Contains a single value extracted from a frequency grid.",
      .title = "Gridded",
  };

  return shortnames;
}
}  // namespace

const std::unordered_map<std::string, WsvShortForm>&
workspace_variables_shortnames() {
  static const auto shortnames = workspace_variables_shortnames_impl();
  return shortnames;
}
