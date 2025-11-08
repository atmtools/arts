#include "workspace_variable_shortnames.h"

namespace {
std::unordered_map<std::string, WsvShortForm>
workspace_variables_shortnames_impl() {
  std::unordered_map<std::string, WsvShortForm> shortnames;

  shortnames["atm"] = {
      .related_to = "The atmosphere",
      .desc       = "Used for variables related to the atmosphere",
  };

  shortnames["surf"] = {
      .related_to = "The surface",
      .desc       = "Used for variables related to the surface",
  };

  shortnames["subsurf"] = {
      .related_to = "The subsurface",
      .desc       = "Used for variables related to the subsurface",
  };

  shortnames["abs"] = {
      .related_to = "Absorption processes",
      .desc       = "Used for variables related to absorption models",
  };

  shortnames["rad"] = {
      .related_to = "Radiative transfer",
      .desc       = "Used for variables related to spectral radiance",
  };

  shortnames["propmat"] = {
      .related_to = "The propagation matrix",
      .desc       = "Used for variables related to the propagation matrix",
  };

  shortnames["tramat"] = {
      .related_to = "The transmittance matrix",
      .desc       = "Used for variables related to the transmittance matrix",
  };

  shortnames["phamat"] = {
      .related_to = "The phase matrix",
      .desc       = "Used for variables related to the phase matrix",
  };

  shortnames["srcvec"] = {
      .related_to = "The source vector",
      .desc       = "Used for variables related to the source vector",
  };

  shortnames["absvec"] = {
      .related_to = "The absorption vector",
      .desc       = "Used for variables related to the absorption vector",
  };

  shortnames["scat"] = {
      .related_to = "Scattering",
      .desc       = "Used for variables related to scattering",
  };

  shortnames["refl"] = {
      .related_to = "Reflectance",
      .desc       = "Used for variables related to reflectance",
  };

  shortnames["jac"] = {
      .related_to = "The Jacobian matrix",
      .desc       = "Used for variables related to the Jacobian matrix",
  };

  shortnames["covmat"] = {
      .related_to = "Covariance matrices",
      .desc       = "Used for variables related to covariance matrices",
  };

  shortnames["leg"] = {
      .related_to = "Legendre polynomials",
      .desc       = "Used for variables related to Legendre polynomials",
  };

  shortnames["nlte"] = {
      .related_to = "Non-local thermodynamic equilibrium",
      .desc =
          "Used for variables related to non-local thermodynamic equilibrium",
  };

  shortnames["cia"] = {
      .related_to = "Collision-induced absorption",
      .desc = "Used for variables related to collision-induced absorption",
  };

  shortnames["ecs"] = {
      .related_to = "Error-corrected sudden",
      .desc       = "Used for variables related to error-corrected sudden",
  };

  shortnames["xfit"] = {
      .related_to = "Cross-section fitting",
      .desc       = "Used for variables related to cross-section fitting",
  };

  shortnames["meas"] = {
      .related_to = "Measurements",
      .desc       = "Used for variables related to measurements",
  };

  shortnames["alt"] = {
      .related_to = "Altitude",
      .desc       = "Used for variables related to altitude position",
  };

  shortnames["lon"] = {
      .related_to = "Longitude",
      .desc       = "Used for variables related to longitude position",
  };

  shortnames["lat"] = {
      .related_to = "Latitude",
      .desc       = "Used for variables related to latitude position",
  };

  shortnames["freq"] = {
      .related_to = "Frequency",
      .desc       = "Used for variables related to a frequency",
  };

  shortnames["za"] = {
      .related_to = "Zenith angle",
      .desc       = "Used for variables related to zenith angle pointing",
  };

  shortnames["az"] = {
      .related_to = "Azimuth angle",
      .desc       = "Used for variables related to azimuth angle pointing",
  };

  shortnames["ray"] = {
      .related_to = "Ray-tracing",
      .desc       = "Used for variables related to ray-tracing",
  };

  shortnames["bkg"] = {
      .related_to = "Background",
      .desc       = "Used for variables related to background",
  };

  shortnames["obs"] = {
      .related_to = "Observers",
      .desc       = "Used for variables related to observers",
  };

  shortnames["pos"] = {
      .related_to = "Position",
      .desc       = "Used for variables related to position",
  };

  shortnames["los"] = {
      .related_to = "Line-of-sight",
      .desc       = "Used for variables related to line-of-sight.",
  };

  shortnames["grid"] = {
      .related_to = "1D Grids",
      .desc =
          "Used for variables related to one-dimensional grids.  These are always sorted and may act as coordinate axes for gridded variables.",
      .style = "Gridded Variables",
  };

  shortnames["field"] = {
      .related_to = "Multi-dimensional field",
      .desc = "Used for variables that can be sampled along one or more grids.",
      .style = "Gridded Variables",
  };

  shortnames["point"] = {
      .related_to = "Samples of a field",
      .desc =
          "Used for variables that are samples of a field at specific grid points.",
      .style = "Gridded Variables",
  };

  shortnames["profile"] = {
      .related_to = "Profiles",
      .desc =
          "Profiles are cuts through a field along a single dimension, e.g., the altitude.",
      .style = "Gridded Variables",
  };

  shortnames["path"] = {
      .related_to = "Values along a trajectory",
      .desc =
          "Paths contain values sampled along a trajectory through the atmosphere.  Unlike, paths are possibly not aligned with any grid.",
      .style = "Gridded Variables",
  };

  shortnames["spec"] = {
      .related_to = "Spectral variables",
      .desc  = "Used for variables that are defined over a frequency grid.",
      .style = "Gridded Variables",
  };

  shortnames["single"] = {
      .related_to = "Singular values on a spectral grid",
      .desc =
          "Used for variables that contain single values defined over a frequency grid.",
      .style = "Gridded Variables",
  };

  return shortnames;
}
}  // namespace

const std::unordered_map<std::string, WsvShortForm>&
workspace_variables_shortnames() {
  static const auto shortnames = workspace_variables_shortnames_impl();
  return shortnames;
}
