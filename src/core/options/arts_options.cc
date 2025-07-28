#include "arts_options.h"

#include <mystring.h>

#include <algorithm>
#include <array>
#include <format>
#include <iterator>
#include <limits>
#include <print>
#include <ranges>
#include <sstream>
#include <utility>
#include <vector>

using namespace std::literals;

using Value = std::vector<std::string>;

namespace stdr = std::ranges;
namespace stdv = std::ranges::views;

namespace {
void fix_static(std::vector<EnumeratedOption>& opts) {
  for (auto& opt : opts) {
    if (opt.desc.empty()) throw std::runtime_error("No empty docstrings");

    if (opt.desc.back() != '\n') opt.desc.push_back('\n');

    for (auto& v : opt.values_and_desc) {
      if (v.size() < 2) throw std::runtime_error("No empty value descriptions");

      auto& l = v.back();

      if (l.empty()) throw std::runtime_error("No empty value descriptions");

      if (l.back() != '\n') v.back().push_back('\n');
    }
  }
}

std::vector<EnumeratedOption> internal_options_create() {
  std::vector<EnumeratedOption> opts;

  opts.emplace_back(EnumeratedOption{
      .name = "ray_path_observer_agendaSetGeometricMaxStep",
      .desc =
          R"(For use with *ray_path_observer_agendaSetGeometric*.  Determines how to densify the geometric path.
)",
      .values_and_desc =
          {Value{"half", "1/2", "Use *ray_pathFillGeometricHalfStep*"},
           Value{"step", "linear", "Use *ray_pathFillGeometricStepwise*"},
           Value{"None", "0", "Do not use any step filling method"}},
      .preferred_print = 1,
  });

  opts.emplace_back(EnumeratedOption{
      .name = "disort_settings_agenda_setup_layer_emission_type",
      .desc =
          R"(For atmospheric emission settings with *disort_settings_agendaSetup*.
)",
      .values_and_desc = {
          Value{"None", "No atmospheric emission"},
          Value{"LinearInTau", "Layer emission is linear in optical depth"},
          Value{
              "LinearInTauNonLTE",
              "Layer emission is linear in optical depth taking non-LTE into account"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "disort_settings_agenda_setup_sun_type",
      .desc =
          R"(For direct radiation settings with *disort_settings_agendaSetup*.
)",
      .values_and_desc = {
          Value{"None", "No direct radiation"},
          Value{"Sun", "Add a *Sun* object as the direct radiation source"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "disort_settings_agenda_setup_space_type",
      .desc =
          R"(For space boundary settings with *disort_settings_agendaSetup*.
)",
      .values_and_desc = {
          Value{"None", "No isotropic emission from space"},
          Value{
              "CosmicMicrowaveBackgroundRadiation",
              "Isotropic emission as per cosmic microwave background radiation"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "disort_settings_agenda_setup_scattering_type",
      .desc =
          R"(For atmospheric scattering settings with *disort_settings_agendaSetup*.
)",
      .values_and_desc = {
          Value{"None", "No atmospheric scattering"},
          Value{"ScatteringSpecies",
                "Use *ArrayOfScatteringSpecies* for scattering"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "disort_settings_agenda_setup_surface_type",
      .desc =
          R"(For surface boundary settings with *disort_settings_agendaSetup*.
)",
      .values_and_desc = {
          Value{"None", "No surface emission or scattering"},
          Value{"Thermal", "Thermal emission, no scattering"},
          Value{"ThermalLambertian", "Thermal emission, Lambertian scattering"},
          Value{"Lambertian", "No surface emission, Lambertian scattering"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "SensorKeyType",
      .desc =
          R"(A key for identifying a sensor property
)",
      .values_and_desc = {
          Value{"f", "Frequency", "freq", "Frequency of the sensor"},
          Value{"za", "Zenith", "zenith", "Zenith angle of the sensor"},
          Value{"aa", "Azimuth", "azimuth", "Azimuth angle of the sensor"},
          Value{"alt", "Altitude", "altitude", "Altitude of the sensor"},
          Value{"lat", "Latitude", "latitude", "Latitude of the sensor"},
          Value{"lon", "Longitude", "longitude", "Longitude of the sensor"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "SensorJacobianModelType",
      .desc =
          R"(How to model the sensor Jacobian model target.
)",
      .values_and_desc = {
          Value{"None", "No model, work purely on sensor data"},
          Value{"PolynomialOffset",
                "The sensor Jacobian is modeled as a polynomial offset"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "HitranLineStrengthOption",
      .desc =
          R"(The way line strength is computed in ARTS when reading Hitran data.

ARTS uses Einstein A-coefficients to compute the line strength. Hitran provides
both the line strength and the Einstein A-coefficient. There is a 1-to-1 conversion
between these two.  However, as with all data, the numbers might differ slightly even
if good cases, so we have provide this selection mechanism to make them match.
)",
      .values_and_desc = {
          Value{"S", "s", "strenght", "s0", "S0", "Line strength"},
          Value{"A", "a", "einstein", "ein", "A0", "Absorption intensity"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "TimeStepType",
      .desc =
          "Time step for the integration.  Is often accompanied by a value.\n",
      .values_and_desc =
          {Value{"hour", "hours", "h", "Time step is in hours"},
           Value{"minute", "minutes", "min", "Time step is in minutes"},
           Value{"second", "seconds", "s", "Time step is in seconds"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "GridType",
      .desc            = R"(Type of Lagrange interpolation weights.
)",
      .values_and_desc = {
          Value{"Standard", "1-to-1 interpolation grid"},
          Value{"Cyclic", "Cyclic interpolation grid"},
          Value{"Log", "Natural logarithm interpolation grid"},
          Value{"Log10", "10-base logarithm interpolation grid"},
          Value{"Log2", "2-base logarithm interpolation grid"},
          Value{
              "SinDeg",
              "Sine in degrees interpolation grid, grid only defined [-90, 90]"},
          Value{
              "SinRad",
              "Sine in radians interpolation grid, grid only defined [-PI/2, PI/2]"},
          Value{
              "CosDeg",
              "Cosine in degrees interpolation grid, grid only defined [0, 180]"},
          Value{
              "CosRad",
              "Cosine in radians interpolation grid, grid only defined [0,  PI]"},
      }});

  opts.emplace_back(
      EnumeratedOption{.name = "SurfaceKey",
                       .desc = R"(A key to identify a surface property.
)",
                       .values_and_desc = {
                           Value{"h", "elevation", "Altitude [m]"},
                           Value{"t", "temperature", "Temperature [K]"},
                       }});

  opts.emplace_back(EnumeratedOption{
      .name            = "AtmKey",
      .desc            = R"(A key to identify an atmospheric property.
)",
      .values_and_desc = {
          Value{"t", "temperature", "Temperature [K]"},
          Value{"p", "pressure", "Pressure [Pa]"},
          Value{"wind_u", "WindU", "Wind field U-component [m/s]"},
          Value{"wind_v", "WindV", "Wind field V-component [m/s]"},
          Value{"wind_w", "WindW", "Wind field W-component [m/s]"},
          Value{"mag_u", "MagU", "Magnetic field U-component [T]"},
          Value{"mag_v", "MagV", "Magnetic field V-component [T]"},
          Value{"mag_w", "MagW", "Magnetic field W-component [T]"},
      }});

  opts.emplace_back(
      EnumeratedOption{.name = "SubsurfaceKey",
                       .desc = R"(A key to identify a subsurface property.
)",
                       .values_and_desc = {
                           Value{"t", "temperature", "Temperature [K]"},
                           Value{"rho", "density", "Density [kg/m^3]"},
                       }});

  opts.emplace_back(EnumeratedOption{
      .name = "InterpolationExtrapolation",
      .desc =
          R"(Instructions about how to handle extrapolation of interpolated data.
)",
      .values_and_desc = {
          Value{"None", "Do not allow interpolation of data outside the grid"},
          Value{"Zero", "Values outside the grid are set to zero"},
          Value{"Nearest",
                "Values outside the grid are set to the nearest grid point"},
          Value{"Linear", "Linearly extrapolate the data outside the grid"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "PartitionFunctionsType",
      .desc =
          R"(Type of partition function data.
)",
      .values_and_desc = {
          Value{"Interp", "Interpolate the data"},
          Value{"Coeff", "Use as polynomial coefficients"},
          Value{
              "StaticInterp",
              "Interpolate the data, the temperature grid is known at compile time."},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "SpeciesEnum",
      .desc =
          R"(The valid species for the ARTS system.
)",
      .values_and_desc =
          {
              Value{"Bath", "AIR", "Any non-descript species"},
              Value{"Water", "H2O", "Water molecule"},
              Value{"CarbonDioxide", "CO2", "Carbon Dioxide molecule"},
              Value{"Ozone", "O3", "Ozone molecule"},
              Value{"NitrogenOxide", "N2O", "Nitrogen Oxide molecule"},
              Value{"CarbonMonoxide", "CO", "Carbon Monoxide molecule"},
              Value{"Methane", "CH4", "Methane molecule"},
              Value{"Oxygen", "O2", "Oxygen molecule"},
              Value{"NitricOxide", "NO", "Nitric Oxide molecule"},
              Value{"SulfurDioxide", "SO2", "Sulfur Dioxide molecule"},
              Value{"NitrogenDioxide", "NO2", "Nitrogen Dioxide molecule"},
              Value{"Ammonia", "NH3", "Ammonia molecule"},
              Value{"NitricAcid", "HNO3", "Nitric Acid molecule"},
              Value{"Hydroxyl", "OH", "Hydroxyl molecule"},
              Value{"HydrogenFluoride", "HF", "Hydrogen Fluoride molecule"},
              Value{"HydrogenChloride", "HCl", "Hydrogen Chloride molecule"},
              Value{"HydrogenBromide", "HBr", "Hydrogen Bromide molecule"},
              Value{"HydrogenIodide", "HI", "Hydrogen Iodide molecule"},
              Value{"ChlorineMonoxide", "ClO", "Chlorine Monoxide molecule"},
              Value{"CarbonylSulfide", "OCS", "Carbonyl Sulfide molecule"},
              Value{"Formaldehyde", "H2CO", "Formaldehyde molecule"},
              Value{"HeavyFormaldehyde", "HDCO", "Heavy Formaldehyde molecule"},
              Value{"VeryHeavyFormaldehyde",
                    "D2CO",
                    "Very Heavy Formaldehyde molecule"},
              Value{"HypochlorousAcid", "HOCl", "Hypochlorous Acid molecule"},
              Value{"Nitrogen", "N2", "Nitrogen molecule"},
              Value{"HydrogenCyanide", "HCN", "Hydrogen Cyanide molecule"},
              Value{"Chloromethane", "CH3Cl", "Chloromethane molecule"},
              Value{"HydrogenPeroxide", "H2O2", "Hydrogen Peroxide molecule"},
              Value{"Acetylene", "C2H2", "Acetylene molecule"},
              Value{"Ethane", "C2H6", "Ethane molecule"},
              Value{"Phosphine", "PH3", "Phosphine molecule"},
              Value{"CarbonylFluoride", "COF2", "Carbonyl Fluoride molecule"},
              Value{
                  "SulfurHexafluoride", "SF6", "Sulfur Hexafluoride molecule"},
              Value{"HydrogenSulfide", "H2S", "Hydrogen Sulfide molecule"},
              Value{"FormicAcid", "HCOOH", "Formic Acid molecule"},
              Value{"LeftHeavyFormicAcid",
                    "DCOOH",
                    "Left Heavy Formic Acid molecule"},
              Value{"RightHeavyFormicAcid",
                    "HCOOD",
                    "Right Heavy Formic Acid molecule"},
              Value{"Hydroperoxyl", "HO2", "Hydroperoxyl molecule"},
              Value{"OxygenAtom", "O", "Oxygen atom"},
              Value{"ChlorineNitrate", "ClONO2", "Chlorine Nitrate molecule"},
              Value{"NitricOxideCation", "NO+", "Nitric Oxide Cation molecule"},
              Value{"HypobromousAcid", "HOBr", "Hypobromous Acid molecule"},
              Value{"Ethylene", "C2H4", "Ethylene molecule"},
              Value{"Methanol", "CH3OH", "Methanol molecule"},
              Value{"Bromomethane", "CH3Br", "Bromomethane molecule"},
              Value{"Acetonitrile", "CH3CN", "Acetonitrile molecule"},
              Value{
                  "HeavyAcetonitrile", "CH2DCN", "Heavy Acetonitrile molecule"},
              Value{"CarbonTetrafluoride",
                    "CF4",
                    "Carbon Tetrafluoride molecule"},
              Value{"Diacetylene", "C4H2", "Diacetylene molecule"},
              Value{"Cyanoacetylene", "HC3N", "Cyanoacetylene molecule"},
              Value{"Hydrogen", "H2", "Hydrogen molecule"},
              Value{"CarbonMonosulfide", "CS", "Carbon Monosulfide molecule"},
              Value{"SulfurTrioxide", "SO3", "Sulfur Trioxide molecule"},
              Value{"Cyanogen", "C2N2", "Cyanogen molecule"},
              Value{"Phosgene", "COCl2", "Phosgene molecule"},
              Value{"SulfurMonoxide", "SO", "Sulfur Monoxide molecule"},
              Value{"CarbonDisulfide", "CS2", "Carbon Disulfide molecule"},
              Value{"Methyl", "CH3", "Methyl molecule"},
              Value{"Cyclopropene", "C3H4", "Cyclopropene molecule"},
              Value{"SulfuricAcid", "H2SO4", "Sulfuric Acid molecule"},
              Value{
                  "HydrogenIsocyanide", "HNC", "Hydrogen Isocyanide molecule"},
              Value{"BromineMonoxide", "BrO", "Bromine Monoxide molecule"},
              Value{"ChlorineDioxide", "OClO", "Chlorine Dioxide molecule"},
              Value{"Propane", "C3H8", "Propane molecule"},
              Value{"Helium", "He", "Helium atom"},
              Value{"ChlorineMonoxideDimer",
                    "Cl2O2",
                    "Chlorine Monoxide Dimer molecule"},
              Value{"HydrogenAtom", "H", "Hydrogen atom"},
              Value{"Argon", "Ar", "Argon atom"},
              Value{"Hexafluoroethane", "C2F6", "Hexafluoroethane molecule"},
              Value{"Perfluoropropane", "C3F8", "Perfluoropropane molecule"},
              Value{"Perfluorobutane", "C4F10", "Perfluorobutane molecule"},
              Value{"Perfluoropentane", "C5F12", "Perfluoropentane molecule"},
              Value{"Perfluorohexane", "C6F14", "Perfluorohexane molecule"},
              Value{"Perfluorooctane", "C8F18", "Perfluorooctane molecule"},
              Value{"Perfluorocyclobutane",
                    "cC4F8",
                    "Perfluorocyclobutane molecule"},
              Value{"CarbonTetrachloride",
                    "CCl4",
                    "Carbon Tetrachloride molecule"},
              Value{"CFC11", "CFC11", "CFC11 molecule"},
              Value{"CFC113", "CFC113", "CFC113 molecule"},
              Value{"CFC114", "CFC114", "CFC114 molecule"},
              Value{"CFC115", "CFC115", "CFC115 molecule"},
              Value{"CFC12", "CFC12", "CFC12 molecule"},
              Value{"Dichloromethane", "CH2Cl2", "Dichloromethane molecule"},
              Value{"Trichloroethane", "CH3CCl3", "Trichloroethane molecule"},
              Value{"Trichloromethane", "CHCl3", "Trichloromethane molecule"},
              Value{"Bromochlorodifluoromethane",
                    "Halon1211",
                    "Bromochlorodifluoromethane molecule"},
              Value{"Bromotrifluoromethane",
                    "Halon1301",
                    "Bromotrifluoromethane molecule"},
              Value{"Dibromotetrafluoroethane",
                    "Halon2402",
                    "Dibromotetrafluoroethane molecule"},
              Value{"HCFC141b", "HCFC141b", "HCFC141b molecule"},
              Value{"HCFC142b", "HCFC142b", "HCFC142b molecule"},
              Value{"HCFC22", "HCFC22", "HCFC22 molecule"},
              Value{"HFC125", "HFC125", "HFC125 molecule"},
              Value{"HFC134a", "HFC134a", "HFC134a molecule"},
              Value{"HFC143a", "HFC143a", "HFC143a molecule"},
              Value{"HFC152a", "HFC152a", "HFC152a molecule"},
              Value{"HFC227ea", "HFC227ea", "HFC227ea molecule"},
              Value{"HFC23", "HFC23", "HFC23 molecule"},
              Value{"HFC236fa", "HFC236fa", "HFC236fa molecule"},
              Value{"HFC245fa", "HFC245fa", "HFC245fa molecule"},
              Value{"HFC32", "HFC32", "HFC32 molecule"},
              Value{"HFC365mfc", "HFC365mfc", "HFC365mfc molecule"},
              Value{"NitrogenTrifluoride",
                    "NF3",
                    "Nitrogen Trifluoride molecule"},
              Value{"SulfurylFluoride", "SO2F2", "Sulfuryl Fluoride molecule"},
              Value{"HFC4310mee", "HFC4310mee", "HFC4310mee molecule"},
              Value{"Germane", "GeH4", "Germane molecule"},
              Value{"Iodomethane", "CH3I", "Iodomethane molecule"},
              Value{"Fluoromethane", "CH3F", "Fluoromethane molecule"},
              Value{"liquidcloud", "liquidcloud", "liquidcloud tag"},
              Value{"icecloud", "icecloud", "icecloud tag"},
              Value{"rain", "rain", "rain tag"},
              Value{"free_electrons", "free_electrons", "free electrons tag"},
              Value{"particles", "particles", "particles tag"},
          },
      .preferred_print = 1,
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "QuantumNumberType",
      .desc            = R"(The type of value for a quantum number.
)",
      .values_and_desc = {
          Value{"alpha", "Quantum number \"alpha\" - Not in VAMDC"},
          Value{"config", "Quantum number \"config\" - Not in VAMDC"},
          Value{"ElecStateLabel", "Quantum number \"ElecStateLabel\""},
          Value{"F", "Quantum number \"F\""},
          Value{"F1", "Quantum number \"F1\""},
          Value{"F10", "Quantum number \"F10\""},
          Value{"F11", "Quantum number \"F11\""},
          Value{"F12", "Quantum number \"F12\""},
          Value{"F2", "Quantum number \"F2\""},
          Value{"F3", "Quantum number \"F3\""},
          Value{"F4", "Quantum number \"F4\""},
          Value{"F5", "Quantum number \"F5\""},
          Value{"F6", "Quantum number \"F6\""},
          Value{"F7", "Quantum number \"F7\""},
          Value{"F8", "Quantum number \"F8\""},
          Value{"F9", "Quantum number \"F9\""},
          Value{"I", "Quantum number \"I\""},
          Value{"J", "Quantum number \"J\""},
          Value{"K", "Quantum number \"K\""},
          Value{"Ka", "Quantum number \"Ka\""},
          Value{"Kc", "Quantum number \"Kc\""},
          Value{"L", "Quantum number \"L\" - Not in VAMDC"},
          Value{"Lambda", "Quantum number \"Lambda\""},
          Value{"N", "Quantum number \"N\""},
          Value{"Omega", "Quantum number \"Omega\""},
          Value{"S", "Quantum number \"S\""},
          Value{"Sigma", "Quantum number \"Sigma\""},
          Value{"SpinComponentLabel", "Quantum number \"SpinComponentLabel\""},
          Value{"asSym", "Quantum number \"asSym\""},
          Value{"elecInv", "Quantum number \"elecInv\""},
          Value{"elecRefl", "Quantum number \"elecRefl\""},
          Value{"elecSym", "Quantum number \"elecSym\""},
          Value{"kronigParity", "Quantum number \"kronigParity\""},
          Value{"l", "Quantum number \"l\""},
          Value{"l1", "Quantum number \"l1\""},
          Value{"l10", "Quantum number \"l10\""},
          Value{"l11", "Quantum number \"l11\""},
          Value{"l12", "Quantum number \"l12\""},
          Value{"l2", "Quantum number \"l2\""},
          Value{"l3", "Quantum number \"l3\""},
          Value{"l4", "Quantum number \"l4\""},
          Value{"l5", "Quantum number \"l5\""},
          Value{"l6", "Quantum number \"l6\""},
          Value{"l7", "Quantum number \"l7\""},
          Value{"l8", "Quantum number \"l8\""},
          Value{"l9", "Quantum number \"l9\""},
          Value{"n", "Quantum number \"n\" - Not in VAMDC"},
          Value{"parity", "Quantum number \"parity\""},
          Value{"r", "Quantum number \"r\""},
          Value{"rotSym", "Quantum number \"rotSym\""},
          Value{"rovibSym", "Quantum number \"rovibSym\""},
          Value{"sym", "Quantum number \"sym\""},
          Value{"tau", "Quantum number \"tau\" - Not in VAMDC"},
          Value{"term", "Quantum number \"term\" - Not in VAMDC"},
          Value{"v", "Quantum number \"v\""},
          Value{"v1", "Quantum number \"v1\""},
          Value{"v10", "Quantum number \"v10\""},
          Value{"v11", "Quantum number \"v11\""},
          Value{"v12", "Quantum number \"v12\""},
          Value{"v2", "Quantum number \"v2\""},
          Value{"v3", "Quantum number \"v3\""},
          Value{"v4", "Quantum number \"v4\""},
          Value{"v5", "Quantum number \"v5\""},
          Value{"v6", "Quantum number \"v6\""},
          Value{"v7", "Quantum number \"v7\""},
          Value{"v8", "Quantum number \"v8\""},
          Value{"v9", "Quantum number \"v9\""},
          Value{"vibInv", "Quantum number \"vibInv\""},
          Value{"vibRefl", "Quantum number \"vibRefl\""},
          Value{"vibSym", "Quantum number \"vibSym\""},
      }});

  opts.emplace_back(EnumeratedOption{
      .name            = "IsoRatioOption",
      .desc            = R"(The type of isotopologue ratio to use.
)",
      .values_and_desc = {
          Value{"Builtin", "Use the built-in isotopologue ratio"},
          Value{"Hitran", "Use the HITRAN isotopologue ratio"},
          Value{"None", "Do not use an isotopologue ratio"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeModelType",
      .desc = R"(The type of line shape model to use.
)",
      .values_and_desc =
          {
              Value{"T0", ":math:`X_0`"},
              Value{"T1", R"(:math:`X_0 \left(\frac{T_0}{T}\right) ^ {X_1}`)"},
              Value{
                  "T2",
                  R"(:math:`X_0 \left(\frac{T_0}{T}\right) ^ {X_1} \left[1 + X_2 \log\left(\frac{T_0}{T}\right)\right]`)"},
              Value{"T3", R"(:math:`X_0 + X_1 \left(T - T_0\right)`)"},
              Value{
                  "T4",
                  R"(:math:`\left[X_0 + X_1 \left(\frac{T_0}{T} - 1\right)\right] \left(\frac{T_0}{T}\right)^{X_2}`)"},
              Value{
                  "T5",
                  R"(:math:`X_0 \left(\frac{T_0}{T}\right)^{\frac{1}{4} + \frac{3}{2}X_1}`)"},
              Value{
                  "AER",
                  R"(:math:`X(200) = X_0`; :math:`X(250) = X_1`; :math:`X(298) = X_2`; :math:`X(340) = X_3`;  Linear interpolation in between)"},
              Value{
                  "DPL",
                  R"(:math:`X_0 \left(\frac{T_0}{T}\right) ^ {X_1} + X_2 \left(\frac{T_0}{T}\right) ^ {X_3}`)"},
              Value{"POLY",
                    R"(:math:`X_0 + X_1 T + X_2 T ^ 2 + X_3 T ^ 3 + \cdots`)"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeModelCoefficient",
      .desc =
          R"(The type of line shape model coefficients.  See :class:`~pyarts.arts.LineShapeModelType` for more information.
)",
      .values_and_desc =
          {
              Value{"X0", "x0", "X_0", "x_0", ":math:`X_0`"},
              Value{"X1", "x1", "X_1", "x_1", ":math:`X_1`"},
              Value{"X2", "x2", "X_2", "x_2", ":math:`X_2`"},
              Value{"X3", "x3", "X_3", "x_3", ":math:`X_3`"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeModelVariable",
      .desc = R"(The type of line shape model variable.
)",
      .values_and_desc =
          {
              Value{"G0", "Pressure broadening speed-independent"},
              Value{"D0", "Pressure f-shifting speed-dependent"},
              Value{"G2", "Pressure broadening speed-dependent"},
              Value{"D2", "Pressure f-shifting speed-independent"},
              Value{"FVC", "Frequency of velocity-changing collisions"},
              Value{"ETA", "Correlation"},
              Value{"Y", "First order line mixing coefficient"},
              Value{"G", "Second order line mixing coefficient"},
              Value{"DV", "Second order line mixing f-shifting"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "PathPositionType",
      .desc = R"(A type of position in a path.
)",
      .values_and_desc =
          {
              Value{"atm", "Atmospheric position"},
              Value{"space", "Space position"},
              Value{"surface", "Surface position"},
              Value{"unknown", "Unknown position"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineByLineVariable",
      .desc = R"(A type of line by line variable.
)",
      .values_and_desc =
          {
              Value{"f0", "Central frequency [Hz]"},
              Value{"e0", "Lower level energy [J]"},
              Value{"a", "Einstein coefficient [1/s]"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineByLineCutoffType",
      .desc = R"(A type of line by line cutoff.
)",
      .values_and_desc =
          {
              Value{"None", "No cutoff"},
              Value{
                  "ByLine",
                  "Line's are cut 1-by-1 around :attr:`~pyarts.arts.AbsorptionLine.f0`"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineByLineLineshape",
      .desc = R"(A type of line shape for line by line calculations.
)",
      .values_and_desc =
          {
              Value{"VP_LTE", "Voigt in local thermodynamic equilibrium"},
              Value{
                  "VP_LTE_MIRROR",
                  "Voigt in local thermodynamic equilibrium with negative frequency lines"},
              Value{
                  "VP_LINE_NLTE",
                  "Voigt in non-local thermodynamic equilibrium with level-by-level data"},
              Value{
                  "VP_ECS_MAKAROV",
                  "Voigt using Makarov's method of error-corrected sudden for line mixing of O2"},
              Value{
                  "VP_ECS_HARTMANN",
                  "Voigt using Hartmann's method of error-corrected sudden for line mixing of CO2"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "PlanetOrMoonType",
      .desc = R"(The type of planetary body that should be considered.
)",
      .values_and_desc =
          {
              Value{"Earth", "Planet is Earth"},
              Value{"Io", "Planet is Io"},
              Value{"Jupiter", "Planet is Jupiter"},
              Value{"Mars", "Planet is Mars"},
              Value{"Venus", "Planet is Venus"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "MissingFieldComponentError",
      .desc = R"(What kind of error is it to miss a field component?
)",
      .values_and_desc =
          {
              Value{"Throw", "It is required to have the field component."},
              Value{
                  "Zero",
                  "It is fine to zero-out the field component if it is missing."},
              Value{
                  "Ignore",
                  "It is fine to ignore that the field component is missing."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "HydrostaticPressureOption",
      .desc =
          R"(These options control how the hydrostatic pressure is computed in ARTS.

There are two main options for how the hydrostatic pressure is computed:

- ``HydrostaticEquation``:
  
  .. math::
    p(z + h) = p(z) \left(1 - \frac{Mg}{RT} h\right) 

- ``HypsometricEquation``:
  
  .. math::
    p(z + h) = p(z) e^{- \frac{Mg}{RT} h } 

where:

- :math:`p(z)` is the pressure at altitude :math:`z`,
- :math:`h` is the altitude difference,
- :math:`g` is the gravitational acceleration at altitude :math:`z`,
- :math:`M` is the mean mass of an atmospheric molecule at altitude :math:`z`,
- :math:`R` is the specific gas constant, and
- :math:`T` is the temperature at altitude :math:`z`.

.. note::
  We do not consider The Eötvös effect.
  We also do not directly support "virtual temperature", as water is baked into :math:`M`.
  Either addition would be welcome as additional options.
)",
      .values_and_desc =
          {
              Value{"HydrostaticEquation",
                    "Piece-wise linearly decaying with altitude."},
              Value{"HypsometricEquation",
                    "Piece-wise exponentailly decaying with altitude."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "AbsorptionBandSortingOption",
      .desc = R"(What kind of sorting should the bands be sorted by?
)",
      .values_and_desc =
          {
              Value{"IntegratedIntensity",
                    "By the total absorption strenght of the band."},
              Value{"FrontFrequency",
                    "By the lowest frequency of the first line."},
              Value{"None", "Do not change the sorting."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "FieldComponent",
      .desc = R"(Selection of a field component
)",
      .values_and_desc =
          {
              Value{"u", "U", "East component"},
              Value{"v", "V", "North component"},
              Value{"w", "W", "Up component"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "SpectralRadianceUnitType",
      .desc = R"(Choice of spectral radiance unit in conversions.

The conversion to one of these units is expected to happen from the internal representation of *spectral_radiance*,
which is [W / m :math:`^{2}` Hz sr].

For the description below, :math:`c` is the speed of light, :math:`k` is the Boltzmann constant,
:math:`f` is the frequency, :math:`h` is the Planck constant, :math:`F(x)` is the conversion
function, and :math:`[I,\; Q,\; U,\; V]` is the internal representation of *spectral_radiance*.

For Rayleigh-Jeans brightness temperature the conversion is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [F(I),\; F(Q),\; F(U),\; F(V)]

.. math::

    F(x) = \frac{c^2 x}{2 k f^2}.

For Planck brightness temperature the conversion is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow
        \left[
        F(I),\;
        F\left( \frac{I + Q}{2} \right) - F\left( \frac{I - Q}{2} \right),\;
        F\left( \frac{I + U}{2} \right) - F\left( \frac{I - U}{2} \right),\;
        F\left( \frac{I + V}{2} \right) - F\left( \frac{I - V}{2} \right)
        \right]

.. math::

    F(x) = \frac{h f}{k}\left[\log\left(1 + \frac{2hf^3}{c^2x}\right)\right]^{-1}.

The spectral radiance per wavelength [W / m :math:`^{2}` m sr] is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [F(I),\; F(Q),\; F(U),\; F(V)]

.. math::

    F(x) = \frac{f^2 x}{c}.

The spectral radiance per wavenumber [W / m :math:`^{2}` m :math:`^{-1}` sr] is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [F(I),\; F(Q),\; F(U),\; F(V)]

.. math::

    F(x) = \frac{x}{c}.

Lastly, the unit option of course just retains the current state [W / m :math:`^{2}` Hz sr].

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [I,\; Q,\; U,\; V].
)",
      .values_and_desc =
          {
              Value{"RJBT",
                    "Tr",
                    R"(Rayleigh-Jeans brightness temperature [K]
)"},
              Value{"PlanckBT",
                    "Tb",
                    R"(Planck brightness temperature [K]
)"},
              Value{"W_m2_m_sr",
                    "W/(m^2 m sr)",
                    "Spectral radiance wavelength [W / m :math:`^{2}` m sr]"},
              Value{
                  "W_m2_m1_sr",
                  "W/(m^2 m-1 sr)",
                  "Spectral radiance wavenumber [W / m :math:`^{2}` m :math:`^{-1}` sr]"},
              Value{"unit",
                    "1",
                    "Unit spectral radiance [W / m :math:`^{2}` Hz sr]"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "FileType",
      .desc = "A choice of file format types.\n",
      .values_and_desc =
          {
              Value{"ascii", "ASCII", "Ascii", "text", "Save as ASCII"},
              Value{"zascii", "ZASCII", "Zip", "zip", "Save as zipped ASCII"},
              Value{"binary", "BINARY", "Binary", "bin", "Save as binary data"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "EarthEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {Value{"Sphere", R"(A spherical Earth.)"},
           Value{"WGS84",
                 R"(The reference ellipsoid used by the GPS system.)"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "IoEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{
                  "Sphere",
                  "A spherical planetesimal.  The radius is taken from a report of the IAU/IAG Working Group."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "EuropaEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{"Sphere", "A spherical planetesimal."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "GanymedeEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{"Sphere", "A spherical planetesimal."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "JupiterEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{
                  "Sphere",
                  "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."},
              Value{
                  "Ellipsoid",
                  "A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "MarsEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{
                  "Sphere",
                  "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."},
              Value{
                  "Ellipsoid",
                  "A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "MoonEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{
                  "Sphere",
                  "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."},
              Value{
                  "Ellipsoid",
                  "A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "VenusEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{
                  "Sphere",
                  "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "SpeciesTagType",
      .desc = "Type of species tag.\n",
      .values_and_desc =
          {
              Value{"Plain", "A plain species, one isotopologue or all."},
              Value{"Predefined", "A predefined model species."},
              Value{"Cia", "A pair of collision-induced species."},
              Value{"XsecFit", "A cross-section fitting of a model species."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "PolarizationChoice",
      .desc = R"(Named polarization states to help create relevant *Stokvec*.

Note that these are just user fri'\n'y suggestions and it is recommended to
create the correct *Stokvec* manually if the desired polarization state is not
represented below.

Also, be aware that the unit of, e.g., *spectral_radiance* (often the last choice
of *SpectralRadianceUnitType*) is in Kelvin, then the code below will give the
measured brightness temperature in Kelvin for these polarization states, but that
if the unit is still in Watts of any kind, then the code below will give 2 times
the select polarized brightness temperatures (but still the correct unpolarized
radiation).
)",
      .values_and_desc =
          {Value{"I", "I", "No polarization state [1, 0, 0, 0]."},
           Value{
               "Q",
               "Q",
               "Difference between vertical and horizontal linear polarization [0, 1, 0, 0]."},
           Value{
               "U",
               "U",
               "Difference between plus and minus 45 degrees linear polarization [0, 0, 1, 0]."},
           Value{
               "V",
               "V",
               "Difference between right and left circular polarization [0, 0, 0, 1]."},
           Value{"Iv", "Iv", "Vertical linear polarization [1, 1, 0, 0]."},
           Value{"Ih", "Ih", "Horizontal linear polarization [1, -1, 0, 0]."},
           Value{"Ip45",
                 "I+45",
                 "Plus 45 degrees linear polarization [1, 0, 1, 0]."},
           Value{"Im45",
                 "I-45",
                 "Minus 45 degrees linear polarization [1, 0, -1, 0]."},
           Value{"Ilhc", "LC", "Left circular polarization [1, 0, 0, -1]."},
           Value{"Irhc", "RC", "Right circular polarization [1, 0, 0, 1]."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "ParticulateProperty",
      .desc =
          R"(Numerical properties used to numerically represent particle populations.
)",
      .values_and_desc =
          {Value{"MassDensity", "m", "Mass density in kg/m^{-3}"},
           Value{"NumberDensity", "n", "Number density in m^{-3}"},
           Value{"DMax", "dmax", "Maximum particle diameter in m."},
           Value{"DVeq", "dveq", "Volume-equivalent diameter in m."},
           Value{"Extinction", "ext", "Extinction in m^{-1}"},
           Value{"SingleScatteringAlbedo", "ssa", "Single scattering albedo"},
           Value{"ShapeParameter",
                 "ShapeParameter",
                 "PSD shape parmeter in arbitary units."},
           Value{"IntercepParameter",
                 "InterceptParameter",
                 "PSD intercept parameter in arbitary units."}},
  });
  opts.emplace_back(EnumeratedOption{
      .name = "SizeParameter",
      .desc =
          R"(Parameters used to represent the size of particles.
)",
      .values_and_desc =
          {Value{"Mass", "m", "Mass in kg"},
           Value{"DMax", "d_max", "Maximum diameter in m"},
           Value{"DVeq", "d_veq", "Volume-equivalent diameter in m"}},
  });

  fix_static(opts);
  return opts;
}
}  // namespace

const std::vector<EnumeratedOption>& internal_options() {
  static std::vector<EnumeratedOption> opts = internal_options_create();
  return opts;
}

std::string_view EnumeratedOption::sz() const {
  const auto n = values_and_desc.size();

  if (n < std::numeric_limits<char>::max()) return "char";
  if (n < std::numeric_limits<short>::max()) return "short";
  if (n < std::numeric_limits<int>::max()) return "int";
  if (n < std::numeric_limits<long>::max()) return "long";
  if (n < std::numeric_limits<long long>::max()) return "long long";

  throw std::runtime_error("Too many values for enum class");
}

std::string add_spaces(const std::string& s, int numspaces) {
  std::string out{s};
  const std::string newline{"\n"};
  const std::string spaces = std::string(numspaces, ' ') + '\n';
  replace(out, newline, spaces);
  return out;
}

std::string EnumeratedOption::docs() const {
  std::string out;

  const auto n = values_and_desc.front().size();

  std::format_to(std::back_inserter(out), "{}\n\nValid options:\n\n", desc);
  for (auto& v : values_and_desc) {
    std::string_view x = "- "sv;
    for (auto& s : v | stdv::take(n - 1)) {
      std::format_to(std::back_inserter(out),
                     R"({}``"{}"``)",
                     std::exchange(x, " or "sv),
                     s);
    }
    std::format_to(
        std::back_inserter(out), ":\n  {}\n", add_spaces(v.back(), 2));
  }

  return out;
}

std::string EnumeratedOption::tail() const {
  std::string out;

  const auto m = values_and_desc.size();
  const auto n = values_and_desc.front().size();

  // Create good enum function
  std::format_to(std::back_inserter(out),
                 R"(
template<> constexpr bool good_enum<{0}>({0} x) noexcept {{
    const auto v = static_cast<std::size_t>(x);
    return v < {1};
}}

template<> struct enumdocs<{0}> {{
  static std::string_view str() noexcept;
  static constexpr std::string_view name = "{0}"sv;
}};

namespace enumtyps {{
    inline constexpr std::array {0}Types = {{)",
                 name,
                 m);

  for (const auto& v : values_and_desc)
    std::format_to(std::back_inserter(out), "{}::{}, ", name, v[0]);

  std::format_to(std::back_inserter(out), R"(}};
}}  // namespace enumtyps

namespace enumstrs {{
)");

  // Create enumstrs arrays
  for (std::size_t i = 0; i < n - 1; i++) {
    std::format_to(std::back_inserter(out),
                   R"(
  template <> struct enum_str_data<{}, {}> {{
    static constexpr std::array strs={{)",
                   name,
                   i);

    for (const auto& v : values_and_desc)
      std::format_to(std::back_inserter(out), R"("{}"sv, )", v[i]);

    std::format_to(std::back_inserter(out), "}};\n  }};\n");
  }

  std::format_to(std::back_inserter(out),
                 R"(
  template <int i={0}>
  inline constexpr auto {1}Names = enum_str_data<{1}, i>::strs;
}}  // namespace enumstrs

template <int i={0}> constexpr std::string_view toString({1} x)  requires(i >= 0 and i < {2}) {{
  if (good_enum(x))
    return enumstrs::{1}Names<i>[static_cast<std::size_t>(x)];
  return "BAD {1}"sv;
}}

template<> constexpr {1} to<{1}>(const std::string_view x) {{
  using namespace enumstrs;
  using namespace enumtyps;)",
                 preferred_print,
                 name,
                 n - 1);

  for (std::size_t i = 0; i < n - 1; i++) {
    std::format_to(std::back_inserter(out),
                   R"(
  if (const auto i = std::ranges::distance(std::ranges::begin({0}Names<{1}>),
                                           std::ranges::find({0}Names<{1}>, x));
                     i < {2})
    return {0}Types[i];)",
                   name,
                   i,
                   m);
  }
  std::format_to(std::back_inserter(out),
                 R"(
  throw std::runtime_error(R"-x-(Bad value.

See https://atmtools.github.io/arts-docs-master/pyarts.arts.{0}.html for valid options.
)-x-");
}}

namespace enumsize {{ inline constexpr std::size_t {0}Size = {2}; }}

std::ostream& operator<<(std::ostream &os, const {0} x);

std::istream& operator>>(std::istream &is, {0}& x);

template<> struct std::formatter<{0}> {{
  using T={0};
    
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() {{ return *this; }}
  [[nodiscard]] constexpr auto& inner_fmt() const {{ return *this; }}

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {{
    return parse_format_tags(tags, ctx);
  }}

  template <class FmtContext>
  FmtContext::iterator format(const T& v, FmtContext& ctx) const {{
    return std::format_to(ctx.out(), "{{0}}{{1}}{{0}}", tags.quote(), toString<{3}>(v));
  }}
}};)",
                 name,
                 docs(),
                 m,
                 preferred_print);

  return out;
}

std::string EnumeratedOption::head() const {
  static constexpr std::array reserved = {"get_options",
                                          "get_options_as_string"};

  for (auto& v : values_and_desc) {
    for (auto& s : v) {
      if (stdr::any_of(reserved, [&s](auto& r) { return r == s; })) {
        throw std::runtime_error("Reserved word used for enum class " + name);
      }
    }
  }

  [*this]() {
    std::vector<std::string> all;
    all.reserve(values_and_desc.size() * (values_and_desc.front().size()));
    for (auto v : values_and_desc) {
      stdr::sort(v);
      auto last = std::unique(v.begin(), v.end());
      all.insert(all.end(), v.begin(), last);
    }
    stdr::sort(all);
    auto last = std::unique(all.begin(), all.end());
    if (last != all.end()) {
      throw std::runtime_error("Duplicate value for enum class \"" + name +
                               "\" \"" + *last + "\"");
    }
  }();

  if (values_and_desc.empty()) {
    throw std::runtime_error("No values for enum class " + name);
  }
  const auto n = values_and_desc.front().size();

  if (n < 2) {
    throw std::runtime_error(
        "No descriptions for enum class " + name +
        std::string{
            ".\nThe Value should be {main_key, alt_keys..., description}, where alt_keys... can be omitted"});
  }

  std::string out;

  std::format_to(
      std::back_inserter(out), "enum class {} : {} {{\n", name, sz());

  for (const auto& v : values_and_desc) {
    if (v.empty()) {
      throw std::runtime_error("Empty value for enum class " + name);
    }

    if (v.size() != n) {
      throw std::runtime_error("Mismatched value count for enum class " + name +
                               " " + v.front());
    }

    std::format_to(std::back_inserter(out), "{}, ", v[0]);
  }

  std::format_to(std::back_inserter(out), "}};\n\n");

  return out;
}

std::string EnumeratedOption::impl() const {
  std::string out;

  std::format_to(std::back_inserter(out),
                 R"(
std::ostream &operator<<(std::ostream &os, const {0} x) {{
  std::print(os, "{{}}", toString<{1}>(x));
  return os;
}}

std::istream &operator>>(std::istream &is, {0}& x) {{
  std::string s;
  is >> s;
  try {{
    x = to<{0}>(s);
  }} catch (const std::exception &e) {{
    throw std::runtime_error(
        std::format("Failed to read {0} from input stream value \"{{}}\":\n{{}}", x, e.what()));
  }}
  return is;
}}

std::string_view enumdocs<{0}>::str() noexcept {{
  return R"-ENUMDOC-({2})-ENUMDOC-"sv;
}}

static_assert(arts_formattable<{0}>, "Not good: {0}");
)",
                 name,
                 preferred_print,
                 docs());

  return out;
}
