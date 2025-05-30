#include "arts_options.h"

#include <algorithm>
#include <array>
#include <limits>
#include <ranges>
#include <sstream>
#include <utility>
#include <vector>

using Value = std::vector<std::string>;

namespace {
std::vector<EnumeratedOption> internal_options_create() {
  std::vector<EnumeratedOption> opts;

  opts.emplace_back(EnumeratedOption{
      .name = "ray_path_observer_agendaSetGeometricMaxStep",
      .desc =
          R"(For use with *ray_path_observer_agendaSetGeometric*.  Determines how to densify the geometric path.
)",
      .values_and_desc = {
          Value{"half", "1/2", "Use *ray_pathFillGeometricHalfStep*"},
          Value{"step", "linear", "Use *ray_pathFillGeometricStepwise*"},
          Value{"None", "0", "Do not use any step filling method"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "disort_settings_agenda_setup_layer_emission_type",
      .desc =
          R"(For atmospheric emission settings with *disort_settings_agendaSetup*.
)",
      .values_and_desc = {
          Value{"None", "No atmospheric emission"},
          Value{"LinearInTau", "Layer emission is linear in optical depth"},
          Value{"LinearInTauNonLTE", "Layer emission is linear in optical depth taking non-LTE into account"},
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
          Value{"f", "Frequency", "Frequency of the sensor"},
          Value{"za", "PointingZenith", "Pointing Zenith of the sensor"},
          Value{"aa", "PointingAzimuth", "Pointing Azimuth of the sensor"},
          Value{"alt", "PointingAltitude", "Pointing Altitude of the sensor"},
          Value{"lat", "PointingLatitude", "Pointing Latitude of the sensor"},
          Value{"lon", "PointingLongitude", "Pointing Longitude of the sensor"},
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
      .values_and_desc = {
          Value{"Bath", "AIR", "Any non-descript species"},
          Value{"Water", "H2O", "Water moluecule"},
          Value{"CarbonDioxide", "CO2", "Carbon Dioxide moluecule"},
          Value{"Ozone", "O3", "Ozone moluecule"},
          Value{"NitrogenOxide", "N2O", "Nitrogen Oxide moluecule"},
          Value{"CarbonMonoxide", "CO", "Carbon Monoxide moluecule"},
          Value{"Methane", "CH4", "Methane moluecule"},
          Value{"Oxygen", "O2", "Oxygen moluecule"},
          Value{"NitricOxide", "NO", "Nitric Oxide moluecule"},
          Value{"SulfurDioxide", "SO2", "Sulfur Dioxide moluecule"},
          Value{"NitrogenDioxide", "NO2", "Nitrogen Dioxide moluecule"},
          Value{"Ammonia", "NH3", "Ammonia moluecule"},
          Value{"NitricAcid", "HNO3", "Nitric Acid moluecule"},
          Value{"Hydroxyl", "OH", "Hydroxyl moluecule"},
          Value{"HydrogenFluoride", "HF", "Hydrogen Fluoride moluecule"},
          Value{"HydrogenChloride", "HCl", "Hydrogen Chloride moluecule"},
          Value{"HydrogenBromide", "HBr", "Hydrogen Bromide moluecule"},
          Value{"HydrogenIodide", "HI", "Hydrogen Iodide moluecule"},
          Value{"ChlorineMonoxide", "ClO", "Chlorine Monoxide moluecule"},
          Value{"CarbonylSulfide", "OCS", "Carbonyl Sulfide moluecule"},
          Value{"Formaldehyde", "H2CO", "Formaldehyde moluecule"},
          Value{"HeavyFormaldehyde", "HDCO", "Heavy Formaldehyde moluecule"},
          Value{"VeryHeavyFormaldehyde",
                "D2CO",
                "Very Heavy Formaldehyde moluecule"},
          Value{"HypochlorousAcid", "HOCl", "Hypochlorous Acid moluecule"},
          Value{"Nitrogen", "N2", "Nitrogen moluecule"},
          Value{"HydrogenCyanide", "HCN", "Hydrogen Cyanide moluecule"},
          Value{"Chloromethane", "CH3Cl", "Chloromethane moluecule"},
          Value{"HydrogenPeroxide", "H2O2", "Hydrogen Peroxide moluecule"},
          Value{"Acetylene", "C2H2", "Acetylene moluecule"},
          Value{"Ethane", "C2H6", "Ethane moluecule"},
          Value{"Phosphine", "PH3", "Phosphine moluecule"},
          Value{"CarbonylFluoride", "COF2", "Carbonyl Fluoride moluecule"},
          Value{"SulfurHexafluoride", "SF6", "Sulfur Hexafluoride moluecule"},
          Value{"HydrogenSulfide", "H2S", "Hydrogen Sulfide moluecule"},
          Value{"FormicAcid", "HCOOH", "Formic Acid moluecule"},
          Value{"LeftHeavyFormicAcid",
                "DCOOH",
                "Left Heavy Formic Acid moluecule"},
          Value{"RightHeavyFormicAcid",
                "HCOOD",
                "Right Heavy Formic Acid moluecule"},
          Value{"Hydroperoxyl", "HO2", "Hydroperoxyl moluecule"},
          Value{"OxygenAtom", "O", "Oxygen atom"},
          Value{"ChlorineNitrate", "ClONO2", "Chlorine Nitrate moluecule"},
          Value{"NitricOxideCation", "NO+", "Nitric Oxide Cation moluecule"},
          Value{"HypobromousAcid", "HOBr", "Hypobromous Acid moluecule"},
          Value{"Ethylene", "C2H4", "Ethylene moluecule"},
          Value{"Methanol", "CH3OH", "Methanol moluecule"},
          Value{"Bromomethane", "CH3Br", "Bromomethane moluecule"},
          Value{"Acetonitrile", "CH3CN", "Acetonitrile moluecule"},
          Value{"HeavyAcetonitrile", "CH2DCN", "Heavy Acetonitrile moluecule"},
          Value{"CarbonTetrafluoride", "CF4", "Carbon Tetrafluoride moluecule"},
          Value{"Diacetylene", "C4H2", "Diacetylene moluecule"},
          Value{"Cyanoacetylene", "HC3N", "Cyanoacetylene moluecule"},
          Value{"Hydrogen", "H2", "Hydrogen moluecule"},
          Value{"CarbonMonosulfide", "CS", "Carbon Monosulfide moluecule"},
          Value{"SulfurTrioxide", "SO3", "Sulfur Trioxide moluecule"},
          Value{"Cyanogen", "C2N2", "Cyanogen moluecule"},
          Value{"Phosgene", "COCl2", "Phosgene moluecule"},
          Value{"SulfurMonoxide", "SO", "Sulfur Monoxide moluecule"},
          Value{"CarbonDisulfide", "CS2", "Carbon Disulfide moluecule"},
          Value{"Methyl", "CH3", "Methyl moluecule"},
          Value{"Cyclopropene", "C3H4", "Cyclopropene moluecule"},
          Value{"SulfuricAcid", "H2SO4", "Sulfuric Acid moluecule"},
          Value{"HydrogenIsocyanide", "HNC", "Hydrogen Isocyanide moluecule"},
          Value{"BromineMonoxide", "BrO", "Bromine Monoxide moluecule"},
          Value{"ChlorineDioxide", "OClO", "Chlorine Dioxide moluecule"},
          Value{"Propane", "C3H8", "Propane moluecule"},
          Value{"Helium", "He", "Helium atom"},
          Value{"ChlorineMonoxideDimer",
                "Cl2O2",
                "Chlorine Monoxide Dimer moluecule"},
          Value{"HydrogenAtom", "H", "Hydrogen atom"},
          Value{"Argon", "Ar", "Argon atom"},
          Value{"Hexafluoroethane", "C2F6", "Hexafluoroethane moluecule"},
          Value{"Perfluoropropane", "C3F8", "Perfluoropropane moluecule"},
          Value{"Perfluorobutane", "C4F10", "Perfluorobutane moluecule"},
          Value{"Perfluoropentane", "C5F12", "Perfluoropentane moluecule"},
          Value{"Perfluorohexane", "C6F14", "Perfluorohexane moluecule"},
          Value{"Perfluorooctane", "C8F18", "Perfluorooctane moluecule"},
          Value{"Perfluorocyclobutane",
                "cC4F8",
                "Perfluorocyclobutane moluecule"},
          Value{
              "CarbonTetrachloride", "CCl4", "Carbon Tetrachloride moluecule"},
          Value{"CFC11", "CFC11", "CFC11 moluecule"},
          Value{"CFC113", "CFC113", "CFC113 moluecule"},
          Value{"CFC114", "CFC114", "CFC114 moluecule"},
          Value{"CFC115", "CFC115", "CFC115 moluecule"},
          Value{"CFC12", "CFC12", "CFC12 moluecule"},
          Value{"Dichloromethane", "CH2Cl2", "Dichloromethane moluecule"},
          Value{"Trichloroethane", "CH3CCl3", "Trichloroethane moluecule"},
          Value{"Trichloromethane", "CHCl3", "Trichloromethane moluecule"},
          Value{"Bromochlorodifluoromethane",
                "Halon1211",
                "Bromochlorodifluoromethane moluecule"},
          Value{"Bromotrifluoromethane",
                "Halon1301",
                "Bromotrifluoromethane moluecule"},
          Value{"Dibromotetrafluoroethane",
                "Halon2402",
                "Dibromotetrafluoroethane moluecule"},
          Value{"HCFC141b", "HCFC141b", "HCFC141b moluecule"},
          Value{"HCFC142b", "HCFC142b", "HCFC142b moluecule"},
          Value{"HCFC22", "HCFC22", "HCFC22 moluecule"},
          Value{"HFC125", "HFC125", "HFC125 moluecule"},
          Value{"HFC134a", "HFC134a", "HFC134a moluecule"},
          Value{"HFC143a", "HFC143a", "HFC143a moluecule"},
          Value{"HFC152a", "HFC152a", "HFC152a moluecule"},
          Value{"HFC227ea", "HFC227ea", "HFC227ea moluecule"},
          Value{"HFC23", "HFC23", "HFC23 moluecule"},
          Value{"HFC236fa", "HFC236fa", "HFC236fa moluecule"},
          Value{"HFC245fa", "HFC245fa", "HFC245fa moluecule"},
          Value{"HFC32", "HFC32", "HFC32 moluecule"},
          Value{"HFC365mfc", "HFC365mfc", "HFC365mfc moluecule"},
          Value{"NitrogenTrifluoride", "NF3", "Nitrogen Trifluoride moluecule"},
          Value{"SulfurylFluoride", "SO2F2", "Sulfuryl Fluoride moluecule"},
          Value{"HFC4310mee", "HFC4310mee", "HFC4310mee moluecule"},
          Value{"Germane", "GeH4", "Germane moluecule"},
          Value{"Iodomethane", "CH3I", "Iodomethane moluecule"},
          Value{"Fluoromethane", "CH3F", "Fluoromethane moluecule"},
          Value{"liquidcloud", "liquidcloud", "liquidcloud tag"},
          Value{"icecloud", "icecloud", "icecloud tag"},
          Value{"rain", "rain", "rain tag"},
          Value{"free_electrons", "free_electrons", "free electrons tag"},
          Value{"particles", "particles", "particles tag"},
      }});

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
              Value{"subsurface", "Sub-surface position"},
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
      .desc = R"(What kind of error is it to miss a field component?
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
      .name            = "EarthEllipsoid",
      .desc            = "Choice of ellipsoid.\n",
      .values_and_desc = {Value{"Sphere", R"(A spherical Earth.

  The radius is set following the value set for the Earth radius.)"},
                          Value{
                              "WGS84",
                              R"(The reference ellipsoid used by the GPS system.

  Should be the standard choice for a non-spherical Earth.)"}},
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

  return opts;
}
}  // namespace

const std::vector<EnumeratedOption>& internal_options() {
  static std::vector<EnumeratedOption> opts = internal_options_create();
  return opts;
}

std::string_view EnumeratedOption::sz() const {
  const auto n = values_and_desc.size();

  if (n < std::numeric_limits<bool>::max()) return "bool";
  if (n < std::numeric_limits<char>::max()) return "char";
  if (n < std::numeric_limits<short>::max()) return "short";
  if (n < std::numeric_limits<int>::max()) return "int";
  if (n < std::numeric_limits<long>::max()) return "long";
  if (n < std::numeric_limits<long long>::max()) return "long long";

  throw std::runtime_error("Too many values for enum class");
}

std::string EnumeratedOption::docs() const {
  std::ostringstream os;

  const auto n = values_and_desc.front().size();

  os << desc << "\n\nValid options:\n\n";
  for (auto& v : values_and_desc) {
    std::string_view x = "- ";
    for (auto& s : v | std::views::take(n - 1)) {
      os << std::exchange(x, " or ") << "``" << '"' << s << '"' << "``";
    }
    os << ": " << v.back() << '\n';
  }

  return os.str();
}

namespace {
int format_ind(const std::string& name) {
  if (name == "SpeciesEnum") return 1;
  return 0;
}
}  // namespace

std::string EnumeratedOption::tail() const {
  std::ostringstream os;

  const auto m = values_and_desc.size();
  const auto n = values_and_desc.front().size();

  // Create good enum function
  os << "template<> constexpr bool good_enum<" << name << ">(" << name
     << " x) noexcept {\n  const auto v = static_cast<std::size_t>(x);\n  return v < "
     << m << ";\n}\n\n";

  // Create enumtyps array
  os << "namespace enumtyps {\n  inline constexpr std::array " << name
     << "Types = {";
  for (const auto& v : values_and_desc) {
    os << name << "::" << v[0] << ", ";
  }
  os << "};\n}  // namespace enumtyps\n\n";

  // Create enumstrs arrays
  os << "namespace enumstrs {\n";
  for (std::size_t i = 0; i < n - 1; i++) {
    os << "  template <> struct enum_str_data<" << name << ", " << i << "> {\n";
    os << "    static constexpr std::array strs={";
    for (const auto& v : values_and_desc) {
      os << '"' << v[i] << "\"sv, ";
    }
    os << "};\n  };\n";
  }
  os << "  template <int i=" << format_ind(name)
     << ">\n  inline constexpr auto " << name << "Names = enum_str_data<"
     << name << ", i>::strs;\n";
  os << "}  // namespace enumstrs\n\n";

  // Create toString functions
  os << "template <int i=" << format_ind(name)
     << "> constexpr const std::string_view toString(" << name
     << " x) requires(i >= 0 and i < " << n - 1
     << ") {\n"
        "  if (good_enum(x))\n    return enumstrs::"
     << name << "Names<i>[static_cast<std::size_t>(x)];\n  return \"BAD "
     << name << "\";\n}\n\n";

  // Create toEnumType functions
  os << "template<> constexpr " << name << " to<" << name
     << ">(const std::string_view x) {\n  using namespace enumstrs;\n";
  for (std::size_t i = 0; i < n - 1; i++) {
    os << "  if (const auto i = std::distance(" << name << "Names<" << i
       << ">.begin(), std::ranges::find(" << name << "Names<" << i
       << ">, x)); i < " << m
       << ")\n"
          "    return enumtyps::"
       << name << "Types[i];\n";
  }
  os << "  throw std::runtime_error(std::string{\"Bad value: \\\"\"} + std::string{x} + R\"-x-(\"\n\n"
     << docs() << ")-x-\");\n}\n\n";

  os << "namespace enumsize { inline constexpr std::size_t " << name
     << "Size = " << m << "; }\n\n";

  // Create ostream operator
  os << "std::ostream &operator<<(std::ostream &os, const " << name
     << " x);\n\n";
  os << "std::istream &operator>>(std::istream &is, " << name << "& x);\n\n";

  // Create xml-io operator
  os << "void xml_read_from_stream(std::istream& is_xml, " << name
     << "& s, bifstream*);\n\n";
  os << "void xml_write_to_stream(std::ostream& os_xml, const " << name
     << "& s, bofstream*, const std::string&);\n\n";

  os << "template<> struct std::formatter<" << name << "> {\n  using T=" << name
     << ";\n";
  os << R"--(
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const T& v, FmtContext& ctx) const {
  constexpr Index IND = )--"
     << format_ind(name) << ";\n";

  os << R"--(
    const auto q = tags.quote();
    return std::format_to(ctx.out(), "{}{}{}", q, toString<IND>(v), q);
  }
)--";
  os << "};\n";

  return os.str();
}

std::string EnumeratedOption::head() const {
  static constexpr std::array reserved = {"get_options",
                                          "get_options_as_string"};

  for (auto& v : values_and_desc) {
    for (auto& s : v) {
      if (std::ranges::any_of(reserved, [&s](auto& r) { return r == s; })) {
        throw std::runtime_error("Reserved word used for enum class " + name);
      }
    }
  }

  [*this]() {
    std::vector<std::string> all;
    all.reserve(values_and_desc.size() * (values_and_desc.front().size()));
    for (auto v : values_and_desc) {
      std::ranges::sort(v);
      auto last = std::unique(v.begin(), v.end());
      all.insert(all.end(), v.begin(), last);
    }
    std::ranges::sort(all);
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

  std::ostringstream os;

  os << "enum class " << name << " : " << sz() << " {\n";
  for (const auto& v : values_and_desc) {
    if (v.empty()) {
      throw std::runtime_error("Empty value for enum class " + name);
    }

    if (v.size() != n) {
      throw std::runtime_error("Mismatched value count for enum class " + name +
                               " " + v.front());
    }

    os << v[0] << ", ";
  }
  os << "};\n\n";

  return os.str();
}

std::string EnumeratedOption::impl() const {
  std::ostringstream os;
  // Create ostream operator
  os << "std::ostream &operator<<(std::ostream &os, const " << name
     << " x) {\n  return os << toString<" << format_ind(name) << ">(x);\n}\n\n";
  os << "std::istream &operator>>(std::istream &is, " << name
     << "& x) {\n  std::string s;\n  is >> s;\n  x = to<" << name
     << ">(s);\n  return is;\n}\n\n";

  os << "void xml_read_from_stream(std::istream& is, " << name
     << "& s, bifstream*) {\n"
        "  std::string x;\n"
        "  is >> x;\n"
        "  if (x != \"<"
     << name << R"(>") throw std::runtime_error("Expected \"<)" << name
     << ">\\\" got: \\\"\" + x + \"\\\"\");\n"
        "  is >> x;\n"
        "  s = to<"
     << name
     << ">(x);\n"
        "  is >> x;\n"
        "  if (x != \"</"
     << name << R"(>") throw std::runtime_error("Expected \"</)" << name
     << ">\\\" got: \\\"\" + x + \"\\\"\");\n"
        "}\n\n";
  os << "void xml_write_to_stream(std::ostream& os, const " << name
     << "& s, bofstream*, const std::string&) {\n"
        "  os << \"<"
     << name << "> \" << toString<" << format_ind(name) << ">(s) << \" </"
     << name << ">\\n\";\n}\n\n";
  os << "static_assert(arts_formattable<" << name << ">, \"Not good: " << name
     << "\");\n\n";
  return os.str();
}
