#include "arts_options.h"

#include <iomanip>
#include <limits>
#include <ranges>
#include <sstream>
#include <utility>
#include <vector>

using Value = std::vector<std::string>;

std::vector<EnumeratedOption> internal_options_create() {
  std::vector<EnumeratedOption> opts;

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
      .name = "GridType",
      .desc = R"(Type of Lagrange interpolation weights.
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
                           Value{"h", "Altitude [m]"},
                           Value{"t", "Temperature [K]"},
                           Value{"wind_u", "Wind field U-component [m/s]"},
                           Value{"wind_v", "Wind field V-component [m/s]"},
                           Value{"wind_w", "Wind field W-component [m/s]"},
                       }});

  opts.emplace_back(
      EnumeratedOption{.name = "AtmKey",
                       .desc = R"(A key to identify an atmospheric property.
)",
                       .values_and_desc = {
                           Value{"t", "Temperature [K]"},
                           Value{"p", "Pressure [Pa]"},
                           Value{"wind_u", "Wind field U-component [m/s]"},
                           Value{"wind_v", "Wind field V-component [m/s]"},
                           Value{"wind_w", "Wind field W-component [m/s]"},
                           Value{"mag_u", "Magnetic field U-component [T]"},
                           Value{"mag_v", "Magnetic field V-component [T]"},
                           Value{"mag_w", "Magnetic field W-component [T]"},
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
          Value{"CFC11", "CFC11 ", "CFC11 moluecule"},
          Value{"CFC113", "CFC113 ", "CFC113 moluecule"},
          Value{"CFC114", "CFC114 ", "CFC114 moluecule"},
          Value{"CFC115", "CFC115 ", "CFC115 moluecule"},
          Value{"CFC12", "CFC12 ", "CFC12 moluecule"},
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
          Value{"HCFC141b", "HCFC141b ", "HCFC141b moluecule"},
          Value{"HCFC142b", "HCFC142b ", "HCFC142b moluecule"},
          Value{"HCFC22", "HCFC22 ", "HCFC22 moluecule"},
          Value{"HFC125", "HFC125 ", "HFC125 moluecule"},
          Value{"HFC134a", "HFC134a ", "HFC134a moluecule"},
          Value{"HFC143a", "HFC143a ", "HFC143a moluecule"},
          Value{"HFC152a", "HFC152a ", "HFC152a moluecule"},
          Value{"HFC227ea", "HFC227ea ", "HFC227ea moluecule"},
          Value{"HFC23", "HFC23 ", "HFC23 moluecule"},
          Value{"HFC236fa", "HFC236fa ", "HFC236fa moluecule"},
          Value{"HFC245fa", "HFC245fa ", "HFC245fa moluecule"},
          Value{"HFC32", "HFC32 ", "HFC32 moluecule"},
          Value{"HFC365mfc", "HFC365mfc ", "HFC365mfc moluecule"},
          Value{"NitrogenTrifluoride", "NF3", "Nitrogen Trifluoride moluecule"},
          Value{"SulfurylFluoride", "SO2F2", "Sulfuryl Fluoride moluecule"},
          Value{"HFC4310mee", "HFC4310mee ", "HFC4310mee moluecule"},
          Value{"Germane", "GeH4", "Germane moluecule"},
          Value{"Iodomethane", "CH3I", "Iodomethane moluecule"},
          Value{"Fluoromethane", "CH3F", "Fluoromethane moluecule"},
          Value{"liquidcloud", "LiquidCloud", "liquidcloud tag"},
          Value{"icecloud", "IceCloud", "icecloud tag"},
          Value{"rain", "Rain", "rain tag"},
          Value{"free_electrons", "FreeElectrons", "free electrons tag"},
          Value{"particles", "Particles", "particles tag"},
      }});

  opts.emplace_back(EnumeratedOption{
      .name = "QuantumNumberType",
      .desc = R"(The type of value for a quantum number.
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
      .name = "IsoRatioOption",
      .desc = R"(The type of isotopologue ratio to use.
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
              Value{"T0", "Constant, X0"},
              Value{"T1", "Standard, X0 * (T0/T) ^ X1"},
              Value{"T2", "X0 * (T0/T) ^ X1 * (1 + X2 * log(T/T0))"},
              Value{"T3", "X0 + X1 * (T - T0)"},
              Value{"T4", "(X0 + X1 * (T0/T - 1)) * (T0/T)^X2"},
              Value{"T5", "X0 * (T0/T)^(0.25 + 1.5*X1)"},
              Value{
                  "AER",
                  "X(200) = X0; X(250) = X1; X(298) = X2; X(340) = X3;  Linear interpolation in between"},
              Value{"DPL", "X0 * (T0/T) ^ X1 + X2 * (T0/T) ^ X3"},
              Value{"POLY", "X0 + X1 * T + X2 * T ^ 2 + X3 * T ^ 3 + ..."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeModelCoefficient",
      .desc =
          R"(The type of line shape model coefficients.  See :class:`~pyarts.arts.LineShapeModelType` for more information.
)",
      .values_and_desc =
          {
              Value{"X0", "Constant, X0"},
              Value{"X1", "Constant, X1"},
              Value{"X2", "Constant, X2"},
              Value{"X3", "Constant, X3"},
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
              Value{"ByLine", "Line's are cut 1-by-1"},
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
      .name = "propagation_matrix_agendaPredefined",
      .desc = R"(The types of predefined *propagation_matrix_agenda*.
)",
      .values_and_desc =
          {
              Value{"Empty", R"(Calls:

- *propmat_clearskyInit*
)"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "spectral_radiance_observer_agendaPredefined",
      .desc = R"(The types of predefined *spectral_radiance_observer_agenda*.
)",
      .values_and_desc =
          {
              Value{"Emission", R"(Calls:

- *propagation_path_observer_agendaExecute*
- *spectral_radianceStandardEmission*)"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "spectral_radiance_space_agendaPredefined",
      .desc = R"(The types of predefined *spectral_radiance_space_agenda*.
)",
      .values_and_desc =
          {
              Value{"UniformCosmicBackground", R"(Calls:

- *spectral_radianceUniformCosmicBackground*
- *spectral_radiance_jacobianEmpty*)"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "spectral_radiance_surface_agendaPredefined",
      .desc = R"(The types of predefined *spectral_radiance_surface_agenda*.
)",
      .values_and_desc =
          {
              Value{"Blackbody", R"(Calls:

- *spectral_radianceSurfaceBlackbody*
)"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "propagation_path_observer_agendaPredefined",
      .desc = R"(The types of predefined *propagation_path_observer_agenda*.
)",
      .values_and_desc =
          {
              Value{"Geometric", R"(Calls:

- *propagation_pathGeometric*

  - With ``pos`` as *spectral_radiance_observer_position*
  - With ``los`` as *spectral_radiance_observer_line_of_sight*
  - With ``as_observer`` as ``1``
)"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "HitranType",
      .desc = R"(The type of HITRAN catalog being read.
)",
      .values_and_desc =
          {
              Value{"Pre2004", "2004 version changed the .par-length"},
              Value{"Post2004", "New par length"},
              Value{
                  "Online",
                  "Online expects a modern .par line followed by Upper then Lower quantum numbers"},
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
      .name = "GuiXScaling",
      .desc = R"(X-scaling in the GUI.
)",
      .values_and_desc =
          {
              Value{"Hz", "No scaling"},
              Value{"GHz", "Scaling by 1e-9"},
              Value{"THz", "Scaling by 1e-12"},
              Value{"Angcm", "Scaling to wavenumbers"},
              Value{"Kaycm", "Scale to Kaysers"},
              Value{"m", "Scale to meter wavelengths"},
              Value{"nm", "Scale to nanometer wavelengths"},
              Value{"Angfreq", "Scale to angular frequency"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "GuiPropmatScaling",
      .desc = R"(Propagation matrix scaling in the GUI.
)",
      .values_and_desc =
          {
              Value{"None", "No scaling"},
              Value{"Normalize",
                    "Scaling by normalizing the maximum one value to 1"},
              Value{"CrossSection",
                    "Scaling by outputting in units of cross section"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "GuiTramatScaling",
      .desc = R"(Transmission matrix scaling in the GUI.
)",
      .values_and_desc =
          {
              Value{"None", "No scaling"},
              Value{"dB", "Scaling by showing absorption in decibel"},
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
              Value{"u", "East component"},
              Value{"v", "North component"},
              Value{"w", "Up component"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "SpectralRadianceUnitType",
      .desc = R"(Choice of spectral radiance unit in conversions.
)",
      .values_and_desc =
          {
              Value{"RJBT",
                    R"(Rayleigh-Jeans brightness temperature [K]
)"},
              Value{"PlanckBT",
                    R"(Planck brightness temperature [K]
)"},
              Value{"W_m2_m_sr", "Spectral radiance [W m^-2 m sr^-1]"},
              Value{"W_m2_m1_sr", "Spectral radiance [W m^-2 m^-1 sr^-1]"},
              Value{"unit", "Unit spectral radiance [1]"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "GuiVMR",
      .desc = R"(The type of VMR scaling in the GUI.
)",
      .values_and_desc =
          {
              Value{"exact", "No scaling"},
              Value{"percent", "Scaling by percent"},
              Value{"ppmv", "Scaling by ppmv"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "AbsorptionMirroringTypeOld",
      .desc = R"(The type of mirroring to use.
)",
      .values_and_desc =
          {
              Value{"None", "No mirroring"},
              Value{"Lorentz", "Mirror, but use Lorentz line shape"},
              Value{"SameAsLineShape", "Mirror using the same line shape"},
              Value{
                  "Manual",
                  "Mirror by having a line in the array of line record with negative F0"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "AbsorptionNormalizationTypeOld",
      .desc = R"(The type of normalization to use.
)",
      .values_and_desc =
          {
              Value{"None", "Do not renormalize the line shape"},
              Value{"VVH",
                    "Renormalize with Van Vleck and Huber specifications"},
              Value{"VVW",
                    "Renormalize with Van Vleck and Weiskopf specifications"},
              Value{"RQ",
                    "Renormalize using Rosenkranz's quadratic specifications"},
              Value{
                  "SFS",
                  "Renormalize using simple frequency scaling of the line strength"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "AbsorptionPopulationTypeOld",
      .desc = R"(The type of population to use.
)",
      .values_and_desc =
          {
              Value{"LTE", "Assume band is in LTE"},
              Value{
                  "NLTE",
                  "Assume band is in NLTE and the upper-to-lower ratio is known"},
              Value{
                  "VibTemps",
                  "Assume band is in NLTE described by vibrational temperatures and LTE at other levels"},
              Value{
                  "ByHITRANRosenkranzRelmat",
                  "Assume band needs to compute relaxation matrix to derive HITRAN Y-coefficients"},
              Value{
                  "ByHITRANFullRelmat",
                  "Assume band needs to compute and directly use the relaxation matrix according to HITRAN"},
              Value{
                  "ByMakarovFullRelmat",
                  "Assume band needs to compute and directly use the relaxation matrix according to Makarov et al 2020"},
              Value{
                  "ByRovibLinearDipoleLineMixing",
                  "Assume band needs to compute and directly use the relaxation matrix according to Hartmann, Boulet, Robert, 2008, 1st edition"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "AbsorptionCutoffTypeOld",
      .desc = R"(The type of line shape to use.
  )",
      .values_and_desc =
          {
              Value{"None", "No cutoff frequency at all"},
              Value{
                  "ByLine",
                  "The cutoff frequency is at SingleLine::F0 plus the cutoff frequency plus the speed independent pressure shift"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeTypeOld",
      .desc = R"(The type of line shape to use.
  )",
      .values_and_desc =
          {
              Value{"DP", "Doppler line shape"},
              Value{"LP", "Lorentz line shape"},
              Value{"VP", "Voigt line shape"},
              Value{"SDVP", "Speed-dependent Voigt line shape"},
              Value{"HTP", "Hartmann-Tran line shape"},
              Value{"SplitLP",
                    "Lorentz line shape split by broadening species"},
              Value{"SplitVP", "Voigt line shape split by broadening species"},
              Value{
                  "SplitSDVP",
                  "Speed-dependent Voigt line shape split by broadening species"},
              Value{"SplitHTP",
                    "Hartmann-Tran line shape split by broadening species"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeTemperatureModelOld",
      .desc = R"(The old type of temperature models property.
)",
      .values_and_desc =
          {
              Value{"None", "No temperature model"},
              Value{"T0", "Constant, X0"},
              Value{"T1", "Standard, X0 * (T0/T) ^ X1"},
              Value{"T2", "X0 * (T0/T) ^ X1 * (1 + X2 * log(T/T0))"},
              Value{"T3", "X0 + X1 * (T - T0)"},
              Value{"T4", "(X0 + X1 * (T0/T - 1)) * (T0/T)^X2"},
              Value{"T5", "X0 * (T0/T)^(0.25 + 1.5*X1)"},
              Value{
                  "LM_AER",
                  "X(200) = X0; X(250) = X1; X(298) = X2; X(340) = X3;  Linear interpolation in between"},
              Value{"DPL", "X0 * (T0/T) ^ X1 + X2 * (T0/T) ^ X3"},
              Value{"POLY", "X0 + X1 * T + X2 * T ^ 2 + X3 * T ^ 3"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeVariableOld",
      .desc = "The old type of line shape variables property.\n",
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
          {Value{"Sphere", R"(A spherical Earth. The radius is set following
the value set for the Earth radius.)"},
           Value{"WGS84", R"(The reference ellipsoid used by the GPS system.
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
              Value{
                  "Sphere",
                  "A spherical planetesimal."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "GanymedeEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {
              Value{
                  "Sphere",
                  "A spherical planetesimal."},
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


  return opts;
}

const std::vector<EnumeratedOption>& internal_options() {
  static std::vector<EnumeratedOption> opts = internal_options_create();
  return opts;
}

std::string_view EnumeratedOption::sz() const {
  if (values_and_desc.size() < std::numeric_limits<char>::max()) {
    return "char";
  }
  if (values_and_desc.size() < std::numeric_limits<short>::max()) {
    return "short";
  }
  if (values_and_desc.size() < std::numeric_limits<int>::max()) {
    return "int";
  }
  if (values_and_desc.size() < std::numeric_limits<long>::max()) {
    return "long";
  }

  throw std::runtime_error("Too many values for enum class");
}

std::string EnumeratedOption::docs() const {
  std::ostringstream os;

  const auto n = values_and_desc.front().size();

  os << "Group name: " << std::quoted(name) << "\n\n"
     << desc << "\n\nValid options:\n\n";
  for (auto& v : values_and_desc) {
    std::string_view x = "- ";
    for (auto& s : v | std::views::take(n - 1)) {
      os << std::exchange(x, " or ") << std::quoted(s);
    }
    os << " - " << v.back() << '\n';
  }

  return os.str();
}

std::string EnumeratedOption::tail() const {
  std::ostringstream os;

  const auto m = values_and_desc.size();
  const auto n = values_and_desc.front().size();

  // Create good enum function
  os << "template<> constexpr bool good_enum<" << name << ">(" << name
     << " x) noexcept {\n  const auto v = static_cast<std::size_t>(x);\n  return v < "
     << m << ";\n}\n\n";

  // Create enumtyps array
  os << "namespace enumtyps {\n  static constexpr std::array " << name
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
  os << "  template <int i=0>\n  inline constexpr auto " << name
     << "Names = enum_str_data<" << name << ", i>::strs;\n";
  os << "}  // namespace enumstrs\n\n";

  // Create toString functions
  os << "template <int i=0> constexpr const std::string_view toString(" << name
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
    for (const auto& v : values_and_desc) {
      all.insert(all.end(), v.begin(), v.end());
    }
    std::ranges::sort(all);
    if (auto ptr = std::adjacent_find(all.begin(), all.end());
        ptr != all.end()) {
      throw std::runtime_error("Duplicate value in enum class " + name +
                               std::string{": "} + *ptr);
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
     << " x) {\n  return os << toString(x);\n}\n\n";
  os << "std::istream &operator>>(std::istream &is, " << name
     << "& x) {\n  std::string s;\n  is >> s;\n  x = to<" << name
     << ">(s);\n  return is;\n}\n\n";

  return os.str();
}