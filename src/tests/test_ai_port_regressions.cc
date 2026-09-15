#include <workspace.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
void check(bool condition, const char* message) {
  if (not condition) throw std::runtime_error(message);
}
void close(Numeric actual, Numeric expected, const char* message, Numeric tolerance = 1e-12) {
  check(std::isfinite(actual) and std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected)),
        message);
}

void singleton_regridding() {
  using namespace scattering;
  auto                t    = std::make_shared<Vector>(Vector{250.0});
  auto                f    = std::make_shared<Vector>(Vector{94e9});
  auto                za   = std::make_shared<Vector>(Vector{45.0});
  auto                aa   = std::make_shared<Vector>(Vector{0.0});
  auto                scat = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{90.0}));
  ScatteringDataGrids grids(std::make_shared<Vector>(Vector{240.0, 260.0}),
                            std::make_shared<Vector>(Vector{90e9, 100e9}),
                            std::make_shared<Vector>(Vector{30.0, 60.0}),
                            std::make_shared<Vector>(Vector{-30.0, 30.0}),
                            std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{80.0, 100.0})));
  const auto          constant = [](const auto& data, auto expected) {
    for (const auto& x : data | by_elem)
      check(std::abs(x - expected) < 1e-12, "Singleton interpolation must preserve constant coefficients");
  };
  AbsorptionVectorData<Numeric, Format::ARO, Representation::Gridded> absorption(t, f, za);
  std::fill_n(absorption.data_handle(), absorption.size(), 2.0);
  constant(absorption.regrid(grids), 2.0);
  ExtinctionMatrixData<Numeric, Format::ARO, Representation::Gridded> extinction(t, f, za);
  std::fill_n(extinction.data_handle(), extinction.size(), 3.0);
  constant(extinction.regrid(grids), 3.0);
  BackscatterMatrixData<Numeric, Format::ARO> backscatter(t, f, za);
  std::fill_n(backscatter.data_handle(), backscatter.size(), 4.0);
  constant(backscatter.regrid(grids), 4.0);
  PhaseMatrixData<Numeric, Format::ARO, Representation::Gridded> phase(t, f, za, aa, scat);
  std::fill_n(phase.data_handle(), phase.size(), 5.0);
  constant(phase.regrid(grids), 5.0);
#ifndef ARTS_NO_SHTNS
  // Only spectral regridding needs SHTNS; retain the gridded checks above
  // and the remaining regressions in builds without that optional backend.
  PhaseMatrixData<Numeric, Format::ARO, Representation::Spectral> spectral(t, f, za, sht::provider.get_instance(1, 1));
  std::fill_n(spectral.data_handle(), spectral.size(), Complex{6.0, 7.0});
  constant(spectral.regrid(grids), Complex{6.0, 7.0});
#endif
}

void monodisperse_cutoff() {
  const ScatteringSpeciesProperty   density{"cloud", ParticulateProperty::NumberDensity};
  const scattering::MonodispersePSD psd{density, 240.0, 260.0};
  AtmPoint                          point{1e5, 250.0};
  point[density] = 7.0;
  for (Numeric temperature : {230.0, 240.0, 250.0, 260.0, 270.0}) {
    point.temperature = temperature;
    const auto result = psd.evaluate_with_derivatives(point, Vector{1e-4}, 1.0, 3.0);
    const bool active = temperature >= 240.0 and temperature <= 260.0;
    close(result.values[0], active ? 7.0 : 0.0, "PSD cutoff value");
    close(result.derivatives.at(density)[0], active ? 1.0 : 0.0, "PSD cutoff derivative");
  }
}

TelsemAtlas atlas_with_skipped_rows() {
  std::ostringstream data;
  data << "3\n";
  for (Index cell = 1; cell <= 3; ++cell) {
    data << cell;
    for (Index i = 0; i < 7; ++i) data << ' ' << (i % 2 ? 0.7 : 0.9);
    for (Index i = 0; i < 7; ++i) data << " 0.01";
    data << (cell == 2 ? " 0 0\n" : " 1 1\n");
  }
  std::istringstream input(data.str());
  TelsemAtlas        atlas;
  atlas.read(input);
  check(atlas.ndat == 2 and atlas.cellnums.size() == 2, "TELSEM skipped records must be compacted");
  check(atlas.contains(1) and not atlas.contains(2) and atlas.contains(3), "TELSEM compacted lookup");
  close(atlas.channel_emissivity[1][0], 0.9, "TELSEM must retain valid rows after a skipped row");
  return atlas;
}

void source_correction() {
  SourceVector source;
  source.J.resize(1, 2);
  source.J = Stokvec{0.0};
  const ArrayOfAscendingGrid frequencies{AscendingGrid{94e9}, AscendingGrid{94e9}};
  const ArrayOfAtmPoint      atmosphere{AtmPoint{0.0, 250.0}, AtmPoint{1e5, 250.0}};
  const ArrayOfPropmatVector total{PropmatVector{Propmat{}}, PropmatVector{Propmat{2.0}}};
  const ArrayOfPropmatVector scattering{PropmatVector{Propmat{}}, PropmatVector{Propmat{1.0}}};
  const ArrayOfStokvecVector absorption{StokvecVector{Stokvec{}}, StokvecVector{Stokvec{0.5}}};
  spectral_rad_srcvec_pathCorrectScattering(source, total, scattering, absorption, frequencies, atmosphere);
  for (Numeric x : source.J[0, 0]) close(x, 0.0, "Vacuum source correction must remain finite and zero");
  close(source.J[0, 1].I() / planck(94e9, 250.0), -0.25, "Absorbing source correction");

  /* Shapes are rejected where a user supplies them, so the call has to go through
   * the workspace to meet the check.  Calling the method directly, as above, is
   * how ARTS itself calls it, and that path is deliberately not checked. */
  source.J.resize(0, 2);
  Workspace ws{WorkspaceInitialization::Empty};
  ws.set("spectral_rad_srcvec_path", source);
  ws.set("spectral_propmat_path", total);
  ws.set("spectral_propmat_scat_path", scattering);
  ws.set("spectral_absvec_scat_path", absorption);
  ws.set("freq_grid_path", frequencies);
  ws.set("atm_path", atmosphere);

  /* Spelled out rather than braced, because a braced empty argument list picks the
   * overload that sets a workspace variable of that name instead of calling it. */
  const std::vector<std::string>                     no_args{};
  const std::unordered_map<std::string, std::string> no_kwargs{};

  bool rejected = false;
  try {
    Method{"spectral_rad_srcvec_pathCorrectScattering", no_args, no_kwargs}(ws);
  } catch (const std::exception&) { rejected = true; }
  check(rejected, "Source correction must reject mismatched frequency dimensions");
}

void surface_polarization() {
  SurfaceField surface;
  surface.ellipsoid                         = Vector2{6371e3, 6371e3};
  surface[SurfaceKey::t]                    = 280.0;
  surface[SurfacePropertyTag{"wind speed"}] = 5.0;
  surface[SurfacePropertyTag{"salinity"}]   = 0.035;
  PropagationPathPoint point;
  const auto           atlas = atlas_with_skipped_rows();
  const auto [lat, lon]      = atlas.get_coordinates(1);
  point.pos                  = Vector3{0.0, lat, lon};
  point.los                  = Vector2{127.0, 0.0};
  MuelmatVector reflection;
  MuelmatMatrix jac;
  const auto    physical = [](const Muelmat& r) {
    const Numeric rv = r[0, 0] + r[0, 1];
    const Numeric rh = r[0, 0] - r[0, 1];
    close(r[2, 2], std::sqrt(rv * rh), "U reflection must use amplitude product");
    close(r[3, 3], std::sqrt(rv * rh), "V reflection must use amplitude product");
    const Stokvec reflected = r * Stokvec{1.0, 0.0, 1.0, 0.0};
    check(reflected.I() + 1e-14 >= std::hypot(reflected.Q(), reflected.U(), reflected.V()),
          "Surface must preserve physical Stokes vectors");
  };
  spectral_surf_reflTelsem(reflection, jac, AscendingGrid{19.35e9}, surface, point, JacobianTargets{}, atlas, -1.0);
  physical(reflection[0]);
  TessemNN net;
  net.nb_inputs     = 5;
  net.nb_cache      = 1;
  net.nb_outputs    = 1;
  net.b1            = Vector(1, 0.0);
  net.b2            = Vector(1, 0.0);
  net.w1            = Matrix(1, 5, 0.0);
  net.w2            = Matrix(1, 1, 0.0);
  net.x_min         = Vector(5, 0.0);
  net.x_max         = Vector(5, 1000.0);
  net.y_min         = Vector{0.2};
  net.y_max         = Vector{0.2};
  TessemNN vertical = net;
  vertical.y_min    = Vector{0.8};
  vertical.y_max    = Vector{0.8};
  spectral_surf_reflTessem(reflection, jac, AscendingGrid{94e9}, surface, point, JacobianTargets{}, net, vertical);
  physical(reflection[0]);
}

void radar_gates() {
  Workspace ws;
  AtmField  atmosphere;
  atmosphere.top_of_atmosphere = 2000.0;
  atmosphere[AtmKey::p]        = 1e5;
  atmosphere[AtmKey::t]        = 250.0;
  SurfaceField surface;
  surface.ellipsoid = Vector2{6371e3, 6371e3};
  const auto gas    = [] {
    Agenda           agenda{"spectral_propmat_agenda"};
    CallbackOperator cb;
    cb.outputs  = {"spectral_propmat", "spectral_nlte_srcvec", "spectral_propmat_jac", "spectral_nlte_srcvec_jac"};
    cb.callback = CallbackOperator::func_t{[](Workspace& local) {
      local.get<PropmatVector>("spectral_propmat")     = PropmatVector{Propmat{0.0}};
      local.get<StokvecVector>("spectral_nlte_srcvec") = StokvecVector{Stokvec{0.0}};
      local.get<PropmatMatrix>("spectral_propmat_jac").resize(0, 1);
      local.get<StokvecMatrix>("spectral_nlte_srcvec_jac").resize(0, 1);
    }};
    agenda.add(Method{"test gas", Wsv{cb}});
    agenda.finalize(true);
    return agenda;
  }();
  const auto path_agenda = [](const Vector& altitudes) {
    Agenda           agenda{"ray_path_observer_agenda"};
    CallbackOperator cb;
    cb.outputs  = {"ray_path"};
    cb.callback = CallbackOperator::func_t{[altitudes](Workspace& local) {
      auto& path = local.get<ArrayOfPropagationPathPoint>("ray_path");
      path.resize(altitudes.size());
      for (Size i = 0; i < altitudes.size(); ++i) {
        path[i].pos      = Vector3{altitudes[i], 0.0, 0.0};
        path[i].los      = Vector2{0.0, 0.0};
        path[i].pos_type = PathPositionType::atm;
        path[i].los_type = PathPositionType::atm;
      }
    }};
    agenda.add(Method{"test path", Wsv{cb}});
    agenda.finalize(true);
    return agenda;
  };
  ArrayOfScatteringSpecies species;
  constexpr Numeric        extinction = 1e-8;
  species.species.emplace_back(HenyeyGreensteinScatterer{
      ExtSSACallback{[extinction](Numeric, const AtmPoint&) { return std::pair{extinction, 1.0}; }}, 0.0});
  MCAntenna antenna;
  antenna.set_pencil_beam();
  const AscendingGrid bins{100.0, 200.0, 400.0};
  StokvecVector       signal, error;
  const auto          run = [&](const Vector& altitudes) {
    MCRadar(ws,
            signal,
            error,
            atmosphere,
            surface,
            species,
            antenna,
            path_agenda(altitudes),
            gas,
            94e9,
            Vector3{1000.0, 0.0, 0.0},
            Vector2{180.0, 0.0},
            Stokvec{1.0},
            bins,
            42,
            1,
            1,
            0.93,
            "1");
    return signal;
  };
  // The coarse segment midpoint (500 m) is beyond the final gate (400 m).
  // Both gates nevertheless overlap the segment and must receive a return.
  const auto coarse = run(Vector{1000.0, 0.0});
  const auto split  = run(Vector{1000.0, 900.0, 800.0, 600.0, 0.0});
  for (Size i = 0; i < signal.size(); ++i) {
    const Numeric midpoint = std::midpoint(bins[i], bins[i + 1]);
    const Numeric expected = extinction / (4.0 * Constant::pi) * std::exp(-2.0 * extinction * midpoint);
    close(coarse[i].I() / expected, 1.0, "Radar must integrate every overlapped gate");
    close(split[i].I() / coarse[i].I(), 1.0, "Gate splitting must preserve the homogeneous return");
    for (Numeric x : error[i]) close(x, 0.0, "Pencil beam is deterministic");
  }
}
}  // namespace

int main() try {
  singleton_regridding();
  monodisperse_cutoff();
  source_correction();
  surface_polarization();
  radar_gates();
  std::cout << "AI port regression tests passed\n";
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
