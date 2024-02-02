#include "workspace_agenda_creator.h"

#include <string>
#include <unordered_map>

SetWsv::SetWsv(std::string n) : name(std::move(n)) {
  if (auto ind = name.find('='); ind not_eq name.npos) {
    *this = SetWsv(name.substr(0, ind), name.substr(ind + 1));
  }
}

AgendaCreator& AgendaCreator::set(const std::string& name,
                                  WorkspaceGroup auto v) {
  a.add(Method{name, std::move(v)});
  return *this;
}

AgendaCreator& AgendaCreator::add(const std::string& name,
                                  std::vector<SetWsv>&& v) {
  std::vector<std::string> args{};
  std::unordered_map<std::string, std::string> kwargs{};

  for (auto& wsv : v) {
    if (wsv.wsv) {
      a.add(Method{"@" + wsv.name, wsv.wsv.value()});
      kwargs[wsv.name] = "@" + wsv.name;
    } else if (wsv.other) {
      kwargs[wsv.name] = wsv.other.value();
    } else {
      args.push_back(wsv.name);
    }
  }

  a.add(Method{name, args, kwargs});

  return *this;
}

AgendaCreator& AgendaCreator::ignore(const std::string& name) {
  a.add(Method{"Ignore", {name}, {}});

  return *this;
}

Agenda AgendaCreator::finalize() && {
  Agenda ag{std::move(a)};
  ag.finalize(true);
  return ag;
};

ENUMCLASS(propagation_matrix_agendaOption, char, Empty)

Agenda get_propagation_matrix_agenda(const std::string& option) {
  AgendaCreator agenda("propmat_clearsky_agenda");

  using enum propagation_matrix_agendaOption;
  switch (topropagation_matrix_agendaOptionOrThrow(option)) {
    case Empty:
      agenda.add("propmat_clearskyInit");
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

ENUMCLASS(spectral_radiance_observer_agendaOption, char, GeometricEmission)

Agenda get_spectral_radiance_observer_agenda(const std::string& option) {
  AgendaCreator agenda("spectral_radiance_observer_agenda");

  using enum spectral_radiance_observer_agendaOption;
  switch (tospectral_radiance_observer_agendaOptionOrThrow(option)) {
    case GeometricEmission:
      agenda.add("propagation_pathGeometric",
                 SetWsv{"pos", "spectral_radiance_observer_position"},
                 SetWsv{"los", "spectral_radiance_observer_line_of_sight"},
                 SetWsv{"as_observer", Index{1}});
      agenda.add("spectral_radianceStandardEmission");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

ENUMCLASS(spectral_radiance_space_agendaOption, char, UniformCosmicBackground)

Agenda get_spectral_radiance_space_agenda(const std::string& option) {
  AgendaCreator agenda("spectral_radiance_space_agenda");

  using enum spectral_radiance_space_agendaOption;
  switch (tospectral_radiance_space_agendaOptionOrThrow(option)) {
    case UniformCosmicBackground:
      agenda.add("spectral_radianceUniformCosmicBackground");
      agenda.add("spectral_radiance_jacobianEmpty");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

ENUMCLASS(spectral_radiance_surface_agendaOption, char, Blackbody)

Agenda get_spectral_radiance_surface_agenda(const std::string& option) {
  AgendaCreator agenda("spectral_radiance_surface_agenda");

  using enum spectral_radiance_surface_agendaOption;
  switch (tospectral_radiance_surface_agendaOptionOrThrow(option)) {
    case Blackbody:
      agenda.add("spectral_radianceSurfaceBlackbody");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}
