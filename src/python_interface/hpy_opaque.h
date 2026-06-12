#pragma once

#include <nanobind/nanobind.h>
#include <workspace.h>

NB_MAKE_OPAQUE(Array<lagrange_interp::lag_t<-1, lagrange_interp::grid_identity>>);
NB_MAKE_OPAQUE(Array<lagrange_interp::lag_t<-1, lagrange_interp::loncross>>);
NB_MAKE_OPAQUE(std::unordered_map<std::string, Wsv>);
NB_MAKE_OPAQUE(QuantumLevel);
NB_MAKE_OPAQUE(QuantumState);
NB_MAKE_OPAQUE(FileType);
NB_MAKE_OPAQUE(std::vector<Jacobian::AtmTarget>);
NB_MAKE_OPAQUE(std::vector<Jacobian::SurfaceTarget>);
NB_MAKE_OPAQUE(std::vector<Jacobian::LineTarget>);
NB_MAKE_OPAQUE(std::vector<Jacobian::SensorTarget>);
NB_MAKE_OPAQUE(std::vector<Jacobian::ErrorTarget>);
NB_MAKE_OPAQUE(lbl::line_shape::species_model::map_t);
NB_MAKE_OPAQUE(lbl::line_shape::model::map_t);
