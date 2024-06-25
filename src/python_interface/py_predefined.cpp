#include <predefined/predef.h>
#include <python_interface.h>

#include <filesystem>
#include <memory>
#include <unordered_map>

#include "debug.h"
#include "hpy_arts.h"
#include "isotopologues.h"
#include "predefined/predef_data.h"
#include "py_macros.h"
#include "species_tags.h"

namespace Python {
void internalCKDMT400(py::module_& m) {
  py::class_<Absorption::PredefinedModel::MT_CKD400::WaterData>(
      m, "MTCKD400WaterData")
      .def(py::init<Absorption::PredefinedModel::MT_CKD400::WaterData>())
      .def_rw("ref_temp",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_temp,
              ":class:`float` Reference temperature")
      .def_rw("ref_press",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_press,
              ":class:`float` Reference pressure")
      .def_rw("ref_h2o_vmr",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_h2o_vmr,
              ":class:`float` Reference water VMR")
      .def_rw(
          "self_absco_ref",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::self_absco_ref,
          ":class:`list` Self absorption")
      .def_rw("for_absco_ref",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::for_absco_ref,
              ":class:`list` Foreign absorption")
      .def_rw("wavenumbers",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::wavenumbers,
              ":class:`list` Wavenumbers")
      .def_rw("self_texp",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::self_texp,
              ":class:`list` Self temperature exponent")
      .def("__getstate__",
           [](const Absorption::PredefinedModel::MT_CKD400::WaterData& t) {
             return std::make_tuple(
                 t.self_absco_ref, t.for_absco_ref, t.wavenumbers, t.self_texp);
           })
      .def("__setstate__",
           [](Absorption::PredefinedModel::MT_CKD400::WaterData* o,
              const std::tuple<std::vector<double>,
                               std::vector<double>,
                               std::vector<double>,
                               std::vector<double>>& state) {
             new (o) Absorption::PredefinedModel::MT_CKD400::WaterData{};
             o->self_absco_ref = std::get<0>(state);
             o->for_absco_ref  = std::get<1>(state);
             o->wavenumbers    = std::get<2>(state);
             o->self_texp      = std::get<3>(state);
           })
      .doc() = "Water data representation for the MT CKD 4.0 model";

  m.def(
      "get_foreign_h2o_ckdmt400",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         PredefinedModelData& data) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD400::compute_foreign_h2o(
            pm,
            f,
            p,
            t,
            x,
            std::get<Absorption::PredefinedModel::MT_CKD400::WaterData>(
                data.data.at("H2O-ForeignContCKDMT400"_isot)));
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::arg("predefined_model_data"),
      R"--(Computes foreign absorption using MT CKD Hitran version

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]
predefined_model_data : ~pyarts.arts.PredefinedModelData
    As WSV

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_self_h2o_ckdmt400",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         PredefinedModelData& data) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD400::compute_self_h2o(
            pm,
            f,
            p,
            t,
            x,
            std::get<Absorption::PredefinedModel::MT_CKD400::WaterData>(
                data.data.at("H2O-SelfContCKDMT400"_isot)));
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::arg("predefined_model_data"),
      R"--(Computes self absorption using MT CKD Hitran version

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]
predefined_model_data : ~pyarts.arts.PredefinedModelData
    As WSV

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalMPM89(py::module_& m) {
  m.def(
      "get_h2o_mpm89",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MPM89::water(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes water absorption using MPM89

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_mpm89",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MPM89::oxygen(pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o") = 0.0,
      R"--(Computes oxygen absorption using MPM89

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]
x_h2o : float , optional
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalMPM93(py::module_& m) {
  m.def(
      "get_n2_mpm93",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MPM93::nitrogen(pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      py::arg("x_h2o") = 0.0,
      R"--(Computes nitrogen absorption using MPM93

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_n2 : float
    Ratio of nitrogen in the atmosphere in the range [0, 1]
x_h2o : float, optional
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalELL07(py::module_& m) {
  m.def(
      "get_water_droplet_ell07",
      [](const Vector& f, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::ELL07::compute(pm, f, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_temperature"),
      py::arg("lwc"),
      R"--(Computes water absorption using PWR98

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_temperature : float
    Temperature value [K]
lwc : float
    Liquid water content [1e-10, ...)

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalPWR98(py::module_& m) {
  m.def(
      "get_h2o_pwr98",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::PWR98::water(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes water absorption using PWR98

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_pwr98",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::PWR98::oxygen(pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o") = 0.0,
      R"--(Computes oxygen absorption using PWR98

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]
x_h2o : float , optional
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalPWR20xx(py::module_& m) {
  m.def(
      "get_h2o_pwr2021",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::PWR20xx::compute_h2o_2021(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes water absorption using PWR2021

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_pwr2021",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::PWR20xx::compute_o2_2021(
            pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o") = 0.0,
      R"--(Computes oxygen absorption using PWR2021

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]
x_h2o : float , optional
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_h2o_pwr2022",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::PWR20xx::compute_h2o_2022(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes water absorption using PWR2022

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_pwr2022",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::PWR20xx::compute_o2_2022(
            pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o") = 0.0,
      R"--(Computes oxygen absorption using PWR2022

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]
x_h2o : float , optional
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_pwr2021",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::PWR20xx::compute_n2(pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      py::arg("x_h2o") = 0.0,
      R"--(Computes nitrogen absorption using PWR2021

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_n2 : float
    Ratio of nitrogen in the atmosphere in the range [0, 1]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalTRE05(py::module_& m) {
  m.def(
      "get_o2_tre05",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::TRE05::oxygen(pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o") = 0.0,
      R"--(Computes oxygen absorption using TRE05

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]
x_h2o : float , optional
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT100(py::module_& m) {
  m.def(
      "get_o2_cia_ckdmt100",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD100::oxygen_cia(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_v0v0_ckdmt100",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric n2)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD100::oxygen_v0v0(pm, f, p, t, x, n2);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_n2"),
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]
x_n2 : float
    Ratio of nitrogen in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_v1v0_ckdmt100",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD100::oxygen_v0v1(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT252(py::module_& m) {
  m.def(
      "get_co2_ckdmt252",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::carbon_dioxide(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_co2"),
      R"--(Computes self absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_co2 : float
    Ratio of carbon dioxide in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_vis_ckdmt252",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::oxygen_vis(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      R"--(Computes self absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_fun_ckdmt252",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         Numeric h2o,
         Numeric o2) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::nitrogen_fun(
            pm, f, p, t, x, h2o, o2);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      py::arg("x_h2o"),
      py::arg("x_o2"),
      R"--(Computes N2 fun absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_n2 : float
    Ratio of nitrogen in the atmosphere in the range [0, 1]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_rot_ckdmt252",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         Numeric h2o,
         Numeric o2) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::nitrogen_rot(
            pm, f, p, t, x, h2o, o2);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      py::arg("x_h2o"),
      py::arg("x_o2"),
      R"--(Computes N2 rot absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_n2 : float
    Ratio of nitrogen in the atmosphere in the range [0, 1]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalSTANDARD(py::module_& m) {
  m.def(
      "get_h2o_self_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::Standard::water_self(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes self water absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_h2o_foreign_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::Standard::water_foreign(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes foreign water absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::Standard::nitrogen(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      R"--(Computes nitrogen absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_n2 : float
    Ratio of nitrogen in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o)
          -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::Standard::oxygen(pm, f, p, t, x, h2o);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o"),
      R"--(Computes oxygen absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_o2 : float
    Ratio of oxygen in the atmosphere in the range [0, 1]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT350(py::module_& m) {
  m.def(
      "get_self_h2o_ckdmt350",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_self_h2o(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_foreign_h2o_ckdmt350",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_foreign_h2o(
            pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes foreign absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT320(py::module_& m) {
  m.def(
      "get_self_h2o_ckdmt320",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::CKDMT320::compute_self_h2o(pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes self absorption using MT CKD version 3.20

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_foreign_h2o_ckdmt320",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropmatVector pm(f.nelem());
        Absorption::PredefinedModel::CKDMT320::compute_foreign_h2o(
            pm, f, p, t, x);
        Vector out(pm.nelem());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      R"--(Computes foreign absorption using MT CKD version 3.20

Parameters
----------
f_grid : ~pyarts.arts.Vector
    Frequency grid [Hz]
rtp_pressure : float
    Pressure value [Pa]
rtp_temperature : float
    Temperature value [K]
x_h2o : float
    Ratio of water in the atmosphere in the range [0, 1]

Returns
-------
abs_coef : ~pyarts.arts.Vector
    Absorption coefficients
)--");
}
void internalNamedModel(py::module_& m) {
  py::class_<Absorption::PredefinedModel::ModelName>(m, "ModelName")
      .def(py::init<Absorption::PredefinedModel::ModelName>())
      .def("__getstate__",
           [](const Absorption::PredefinedModel::ModelName&) {
             return py::make_tuple();
           })
      .def("__setstate__",
           [](Absorption::PredefinedModel::ModelName* s, const std::tuple<>&) {
             new (s) Absorption::PredefinedModel::ModelName{};
           })
      .doc() = "Named model, contains no data";
}

void py_predefined(py::module_& m) try {
  //! The predef python namespace where all INTERNAL data lives
  auto predef  = m.def_submodule("predef");
  predef.doc() = "Contains predefined absorption models";

  //! ARTS Workspace class, must live on the main (m) namespace
  py::class_<PredefinedModelData> pdmd(m, "PredefinedModelData");
  workspace_group_interface(pdmd);
  pdmd.def_static("fromcatalog",
                  [](const char* const basename,
                     const ArrayOfArrayOfSpeciesTag& specs) {
                    auto out = std::make_shared<PredefinedModelData>();
                    for (auto& spec : specs) {
                      for (auto& mod : spec) {
                        if (mod.type == SpeciesTagType::Predefined) {
                          String filename =
                              (std::filesystem::path(basename) /
                               (mod.Isotopologue().FullName() + ".xml"))
                                  .string();
                          if (find_xml_file_existence(filename)) {
                            PredefinedModelData data;
                            xml_read_from_file(filename, data);
                            std::visit(
                                [&](auto& model) {
                                  out->data[mod.Isotopologue()] = model;
                                },
                                data.data.at(mod.Isotopologue()));
                          } else {
                            out->data[mod.Isotopologue()] =
                                Absorption::PredefinedModel::ModelName{};
                          }
                        }
                      }
                    }
                    return out;
                  })
      .def("__getstate__",
           [](const PredefinedModelData& t) {
             const std::unordered_map<SpeciesIsotope,
                                      Absorption::PredefinedModel::ModelVariant>
                 x{t.data.begin(), t.data.end()};
             return std::make_tuple(x);
           })
      .def("__setstate__",
           [](PredefinedModelData* pm,
              const std::tuple<std::unordered_map<
                  SpeciesIsotope,
                  Absorption::PredefinedModel::ModelVariant>>& state) {
             new (pm) PredefinedModelData{};
             const std::unordered_map<SpeciesIsotope,
                                      Absorption::PredefinedModel::ModelVariant>
                 x = std::get<0>(state);
             for (auto& [k, v] : x) {
               pm->data[k] = v;
             }
           });

  //! All internal functionality, included methods of named classes go on the predef namespace
  internalNamedModel(predef);
  internalCKDMT350(predef);
  internalCKDMT320(predef);
  internalCKDMT252(predef);
  internalCKDMT100(predef);
  internalCKDMT400(predef);
  internalMPM89(predef);
  internalMPM93(predef);
  internalPWR98(predef);
  internalPWR20xx(predef);
  internalTRE05(predef);
  internalELL07(predef);
  internalSTANDARD(predef);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize predefined\n", e.what()));
}
}  // namespace Python
