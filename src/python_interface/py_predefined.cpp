#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <predefined/predef.h>
#include <python_interface.h>

#include <filesystem>
#include <memory>
#include <unordered_map>

#include "debug.h"
#include "hpy_arts.h"
#include "isotopologues.h"
#include "predefined/predef_data.h"
#include "species_tags.h"

namespace Python {
void internalCKDMT400(py::module_& m) {
  py::class_<Absorption::PredefinedModel::MT_CKD400::WaterData> mm(
      m, "MTCKD400WaterData");
  generic_interface(mm);
  mm.def(py::init<Absorption::PredefinedModel::MT_CKD400::WaterData>())
      .def_rw("ref_temp",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_temp,
              "Reference temperature\n\n.. :class:`Numeric`")
      .def_rw("ref_press",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_press,
              "Reference pressure\n\n.. :class:`Numeric`")
      .def_rw("ref_h2o_vmr",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_h2o_vmr,
              "Reference water VMR\n\n.. :class:`Numeric`")
      .def_rw(
          "self_absco_ref",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::self_absco_ref,
          "Self absorption\n\n.. :class:`Vector`")
      .def_rw("for_absco_ref",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::for_absco_ref,
              "Foreign absorption\n\n.. :class:`Vector`")
      .def_rw("wavenumbers",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::wavenumbers,
              "Wavenumbers\n\n.. :class:`Vector`")
      .def_rw("self_texp",
              &Absorption::PredefinedModel::MT_CKD400::WaterData::self_texp,
              "Self temperature exponent\n\n.. :class:`Vector`")

      .doc() = "Water data representation for the MT CKD 4.0 model";

  m.def(
      "get_foreign_h2o_ckdmt400",
      [](const Vector& f,
         const AtmPoint& atm,
         PredefinedModelData& data) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD400::compute_foreign_h2o(
            pm,
            f,
            atm,
            std::get<Absorption::PredefinedModel::MT_CKD400::WaterData>(
                data.at("H2O-ForeignContCKDMT400"_isot).data));
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      "predefined_model_data"_a,
      R"--(Computes foreign absorption using MT CKD Hitran version

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state
predefined_model_data : ~pyarts3.arts.PredefinedModelData
    As WSV

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_self_h2o_ckdmt400",
      [](const Vector& f,
         const AtmPoint& atm,
         PredefinedModelData& data) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD400::compute_self_h2o(
            pm,
            f,
            atm,
            std::get<Absorption::PredefinedModel::MT_CKD400::WaterData>(
                data.at("H2O-SelfContCKDMT400"_isot).data));
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      "predefined_model_data"_a,
      R"--(Computes self absorption using MT CKD Hitran version

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state
predefined_model_data : ~pyarts3.arts.PredefinedModelData
    As WSV

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalMPM89(py::module_& m) {
  m.def(
      "get_h2o_mpm89",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MPM89::water(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes water absorption using MPM89

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_mpm89",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MPM89::oxygen(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes oxygen absorption using MPM89

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalMPM93(py::module_& m) {
  m.def(
      "get_n2_mpm93",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MPM93::nitrogen(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes nitrogen absorption using MPM93

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalELL07(py::module_& m) {
  m.def(
      "get_water_droplet_ell07",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::ELL07::compute(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes water absorption using PWR98

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalPWR98(py::module_& m) {
  m.def(
      "get_h2o_pwr98",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::PWR98::water(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes water absorption using PWR98

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_pwr98",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::PWR98::oxygen(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes oxygen absorption using PWR98

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalPWR20xx(py::module_& m) {
  m.def(
      "get_h2o_pwr2021",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::PWR20xx::compute_h2o_2021(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes water absorption using PWR2021

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_pwr2021",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::PWR20xx::compute_o2_2021(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes oxygen absorption using PWR2021

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_h2o_pwr2022",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::PWR20xx::compute_h2o_2022(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes water absorption using PWR2022

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_pwr2022",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::PWR20xx::compute_o2_2022(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes oxygen absorption using PWR2022

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_pwr2021",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::PWR20xx::compute_n2(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes nitrogen absorption using PWR2021

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalTRE05(py::module_& m) {
  m.def(
      "get_o2_tre05",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::TRE05::oxygen(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes oxygen absorption using TRE05

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT100(py::module_& m) {
  m.def(
      "get_o2_cia_ckdmt100",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD100::oxygen_cia(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_v0v0_ckdmt100",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD100::oxygen_v0v0(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_v1v0_ckdmt100",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD100::oxygen_v0v1(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT252(py::module_& m) {
  m.def(
      "get_co2_ckdmt252",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD252::carbon_dioxide(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_vis_ckdmt252",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD252::oxygen_vis(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_fun_ckdmt252",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD252::nitrogen_fun(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes N2 fun absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_rot_ckdmt252",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::MT_CKD252::nitrogen_rot(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes N2 rot absorption using MT CKD version 2.52

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalSTANDARD(py::module_& m) {
  m.def(
      "get_h2o_self_standard",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::Standard::water_self(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self water absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_h2o_foreign_standard",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::Standard::water_foreign(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes foreign water absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_n2_standard",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::Standard::nitrogen(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes nitrogen absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_o2_standard",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::Standard::oxygen(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes oxygen absorption using Rosenkranz standard

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT350(py::module_& m) {
  m.def(
      "get_self_h2o_ckdmt350",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::CKDMT350::compute_self_h2o(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_foreign_h2o_ckdmt350",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::CKDMT350::compute_foreign_h2o(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes foreign absorption using MT CKD version 3.50

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}

void internalCKDMT320(py::module_& m) {
  m.def(
      "get_self_h2o_ckdmt320",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::CKDMT320::compute_self_h2o(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes self absorption using MT CKD version 3.20

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");

  m.def(
      "get_foreign_h2o_ckdmt320",
      [](const Vector& f, const AtmPoint& atm) -> Vector {
        PropmatVector pm(f.size());
        Absorption::PredefinedModel::CKDMT320::compute_foreign_h2o(pm, f, atm);
        Vector out(pm.size());
        std::transform(pm.begin(), pm.end(), out.begin(), [](auto& prop) {
          return prop.A();
        });
        return out;
      },
      "f_grid"_a,
      "atm"_a,
      R"--(Computes foreign absorption using MT CKD version 3.20

Parameters
----------
f_grid : ~pyarts3.arts.Vector
    Frequency grid [Hz]
atm : AtmPoint
    The atmospheric state

Returns
-------
abs_coef : ~pyarts3.arts.Vector
    Absorption coefficients
)--");
}
void internalNamedModel(py::module_& m) {
  py::class_<Absorption::PredefinedModel::ModelName> mm(m, "ModelName");
  generic_interface(mm);
  mm.def(py::init<Absorption::PredefinedModel::ModelName>())

      .doc() = "Named model, contains no data";
}

void py_predefined(py::module_& m) try {
  //! The predef python namespace where all INTERNAL data lives
  auto predef  = m.def_submodule("predef");
  predef.doc() = "Contains predefined absorption models";

  py::class_<PredefinedModelDataVariant> var(m, "PredefinedModelDataVariant");
  var.def_rw(
      "data",
      &PredefinedModelDataVariant::data,
      "The data\n\n.. :class:`~pyarts3.arts.predef.ModelName`\n\n.. :class:`~pyarts3.arts.predef.MTCKD400WaterData`");
  generic_interface(var);

  //! ARTS Workspace class, must live on the main (m) namespace
  auto pdmd = py::bind_map<PredefinedModelData>(m, "PredefinedModelData");
  generic_interface(pdmd);
  pdmd.def_static(
          "fromcatalog",
          [](const char* const basename, const ArrayOfSpeciesTag& specs) {
            auto out = PredefinedModelData{};
            for (auto& mod : specs) {
              if (mod.type == SpeciesTagType::Predefined) {
                String filename = (std::filesystem::path(basename) /
                                   (mod.Isotopologue().FullName() + ".xml"))
                                      .string();
                if (find_xml_file_existence(filename)) {
                  PredefinedModelData data;
                  xml_read_from_file(filename, data);
                  std::visit(
                      [&](auto& model) {
                        out[mod.Isotopologue()].data = model;
                      },
                      data.at(mod.Isotopologue()).data);
                } else {
                  out[mod.Isotopologue()].data =
                      Absorption::PredefinedModel::ModelName{};
                }
              }
            }
            return out;
          },
          "basename"_a,
          "specs"_a,
          "Reads predefined models from catalog")
      .def(
          "spectral_propmat",
          [](const PredefinedModelData& self,
             const AscendingGrid& f,
             const AtmPoint& atm,
             const SpeciesEnum& spec,
             const py::kwargs&) {
            PropmatVector spectral_propmat(f.size());
            PropmatMatrix spectral_propmat_jac(0, f.size());
            JacobianTargets jac_targets{};

            spectral_propmatAddPredefined(spectral_propmat,
                                          spectral_propmat_jac,
                                          self,
                                          spec,
                                          jac_targets,
                                          f,
                                          atm);

            return spectral_propmat;
          },
          "f"_a,
          "atm"_a,
          "spec"_a   = SpeciesEnum::Bath,
          "kwargs"_a = py::kwargs{},
          R"--(Computes the predefined model absorption in 1/m

Parameters
----------
f : AscendingGrid
    Frequency grid [Hz]
atm : AtmPoint
    Atmospheric point
    spec : SpeciesEnum, optional
    Species to use

Returns
-------
spectral_propmat : PropmatVector
    Propagation matrix by frequency [1/m]

)--");

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
      std::format("DEV ERROR:\nCannot initialize predefined\n{}", e.what()));
}
}  // namespace Python
