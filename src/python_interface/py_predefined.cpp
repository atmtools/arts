#include <predefined/predef.h>
#include <py_auto_interface.h>
#include <pybind11/cast.h>
#include <pybind11/pytypes.h>

#include "predefined/predef_data.h"
#include "py_macros.h"

namespace Python {
void internalCKDMT400(
    py::module_& m,
    py::class_<Absorption::PredefinedModel::MT_CKD400::WaterData>& c) {
  c.def(py::init([]() {
     return std::make_unique<Absorption::PredefinedModel::MT_CKD400::WaterData>();
   }))
      .def_readonly_static(
          "key", &Absorption::PredefinedModel::MT_CKD400::WaterData::key)
      .def_readwrite(
          "ref_temp",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_temp)
      .def_readwrite(
          "ref_press",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_press)
      .def_readwrite(
          "ref_h2o_vmr",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::ref_h2o_vmr)
      .def_readwrite(
          "self_absco_ref",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::self_absco_ref)
      .def_readwrite(
          "for_absco_ref",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::for_absco_ref)
      .def_readwrite(
          "wavenumbers",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::wavenumbers)
      .def_readwrite(
          "self_texp",
          &Absorption::PredefinedModel::MT_CKD400::WaterData::self_texp)
      .def(py::pickle(
          [](const Absorption::PredefinedModel::MT_CKD400::WaterData &t) {
            return py::make_tuple(t.self_absco_ref, t.for_absco_ref,
                                  t.wavenumbers, t.self_texp);
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")
            auto out = std::make_unique<Absorption::PredefinedModel::MT_CKD400::WaterData>();
            out->self_absco_ref = t[0].cast<std::vector<double>>();
            out->for_absco_ref = t[1].cast<std::vector<double>>();
            out->wavenumbers = t[2].cast<std::vector<double>>();
            out->self_texp = t[3].cast<std::vector<double>>();
            return out;
          }));

  m.def(
      "get_foreign_h2o_ckdmt400",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         PredefinedModelData& data) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD400::compute_foreign_h2o(
            pm,
            f,
            p,
            t,
            x,
            data.get<Absorption::PredefinedModel::MT_CKD400::WaterData>());
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::arg("predefined_model_data"),
      py::doc(R"--(Computes foreign absorption using MT CKD Hitran version

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
    predefined_model_data : PredefinedModelData
       As WSV
)--"));

  m.def(
      "get_self_h2o_ckdmt400",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         PredefinedModelData& data) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD400::compute_self_h2o(
            pm,
            f,
            p,
            t,
            x,
            data.get<Absorption::PredefinedModel::MT_CKD400::WaterData>());
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::arg("predefined_model_data"),
      py::doc(R"--(Computes self absorption using MT CKD Hitran version

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
    predefined_model_data : PredefinedModelData
       As WSV
)--"));
}

void internalMPM89(py::module_& m) {
  m.def(
      "get_h2o_mpm89",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MPM89::water(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes water absorption using MPM89

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_o2_mpm89",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MPM89::oxygen(
            pm, f, p, t, x, h2o);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o")=0.0,
      py::doc(R"--(Computes oxygen absorption using MPM89

Parameters:
    f_grid : Vector
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
    x_h2o : Numeric , optional
        Ratio of water in the atmosphere in the range [0, 1]
)--"));
}

void internalELL07(py::module_& m) {
  m.def(
      "get_water_droplet_ell07",
      [](const Vector& f, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::ELL07::compute(pm, f, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_temperature"),
      py::arg("lwc"),
      py::doc(R"--(Computes water absorption using PWR98

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_temperature : Numeric
        Temperature value [K]
    lwc : Numeric
        Liquid water content [1e-10, ...)
)--"));
}

void internalPWR98(py::module_& m) {
  m.def(
      "get_h2o_pwr98",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::PWR98::water(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes water absorption using PWR98

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_o2_pwr98",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::PWR98::oxygen(
            pm, f, p, t, x, h2o);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o")=0.0,
      py::doc(R"--(Computes oxygen absorption using PWR98

Parameters:
    f_grid : Vector
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
    x_h2o : Numeric , optional
        Ratio of water in the atmosphere in the range [0, 1]
)--"));
}



void internalTRE05(py::module_& m) {
  m.def(
      "get_o2_tre05",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::TRE05::oxygen(
            pm, f, p, t, x, h2o);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o")=0.0,
      py::doc(R"--(Computes oxygen absorption using TRE05

Parameters:
    f_grid : Vector
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
    x_h2o : Numeric , optional
        Ratio of water in the atmosphere in the range [0, 1]
)--"));
}

void internalCKDMT100(py::module_& m) {
  m.def(
      "get_o2_cia_ckdmt100",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD100::oxygen_cia(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::doc(R"--(Computes self absorption using MT CKD version 3.50

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_o2_v0v0_ckdmt100",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric n2) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD100::oxygen_v0v0(pm, f, p, t, x, n2);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_n2"),
      py::doc(R"--(Computes self absorption using MT CKD version 3.50

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
    x_n2 : Numeric
        Ratio of nitrogen in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_o2_v1v0_ckdmt100",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD100::oxygen_v0v1(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::doc(R"--(Computes self absorption using MT CKD version 3.50

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
)--"));
}

void internalCKDMT252(py::module_& m) {
  m.def(
      "get_co2_ckdmt252",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::carbon_dioxide(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_co2"),
      py::doc(R"--(Computes self absorption using MT CKD version 2.52

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_co2 : Numeric
        Ratio of carbon dioxide in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_o2_vis_ckdmt252",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::oxygen_vis(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::doc(R"--(Computes self absorption using MT CKD version 2.52

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_n2_fun_ckdmt252",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o, Numeric o2) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::nitrogen_fun(pm, f, p, t, x, h2o, o2);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      py::arg("x_h2o"),
      py::arg("x_o2"),
      py::doc(R"--(Computes N2 fun absorption using MT CKD version 2.52

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_n2 : Numeric
        Ratio of nitrogen in the atmosphere in the range [0, 1]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_n2_rot_ckdmt252",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o, Numeric o2) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::MT_CKD252::nitrogen_rot(pm, f, p, t, x, h2o, o2);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      py::arg("x_h2o"),
      py::arg("x_o2"),
      py::doc(R"--(Computes N2 rot absorption using MT CKD version 2.52

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_n2 : Numeric
        Ratio of nitrogen in the atmosphere in the range [0, 1]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
)--"));
}

void internalSTANDARD(py::module_& m) {
  m.def(
      "get_h2o_self_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::Standard::water_self(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes self water absorption using Rosenkranz standard

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_h2o_foreign_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::Standard::water_foreign(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes foreign water absorption using Rosenkranz standard

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_n2_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::Standard::nitrogen(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_n2"),
      py::doc(R"--(Computes nitrogen absorption using Rosenkranz standard

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_n2 : Numeric
        Ratio of nitrogen in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_o2_standard",
      [](const Vector& f, Numeric p, Numeric t, Numeric x, Numeric h2o) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::Standard::oxygen(pm, f, p, t, x, h2o);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_o2"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes oxygen absorption using Rosenkranz standard

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_o2 : Numeric
        Ratio of oxygen in the atmosphere in the range [0, 1]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));
}

void internalCKDMT350(py::module_& m) {
  m.def(
      "get_self_h2o_ckdmt350",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_self_h2o(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes self absorption using MT CKD version 3.50

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_foreign_h2o_ckdmt350",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_foreign_h2o(
            pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes foreign absorption using MT CKD version 3.50

Parameters:
    f_grid : Vector
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));
}

struct PyInternalClasses {
  py::class_<Absorption::PredefinedModel::MT_CKD400::WaterData> mod1;
};

PyInternalClasses internalNamedClasses(py::module_& predef) {
  auto mtckd400 = predef.def_submodule("MTCKD400");
  py::class_<Absorption::PredefinedModel::MT_CKD400::WaterData>
      mtckd400_water_data(mtckd400, "WaterData");

    return {mtckd400_water_data};
}

void py_predefined(py::module_& m) {
  //! The predef python namespace where all INTERNAL data lives
  auto predef = m.def_submodule("predef");
  predef.doc() = "Contains predefined absorption models";

  //! Must name internal classes first as PredefinedModelData can return them
  auto [mtckd400_water_data] = internalNamedClasses(predef);

  //! A predefined model key
  py::class_<PredefinedModelDataKey>(predef, "PredefinedModelDataKey")
      .def(py::init([]() { return std::make_unique<PredefinedModelDataKey>(); }))
      .def(py::init([](const std::string& c) {
        return Absorption::PredefinedModel::toDataKeyOrThrow(c);
      }))
      .PythonInterfaceCopyValue(PredefinedModelDataKey)
      .PythonInterfaceBasicRepresentation(PredefinedModelDataKey)
      .def(py::pickle(
          [](const PredefinedModelDataKey& t) {
            return py::make_tuple(
                std::string(Absorption::PredefinedModel::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<PredefinedModelDataKey>(
                Absorption::PredefinedModel::toDataKey(
                    t[0].cast<std::string>()));
          }));
  py::implicitly_convertible<std::string, PredefinedModelDataKey>();

  //! ARTS Workspace class, must live on the main (m) namespace
  py::class_<PredefinedModelData>(m, "PredefinedModelData")
      .def(py::init([]() { return std::make_unique<PredefinedModelData>(); }))
      .PythonInterfaceWorkspaceVariableConversion(PredefinedModelData)
      .PythonInterfaceFileIO(PredefinedModelData)
      .PythonInterfaceCopyValue(PredefinedModelData)
      .PythonInterfaceBasicRepresentation(PredefinedModelData)
      .def("set",
           [](PredefinedModelData& x,
              Absorption::PredefinedModel::MT_CKD400::WaterData d) {
             x.set(std::move(d));
           })
      .def("get_h2o_data_mtckd400",
           [](PredefinedModelData& x) {
             return x
                 .get<Absorption::PredefinedModel::MT_CKD400::WaterData>();
           })
      .def(py::pickle(
          [](const PredefinedModelData& t) { return py::make_tuple(t.data); },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            auto out = std::make_unique<PredefinedModelData>();;
            out->data = t[0].cast<PredefinedModelData::DataMap>();
            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(PredefinedModelData);
  
  //! All internal functionality, included methods of named classes go on the predef namespace
  internalCKDMT350(predef);
  internalCKDMT252(predef);
  internalCKDMT100(predef);
  internalCKDMT400(predef, mtckd400_water_data);
  internalMPM89(predef);
  internalPWR98(predef);
  internalTRE05(predef);
  internalELL07(predef);
  internalSTANDARD(predef);
}
}  // namespace Python