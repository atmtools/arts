#include "lbl_temperature_model.h"

#include <debug.h>
#include <double_imanip.h>
#include <lbl_data.h>

#include <limits>
#include <utility>

#include "arts_conversions.h"

namespace lbl::temperature {
static_assert(LineShapeModelTypeSize == enumsize::LineShapeModelTypeSize,
              "Invalid number of models");

data::data(LineShapeModelType type, Vector X) : t(type), x(std::move(X)) {
  ARTS_USER_ERROR_IF(not good_enum(t), "Invalid model type")
  ARTS_USER_ERROR_IF(Size n = model_size(t);
                     n != std::numeric_limits<Size>::max() and
                         n != static_cast<Size>(x.size()),
                     "Invalid number of parameters for model {}"
                     "\nExpected {}"
                     " parameters but got {} parameters.",
                     t,
                     n,
                     x.size())
}

namespace model {
using std::log;
using std::pow;

Numeric T1(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1);
}

Numeric dT1_dX0(Numeric, Numeric X1, Numeric T0, Numeric T) {
  return pow(T0 / T, X1);
}

Numeric dT1_dX1(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * log(T0 / T);
}

Numeric dT1_dT0(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * X1 * pow(T0 / T, X1) / T0;
}

Numeric dT1_dT(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return -X0 * X1 * pow(T0 / T, X1) / T;
}

Numeric T2(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * (1 + X2 * log(T / T0));
}

Numeric dT2_dX0(Numeric, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X1) * (1 + X2 * log(T / T0));
}

Numeric dT2_dX1(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) * log(T0 / T);
}

Numeric dT2_dX2(Numeric X0, Numeric X1, Numeric, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * log(T / T0);
}

Numeric dT2_dT0(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X0 * X1 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) / T0 -
         X0 * X2 * pow(T0 / T, X1) / T0;
}

Numeric dT2_dT(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return -X0 * X1 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) / T +
         X0 * X2 * pow(T0 / T, X1) / T;
}

Numeric T4(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return (X0 + X1 * (T0 / T - 1)) * pow(T0 / T, X2);
}

Numeric dT4_dX0(Numeric, Numeric, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X2);
}

Numeric dT4_dX1(Numeric, Numeric, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X2) * (T0 / T - 1.);
}

Numeric dT4_dX2(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1)) * log(T0 / T);
}

Numeric dT4_dT0(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X2 * pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1.)) / T0 +
         X1 * pow(T0 / T, X2) / T;
}

Numeric dT4_dT(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return -X2 * pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1.)) / T -
         T0 * X1 * pow(T0 / T, X2) / (T * T);
}

Numeric T5(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, 0.25 + 1.5 * X1);
}

Numeric dT5_dX0(Numeric, Numeric X1, Numeric T0, Numeric T) {
  return pow(T0 / T, 1.5 * X1 + 0.25);
}

Numeric dT5_dX1(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return 1.5 * X0 * pow(T0 / T, 1.5 * X1 + 0.25) * log(T0 / T);
}

Numeric dT5_dT0(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, 1.5 * X1 + 0.25) * (1.5 * X1 + 0.25) / T0;
}

Numeric dT5_dT(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return -X0 * pow(T0 / T, 1.5 * X1 + 0.25) * (1.5 * X1 + 0.25) / T;
}

Numeric DPL(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) + X2 * pow(T0 / T, X3);
}

Numeric dDPL_dX0(Numeric, Numeric X1, Numeric, Numeric, Numeric T0, Numeric T) {
  return pow(T0 / T, X1);
}

Numeric dDPL_dX1(
    Numeric X0, Numeric X1, Numeric, Numeric, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * log(T0 / T);
}

Numeric dDPL_dX2(Numeric, Numeric, Numeric, Numeric X3, Numeric T0, Numeric T) {
  return pow(T0 / T, X3);
}

Numeric dDPL_dX3(
    Numeric, Numeric, Numeric X2, Numeric X3, Numeric T0, Numeric T) {
  return X2 * pow(T0 / T, X3) * log(T0 / T);
}

Numeric dDPL_dT0(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T) {
  return X0 * X1 * pow(T0 / T, X1) / T0 + X2 * X3 * pow(T0 / T, X3) / T0;
}

Numeric dDPL_dT(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T) {
  return -X0 * X1 * pow(T0 / T, X1) / T + -X2 * X3 * pow(T0 / T, X3) / T;
}

Numeric POLY(const ConstVectorView& x, Numeric T) {
  Numeric poly_fac = 1.0;
  Numeric poly_sum = 0.0;
  for (auto X : x) {
    poly_sum += X * poly_fac;
    poly_fac *= T;
  }
  return poly_sum;
}

Numeric dPOLY_dT(const ConstVectorView& x, Numeric T) {
  Numeric poly_fac = 1.0;
  Numeric poly_sum = 0.0;
  for (Size i = 1; i < x.size(); ++i) {
    poly_sum += static_cast<Numeric>(i) * x[i] * poly_fac;
    poly_fac *= T;
  }
  return poly_sum;
}
}  // namespace model

Numeric data::operator()(Numeric T0, Numeric T) const {
#define SWITCHCASE(mod)         \
  case LineShapeModelType::mod: \
    return operator()<LineShapeModelType::mod>(T0, T)

  switch (t) {
    SWITCHCASE(T0);
    SWITCHCASE(T1);
    SWITCHCASE(T2);
    SWITCHCASE(T3);
    SWITCHCASE(T4);
    SWITCHCASE(T5);
    SWITCHCASE(AER);
    SWITCHCASE(DPL);
    SWITCHCASE(POLY);
  }

  return NAN;

#undef SWITCHCASE
}

#define SWITCHCASE(name, mod) \
  case LineShapeModelType::mod: return d##name<LineShapeModelType::mod>(T0, T)

#define DERIVATIVES(name)                              \
  Numeric data::d##name(Numeric T0, Numeric T) const { \
    switch (t) {                                       \
      SWITCHCASE(name, T0);                            \
      SWITCHCASE(name, T1);                            \
      SWITCHCASE(name, T2);                            \
      SWITCHCASE(name, T3);                            \
      SWITCHCASE(name, T4);                            \
      SWITCHCASE(name, T5);                            \
      SWITCHCASE(name, AER);                           \
      SWITCHCASE(name, DPL);                           \
      SWITCHCASE(name, POLY);                          \
    }                                                  \
                                                       \
    return std::numeric_limits<Numeric>::quiet_NaN();  \
  }

DERIVATIVES(X0)
DERIVATIVES(X1)
DERIVATIVES(X2)
DERIVATIVES(X3)
DERIVATIVES(T0)
DERIVATIVES(T)

#undef DERIVATIVES
#undef SWITCHCASE

std::istream& operator>>(std::istream& is, temperature::data& x) try {
  String name;
  is >> name;
  x.t = to<LineShapeModelType>(name);

  Size n = model_size(x.t);
  if (n == std::numeric_limits<Size>::max()) {
    is >> n;
  }
  x.x.resize(n);

  for (auto& v : x.x) is >> double_imanip() >> v;
  return is;
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading LineShapeTemperatureModel:\n{}", e.what()));
}

LineShapeModelType data::Type() const { return t; }

const Vector& data::X() const { return x; }

Numeric& data::X(LineShapeModelCoefficient coeff) {
  switch (coeff) {
    case LineShapeModelCoefficient::X0:
      ARTS_USER_ERROR_IF(x.size() < 1, "No X0 value in data.");
      return x[0];
    case LineShapeModelCoefficient::X1:
      ARTS_USER_ERROR_IF(x.size() < 2, "No X1 value in data.");
      return x[1];
    case LineShapeModelCoefficient::X2:
      ARTS_USER_ERROR_IF(x.size() < 3, "No X2 value in data.");
      return x[2];
    case LineShapeModelCoefficient::X3:
      ARTS_USER_ERROR_IF(x.size() < 4, "No X3 value in data.");
      return x[3];
  }
  std::unreachable();
}

const Numeric& data::X(LineShapeModelCoefficient coeff) const {
  return const_cast<data*>(this)->X(coeff);
}

bool data::is_zero() const {
  switch (t) {
    using enum LineShapeModelType;
    case T0:  [[fallthrough]];
    case T1:  [[fallthrough]];
    case T2:  [[fallthrough]];
    case T5:  return x[0] == 0;
    case T3:  [[fallthrough]];
    case T4:  return x[0] == 0 and x[1] == 0;
    case AER: [[fallthrough]];
    case POLY:
      for (auto y : x)
        if (y != 0) return false;
      return true;
    case DPL: return x[0] == 0 and x[2] == 0;
  }
  return true;
}
}  // namespace lbl::temperature

void xml_io_stream<lbl::temperature::data>::write(
    std::ostream& os,
    const lbl::temperature::data& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.Type(), pbofs);
  xml_write_to_stream(os, x.X(), pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<lbl::temperature::data>::read(std::istream& is,
                                                 lbl::temperature::data& x,
                                                 bifstream* pbifs) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  LineShapeModelType t;
  Vector v;

  xml_read_from_stream(is, t, pbifs);
  xml_read_from_stream(is, v, pbifs);

  x = lbl::temperature::data{t, v};

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading {}:\n{}", type_name, e.what()));
}

namespace {
std::string metric_unit_x0(const Numeric x,
                           const std::optional<LineShapeModelVariable>& type) {
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;

  if (type) {
    const auto [c, v] = Conversion::metric_prefix(x);
    switch (*type) {
      case G0:
      case D0:
      case G2:
      case D2:
      case FVC:
        return std::format("{:.2f}{}{}Hz/Pa", v, c == ' ' ? ""sv : " "sv, c);
      case ETA: return std::format("{}", x);
      case DV:
        return std::format("{:.2f}{}{}Hz/Pa^2", v, c == ' ' ? ""sv : " "sv, c);
      case G:
        return std::format("{:.2f}{}{}/Pa^2", v, c == ' ' ? ""sv : " "sv, c);
      case Y:
        return std::format("{:.2f}{}{}/Pa", v, c == ' ' ? ""sv : " "sv, c);
      case unused: break;
    }
  }

  return std::format("{}", x);
}

std::string metric_unit_x1(const Numeric x,
                           LineShapeModelType t,
                           const std::optional<LineShapeModelVariable>& type) {
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;

  if (type) {
    const auto [c, v] = Conversion::metric_prefix(x);
    switch (*type) {
      case G0:
      case D0:
      case G2:
      case D2:
      case FVC:
        switch (t) {
          case T0:
          case T1:
          case T2:
          case T5:
          case DPL:
          case POLY: return std::format("{}", x);
          case T3:
          case T4:
          case AER:
            return std::format(
                "{:.2f}{}{}Hz/Pa", v, c == ' ' ? ""sv : " "sv, c);
        }
        break;
      case ETA: return std::format("{}", x);
      case DV:
        switch (t) {
          case T0:
          case T1:
          case T2:
          case T5:
          case DPL:
          case POLY: return std::format("{}", x);
          case T3:
            return std::format(
                "{:.2f}{}{}Hz/Pa^2/K", v, c == ' ' ? ""sv : " "sv, c);
          case T4:
          case AER:
            return std::format(
                "{:.2f}{}{}Hz/Pa^2", v, c == ' ' ? ""sv : " "sv, c);
        }
        break;
      case G:
        switch (t) {
          case T0:
          case T1:
          case T2:
          case T5:
          case DPL:
          case POLY: return std::format("{}", x);
          case T3:
            return std::format(
                "{:.2f}{}{}/Pa^2/K", v, c == ' ' ? ""sv : " "sv, c);
          case T4:
          case AER:
            return std::format(
                "{:.2f}{}{}/Pa^2", v, c == ' ' ? ""sv : " "sv, c);
        }
        break;
      case Y:
        switch (t) {
          case T0:
          case T1:
          case T2:
          case T5:
          case DPL:
          case POLY: return std::format("{}", x);
          case T3:
            return std::format(
                "{:.2f}{}{}/Pa/K", v, c == ' ' ? ""sv : " "sv, c);
          case T4:
          case AER:
            return std::format("{:.2f}{}{}/Pa", v, c == ' ' ? ""sv : " "sv, c);
        }
      case unused: break;
    }
  }

  return std::format("{}", x);
}
}  // namespace

std::string to_educational_string(const lbl::temperature::data& data,
                                  std::optional<LineShapeModelVariable> type) {
  switch (data.Type()) {
    using enum LineShapeModelType;
    using enum LineShapeModelCoefficient;
    case T0: return std::format("{}", metric_unit_x0(data.X(X0), type));
    case T1:
      return std::format("{} * pow(T0 / T, {})",
                         metric_unit_x0(data.X(X0), type),
                         metric_unit_x1(data.X(X1), data.Type(), type));
    case T2:
      return std::format("{} * pow(T0 / T, {}) * (1 + {} * log(T / T0))",
                         metric_unit_x0(data.X(X0), type),
                         metric_unit_x1(data.X(X1), data.Type(), type),
                         data.X(X2));
    case T3:
      return std::format("{} + {} * (T - T0)",
                         metric_unit_x0(data.X(X0), type),
                         metric_unit_x1(data.X(X1), data.Type(), type));
    case T4:
      return std::format("({0} + {1} * (T0 / T - 1)) * pow(T0 / T, {2})",
                         metric_unit_x0(data.X(X0), type),
                         metric_unit_x1(data.X(X1), data.Type(), type),
                         data.X(X2));
    case T5:
      return std::format("{} * pow(T0 / T, 0.25 + 1.5 * {})",
                         metric_unit_x0(data.X(X0), type),
                         metric_unit_x1(data.X(X1), data.Type(), type));
    case AER:
      return std::format(
          "(T < 250.0) ? "
          "{0} + (T - 200.0) * ({1} - {0}) / (250.0 - 200.0) : "
          "(T > 296.0) ? "
          "{2} + (T - 296.0) * ({3} - {2}) / (340.0 - 296.0) : "
          "{1} + (T - 250.0) * ({2} - {1}) / (296.0 - 250.0)",
          metric_unit_x0(data.X(X0), type),
          metric_unit_x0(data.X(X1), type),
          metric_unit_x0(data.X(X2), type),
          metric_unit_x0(data.X(X3), type));
    case DPL:
      return std::format("{} * pow(T0 / T, {}) + {} * pow(T0 / T, {})",
                         metric_unit_x0(data.X(X0), type),
                         metric_unit_x1(data.X(X1), data.Type(), type),
                         metric_unit_x0(data.X(X2), type),
                         metric_unit_x1(data.X(X3), data.Type(), type));
    case POLY: {
      std::string s{};
      Size i{};
      std::string res{""};
      std::string_view x = ""sv;
      for (auto& v : data.X()) {
        res += std::format("{}{}{}{}",
                           std::exchange(x, " + "sv),
                           metric_unit_x0(v, type),
                           i == 0   ? ""s
                           : i == 1 ? "/K"s
                                    : std::format("/K^{}", i),
                           s);
        s   += " * T";
        i++;
      }
      return res;
    }
  }
  std::unreachable();
}
