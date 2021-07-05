#include "minimize.h"

namespace Minimize {
int Polynom::operator()(const T4::InputType& p, T4::ValueType& f) const {
  for (Index i=0; i<m_values; i++) {
    f[i] = p[0] - Y[i];
    Numeric x = 1;
    for (int j=1; j<m_inputs; j++) {
      x *= X[i];
      f[i] += p[j] * x;
    }
  }
  return 0;
}

int Polynom::df(const T4::InputType&, T4::JacobianType& J) const {
  for (Index i=0; i<m_values; i++) {
    J(i, 0) = 1;
    Numeric x = 1;
    for (int j=1; j<m_inputs; j++) {
      x *= X[i];
      J(i, j) = x;
    }
  }
  return 0;
}

Polynom::InputType Polynom::x0() const { 
  InputType out(m_inputs);
  for (int j=0; j<m_inputs; j++) {
    out[j] = 1;
  }
  return out;
}

int T4::operator()(const T4::InputType& p, T4::ValueType& f) const {
  for (Index i=0; i<m_values; i++) {
    const Numeric G = T0 / T[i];
    const Numeric GX = std::pow(G, p[2]);
    f[i] = (p[0] + p[1] * (G - 1)) * GX - Y[i];
  }
  return 0;
}

int T4::df(const T4::InputType& p, T4::JacobianType& J) const {
  for (Index i=0; i<m_values; i++) {
    const Numeric G = T0 / T[i];
    const Numeric GX = std::pow(G, p[2]);
    J(i, 0) = GX;
    J(i, 1) = (G - 1) * GX;
    J(i, 2) = (p[0] + p[1] * (G - 1)) * GX * std::log(G);
  }
  return 0;
}

T4::InputType T4::x0() const { 
  const Numeric mean_y = mean(Y);
  InputType out(m_inputs);
  out << mean_y, -0.01*mean_y, mean_y < 0 ? -EXP0 : EXP0;
  return out;
}

int DPL::operator()(const DPL::InputType& p, DPL::ValueType& f) const {
  for (Index i=0; i<m_values; i++) {
    const Numeric G = T0 / T[i];
    const Numeric GX1 = std::pow(G, p[1]);
    const Numeric GX3 = std::pow(G, p[3]);
    f[i] = p[0] * GX1 + p[2] * GX3 - Y[i];
  }
  return 0;
}

int DPL::df(const DPL::InputType& p, DPL::JacobianType& J) const {
  for (Index i=0; i<m_values; i++) {
    const Numeric G = T0 / T[i];
    const Numeric lG = std::log(G);
    const Numeric GX1 = std::pow(G, p[1]);
    const Numeric GX3 = std::pow(G, p[3]);
    J(i, 0) = GX1;
    J(i, 1) = p[0] * GX1 * lG;
    J(i, 2) = GX3;
    J(i, 3) = p[2] * GX3 * lG;
  }
  return 0;
}

DPL::InputType DPL::x0() const { 
  const Numeric mean_y = mean(Y);
  InputType out(m_inputs);
  out << mean_y, mean_y < 0 ? -EXP0 : EXP0, -0.01*mean_y, mean_y < 0 ? EXP0 : -EXP0;
  return out;
}
}
