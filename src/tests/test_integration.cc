/*!
  \file   test_integration.cc
  \author Claas Teichmann <claas@sat.physik.uni-bremen.de>
  \date   2003/05/27
  
  \brief  Testfile for the AngIntegrate_trapezoid function from math_funcs.cc
  
*/

#include <matpack.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "array.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "test_perf.h"

inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);
inline constexpr Numeric PI      = Constant::pi;

void init_xy(float stepsize,
             int frequency,
             Matrix& Integrand,
             Vector& za_grid,
             Vector& aa_grid);

void init_x(
    int vsize, float stepsize, int frequency, Vector& Integrand, Vector& Theta);

Numeric AngIntegrate_trapezoid_original(MatrixView Integrand,
                                        ConstVectorView za_grid,
                                        ConstVectorView aa_grid);

Numeric AngIntegrate_trapezoid_opt(MatrixView Integrand,
                                   ConstVectorView za_grid,
                                   ConstVectorView aa_grid);

Numeric AngIntegrate_trapezoid_fixedstep(MatrixView Integrand,
                                         ConstVectorView za_grid,
                                         ConstVectorView aa_grid,
                                         Numeric stepsize);

Numeric AngIntegrate_trapezoid_fixedstep_opt(MatrixView Integrand,
                                             ConstVectorView za_grid,
                                             ConstVectorView aa_grid,
                                             Numeric stepsize);

Numeric AngIntegrate_trapezoid_fixedstep_opt2(MatrixView Integrand,
                                              ConstVectorView za_grid,
                                              ConstVectorView aa_grid,
                                              Numeric stepsize);

Numeric AngIntegrate_trapezoid_original(ConstVectorView Integrand,
                                        ConstVectorView za_grid);

Numeric AngIntegrate_trapezoid_fixedstep(ConstVectorView Integrand,
                                         ConstVectorView za_grid,
                                         Numeric stepsize);

Numeric test_xy(int z_size, int a_size, float stepsize, int frequency);

Numeric test_xy_opt(int z_size, int a_size, float stepsize, int frequency);

Numeric test_xy_fixedstep(int z_size,
                          int a_size,
                          float stepsize,
                          int frequency);

Numeric test_xy_fixedstep_opt(int z_size,
                              int a_size,
                              float stepsize,
                              int frequency);

Numeric test_xy_fixedstep_opt2(int z_size,
                               int a_size,
                               float stepsize,
                               int frequency);

Numeric test_AngIntegrate_trapezoid_opti(int z_size,
                                         int a_size,
                                         float stepsize,
                                         int frequency);

Numeric test_x(int vsize, float stepsize, int frequency);

Numeric test_x_fixedstep(int vsize, int frequency);

int main(int argc, char* argv[]) {
  if (argc == 1) {
    std::cerr << argv[0] << " requires one parameter" << '\n';
    exit(1);
  }

  std::cout << "Uebergabewert von argv : " << argv[1] << '\n';

  int frequency = (int)strtol(argv[1], NULL, 10);  // Zahl zur Basis 10
  std::cout << "Wert von frequency    : " << frequency << '\n';

  //  test_x(1801, 0.1, frequency);
  //  test_x_fixedstep(1801, frequency);
  //  test_x(181, 1, frequency);
  //  test_x(19, 10, frequency);

  test_xy(181, 361, 1.0, frequency);
  //  test_xy_opt(1801, 3601, 0.1, frequency);
  //  test_xy_fixedstep(1801, 3601, 0.1, frequency);
  //  test_xy_fixedstep_opt(1801, 3601, 0.1, frequency);
  //  test_xy_fixedstep_opt(181, 361, 1.0, frequency);
  //  test_xy_fixedstep_opt(19, 37, 10, frequency);
  //  test_xy_fixedstep_opt2(1801, 3601, 0.1, frequency);
  //  test_xy_fixedstep_opt2(181, 361, 1.0, frequency);
  test_AngIntegrate_trapezoid_opti(181, 361, 1.0, frequency);
  //  test_xy_fixedstep_opt2(19, 37, 10, frequency);
}

//! init_x
/*!
  This function fills a Vector with funcition values
  The function is y=x
 */
void init_x(int vsize,
            float stepsize,
            int frequency,
            Vector& Integrand,
            Vector& Theta) {
  std::cout << "----------------init_x---------------\n";

  //  Integrand.resize(vsize); // function to be integrated
  //  Theta.resize(vsize);     // Theta values

  for (Size i = 0; i < Integrand.size(); i++) Integrand[i] = (float)i * stepsize;

  //Theta is between 0 and 180
  for (Size i = 0; i < Theta.size(); i++) Theta[i] = (float)i * stepsize;

  std::cout << "function Y = X" << '\n'
            << "vsize = " << vsize << '\n'
            << "stepsize = " << stepsize << '\n'
            << "frequency = " << frequency << '\n'
            << "Integrand: von " << Integrand[0] << " bis "
            << Integrand[Integrand.size() - 1] << '\n'
            << "Theta: von " << Theta[0] << " bis " << Theta[Theta.size() - 1]
            << '\n';
}

//! init_xy
/*!
  This function fills Matrix with funcition values
  The funcion is a circle around 0
 */
void init_xy(float stepsize,
             int frequency,
             Matrix& Integrand,
             Vector& za_grid,
             Vector& aa_grid) {
  std::cout << ">>>>>-----------init_xy---------------\n";
  Index n_za = za_grid.size();
  Index n_aa = aa_grid.size();

  // The r=1 so we get a circle
  for (Index i = 0; i < n_za; i++)
    for (Index j = 0; j < n_aa; j++) Integrand[i, j] = 1;

  //za_grid (Theta) is between 0 and 180
  for (Index i = 0; i < n_za; i++) za_grid[i] = (float)i * stepsize;

  //aa_grid (Phi) is between 0 and 360
  for (Index i = 0; i < n_aa; i++) aa_grid[i] = (float)i * stepsize;

  std::cout << "function x^2 + y^2 + z^2 = 1" << '\n'
            << "n_za = " << n_za << '\n'
            << "n_aa = " << n_aa << '\n'
            << "stepsize = " << stepsize << '\n'
            << "frequency = " << frequency << '\n'
            << "Integrand(*,0): von " << Integrand[0, 0] << " bis "
            << Integrand[n_za - 1, 0] << '\n'
            << "Integrand(0,*): von " << Integrand[0, 0] << " bis "
            << Integrand[0, n_aa - 1] << '\n'
            << "za_grid (Theta): von " << za_grid[0] << " bis "
            << za_grid[za_grid.size() - 1] << '\n'
            << "aa_grid (Phi)  : von " << aa_grid[0] << " bis "
            << aa_grid[aa_grid.size() - 1] << '\n';
  std::cout << "---------------init_xy---------------<<<<<\n";
}

//! AngIntegrate_trapezoid_original
/*! 
    The original function from math_funcs.cc
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.

    \param Integrand The Matrix to be integrated
    \param za_grid Input : The zenith angle grid 
    \param aa_grid Input : The azimuth angle grid 
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid_original(MatrixView Integrand,
                                        ConstVectorView za_grid,
                                        ConstVectorView aa_grid) {
  Index n = za_grid.size();
  Index m = aa_grid.size();
  Vector res1(n);
  assert(same_shape<2>({n, m}, Integrand));

  for (Index i = 0; i < n; ++i) {
    res1[i] = 0.0;

    for (Index j = 0; j < m - 1; ++j) {
      res1[i] += 0.5 * DEG2RAD * (Integrand[i, j] + Integrand[i, j + 1]) *
                 (aa_grid[j + 1] - aa_grid[j]) * sin(za_grid[i] * DEG2RAD);
    }
  }
  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i) {
    res +=
        0.5 * DEG2RAD * (res1[i] + res1[i + 1]) * (za_grid[i + 1] - za_grid[i]);
  }

  //std::cout<<res<<"\n";
  return res;
}
//! AngIntegrate_trapezoid_opt
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.

    \param Integrand The Matrix to be integrated
    \param za_grid Input : The zenith angle grid 
    \param aa_grid Input : The azimuth angle grid 
    
    \return The resulting integral
*/

Numeric AngIntegrate_trapezoid_opt(MatrixView Integrand,
                                   ConstVectorView za_grid,
                                   ConstVectorView aa_grid) {
  Index n = za_grid.size();
  Index m = aa_grid.size();
  Vector res1(n);
  assert(same_shape<2>({n, m}, Integrand));

  for (Index i = 0; i < n; ++i) {
    res1[i] = 0.0;

    for (Index j = 0; j < m - 1; ++j) {
      res1[i] += 0.5 * DEG2RAD * (Integrand[i, j] + Integrand[i, j + 1]) *
                 (aa_grid[j + 1] - aa_grid[j]) * sin(za_grid[i] * DEG2RAD);
    }
  }
  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i) {
    res +=
        0.5 * DEG2RAD * (res1[i] + res1[i + 1]) * (za_grid[i + 1] - za_grid[i]);
  }

  //std::cout<<res<<"\n";
  return res;
}

//! AngIntegrate_trapezoid_fixedstep
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.
    Here we use a fixed stepsize, e.g. 0.1

    \param Integrand The Matrix to be integrated
    \param za_grid Input : The zenith angle grid 
    \param aa_grid Input : The azimuth angle grid 
    \param stepsize Input : The grid stepsize
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid_fixedstep(MatrixView Integrand,
                                         ConstVectorView za_grid,
                                         ConstVectorView aa_grid,
                                         Numeric stepsize) {
  Index n = za_grid.size();
  Index m = aa_grid.size();
  Vector res1(n);
  assert(same_shape<2>({n, m}, Integrand));

  for (Index i = 0; i < n; ++i) {
    res1[i] = 0.0;

    for (Index j = 0; j < m - 1; ++j) {
      res1[i] += 0.5 * DEG2RAD * (Integrand[i, j] + Integrand[i, j + 1]) *
                 stepsize * sin(za_grid[i] * DEG2RAD);
    }
  }
  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i) {
    res += 0.5 * DEG2RAD * (res1[i] + res1[i + 1]) * stepsize;
  }

  //std::cout<<res<<"\n";
  return res;
}

//! AngIntegrate_trapezoid_fixedstep_opt
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.
    Here we use a fixed stepsize, e.g. 0.1

    \param Integrand The Matrix to be integrated
    \param za_grid Input : The zenith angle grid 
    \param aa_grid Input : The azimuth angle grid 
    \param stepsize Input : The grid stepsize
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid_fixedstep_opt(MatrixView Integrand,
                                             ConstVectorView za_grid,
                                             ConstVectorView aa_grid,
                                             Numeric stepsize) {
  Index n = za_grid.size();
  Index m = aa_grid.size();
  Vector res1(n);
  assert(same_shape<2>({n, m}, Integrand));

  for (Index i = 0; i < n; ++i) {
    res1[i] = 0.0;

    res1[i] += Integrand[i, 0];
    for (Index j = 1; j < m - 1; j++) {
      res1[i] += Integrand[i, j] * 2;
    }
    res1[i] += Integrand[i, m - 1];
    res1[i] *= 0.5 * DEG2RAD * stepsize * sin(za_grid[i] * DEG2RAD);
  }
  Numeric res  = 0.0;
  res         += res1[0];
  for (Index i = 1; i < n - 1; i++) {
    res += res1[i] * 2;
  }
  res += res1[n - 1];
  res *= 0.5 * DEG2RAD * stepsize;

  //std::cout<<res<<"\n";
  return res;
}

//! AngIntegrate_trapezoid_fixedstep_opt2
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.
    Here we use a fixed stepsize, e.g. 0.1
    with some more optimizations to see whether it is faster or not

    \param Integrand The Matrix to be integrated
    \param za_grid Input : The zenith angle grid 
    \param aa_grid Input : The azimuth angle grid 
    \param stepsize Input: The grid stepsize
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid_fixedstep_opt2(MatrixView Integrand,
                                              ConstVectorView za_grid,
                                              ConstVectorView aa_grid,
                                              Numeric stepsize) {
  Index n = za_grid.size();
  Index m = aa_grid.size();
  Vector res1(n);
  assert(same_shape<2>({n, m}, Integrand));

  Numeric temp = 0.0;

  for (Index i = 0; i < n; ++i) {
    temp = Integrand[i, 0];
    for (Index j = 1; j < m - 1; j++) {
      temp += Integrand[i, j] * 2;
    }
    temp    += Integrand[i, m - 1];
    temp    *= 0.5 * DEG2RAD * stepsize * sin(za_grid[i] * DEG2RAD);
    res1[i]  = temp;
  }

  Numeric res = res1[0];
  for (Index i = 1; i < n - 1; i++) {
    res += res1[i] * 2;
  }
  res += res1[n - 1];
  res *= 0.5 * DEG2RAD * stepsize;

  //std::cout<<res<<"\n";
  return res;
}

//! AngIntegration_trapezoid_original
/*!
  The original function from math_funcs.cc
*/
Numeric AngIntegrate_trapezoid_original(ConstVectorView Integrand,
                                        ConstVectorView za_grid) {
  Index n = za_grid.size();
  assert(same_shape<1>({n}, Integrand));

  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i) {
    // in this place 0.5 * 2 * PI is calculated:
    res += PI * DEG2RAD *
           (Integrand[i] * sin(za_grid[i] * DEG2RAD) +
            Integrand[i + 1] * sin(za_grid[i + 1] * DEG2RAD)) *
           (za_grid[i + 1] - za_grid[i]);
  }

  //std::cout<<res<<"\n";
  return res;
}

//! AngIntegration_trapezoid_original
/*!
  The original function from math_funcs.cc
*/
Numeric AngIntegrate_trapezoid_fixedstep(ConstVectorView Integrand,
                                         ConstVectorView za_grid,
                                         Numeric stepsize) {
  Index n = za_grid.size();
  assert(same_shape<1>({n}, Integrand));

  Numeric res = 0.0;
  // std::cout << "Stepsize: " << stepsize << '\n';
  res += (Integrand[0] * sin(za_grid[0] * DEG2RAD));
  for (Index i = 1; i < n - 1; ++i) {
    res += (Integrand[i] * sin(za_grid[i] * DEG2RAD) * 2);
    // std::cout << i << '\n';
  }
  res += ((Integrand[n - 1] * sin(za_grid[n - 1] * DEG2RAD)));
  // std::cout << n-1 << '\n';
  // normally ther would be a 2* here, but it's already in the equations above
  res *= PI * DEG2RAD * stepsize;

  //std::cout<<res<<"\n";
  return res;
}

//! test_xy
/*! 
    Performs an integration over a simple Function given by a Matrix
    uses AngIntegrate_trapezoid_original to integrate the funcion

    \param z_size The size of the zenith dimension, e.g. 1801
    \param a_size The size of the azimuth dimension, e.g. 3601
    \param stepsize The size of the steps, e.g. 0.1
    \param frequency Only for test purposes: frequency of integral computation

    \return The resulting integral
*/
Numeric test_xy(int z_size, int a_size, float stepsize, int frequency) {
  std::cout << ">>>>>-----------test_xy---------------\n";
  Matrix Integrand(z_size, a_size);  // function to be integrated
  Vector za_grid(z_size);            // zenith (Theta) values
  Vector aa_grid(a_size);            // azimuth (Phi) values

  init_xy(stepsize, frequency, Integrand, za_grid, aa_grid);

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_original(Integrand, za_grid, aa_grid);
  });

  double error = result / (4 * PI) - 1;

  std::cout.precision(15);
  std::cout << "stepsize is    : " << stepsize << '\n'
            << "z_size         : " << z_size << '\n'
            << "a_size         : " << a_size << '\n'
            << "1 is           : " << result / (4 * PI) << '\n'
            << "The result is  : " << result << '\n'
            << "The error is   : " << error * 100 << " %\n"
            << "Number of loops: " << frequency << '\n'
            << "elapsed time   : " << t.dt << '\n'
            << "----------------test_xy----------<<<<<\n";

  return result;
}

//! test_xy_opt
/*! 
    Performs an integration over a simple Function given by a Matrix
    uses AngIntegrate_trapezoid_opt to integrate the funcion

    \param z_size The size of the zenith dimension, e.g. 1801
    \param a_size The size of the azimuth dimension, e.g. 3601
    \param stepsize The size of the steps, e.g. 0.1
    \param frequency Only for test purposes: frequency of integral computation

    \return The resulting integral
*/
Numeric test_xy_opt(int z_size, int a_size, float stepsize, int frequency) {
  std::cout << ">>>>>-----------test_xy_opt---------------\n";
  Matrix Integrand(z_size, a_size);  // function to be integrated
  Vector za_grid(z_size);            // zenith (Theta) values
  Vector aa_grid(a_size);            // azimuth (Phi) values

  init_xy(stepsize, frequency, Integrand, za_grid, aa_grid);

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_opt(Integrand, za_grid, aa_grid);
  });

  double error = result / (4 * PI) - 1;

  std::cout.precision(15);
  std::cout << "stepsize is    : " << stepsize << '\n'
            << "z_size         : " << z_size << '\n'
            << "a_size         : " << a_size << '\n'
            << "1 is           : " << result / (4 * PI) << '\n'
            << "The result is  : " << result << '\n'
            << "The error is   : " << error * 100 << " %\n"
            << "Number of loops: " << frequency << '\n'
            << "elapsed time   : " << t.dt << '\n'
            << "----------------test_xy_opt----------<<<<<\n";

  return result;
}

//! test_xy_fixedstep
/*! 
    Performs an integration over a simple Function given by a Matrix
    uses AngIntegrate_trapezoid_fixedstep to integrate the funcion

    \param z_size The size of the zenith dimension, e.g. 1801
    \param a_size The size of the azimuth dimension, e.g. 3601
    \param stepsize The size of the steps, e.g. 0.1
    \param frequency Only for test purposes: frequency of integral computation

    \return The resulting integral
*/
Numeric test_xy_fixedstep(int z_size,
                          int a_size,
                          float stepsize,
                          int frequency) {
  std::cout << ">>>>>-----------test_xy_fixedstep---------------\n";
  Matrix Integrand(z_size, a_size);  // function to be integrated
  Vector za_grid(z_size);            // zenith (Theta) values
  Vector aa_grid(a_size);            // azimuth (Phi) values

  init_xy(stepsize, frequency, Integrand, za_grid, aa_grid);

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_fixedstep(
          Integrand, za_grid, aa_grid, stepsize);
  });

  double error = result / (4 * PI) - 1;

  std::cout.precision(15);
  std::cout << "stepsize is    : " << stepsize << '\n'
            << "z_size         : " << z_size << '\n'
            << "a_size         : " << a_size << '\n'
            << "1 is           : " << result / (4 * PI) << '\n'
            << "The result is  : " << result << '\n'
            << "The error is   : " << error * 100 << " %\n"
            << "Number of loops: " << frequency << '\n'
            << "elapsed time   : " << t.dt << '\n'
            << "----------------test_xy_fixedstep----------<<<<<\n";

  return result;
}

//! test_xy_fixedstep_opt
/*! 
    Performs an integration over a simple Function given by a Matrix
    uses AngIntegrate_trapezoid_fixedstep to integrate the funcion

    \param z_size The size of the zenith dimension, e.g. 1801
    \param a_size The size of the azimuth dimension, e.g. 3601
    \param stepsize The size of the steps, e.g. 0.1
    \param frequency Only for test purposes: frequency of integral computation

    \return The resulting integral
*/
Numeric test_xy_fixedstep_opt(int z_size,
                              int a_size,
                              float stepsize,
                              int frequency) {
  std::cout << ">>>>>-----------test_xy_fixedstep_opt---------------\n";
  Matrix Integrand(z_size, a_size);  // function to be integrated
  Vector za_grid(z_size);            // zenith (Theta) values
  Vector aa_grid(a_size);            // azimuth (Phi) values

  init_xy(stepsize, frequency, Integrand, za_grid, aa_grid);

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_fixedstep_opt(
          Integrand, za_grid, aa_grid, stepsize);
  });

  double error = result / (4 * PI) - 1;

  std::cout.precision(15);
  std::cout << "stepsize is    : " << stepsize << '\n'
            << "z_size         : " << z_size << '\n'
            << "a_size         : " << a_size << '\n'
            << "1 is           : " << result / (4 * PI) << '\n'
            << "The result is  : " << result << '\n'
            << "The error is   : " << error * 100 << " %\n"
            << "Number of loops: " << frequency << '\n'
            << "elapsed time   : " << t.dt << '\n'
            << "----------------test_xy_fixedstep_opt----------<<<<<\n";

  return result;
}

//! test_xy_fixedstep_opt2
/*! 
    Performs an integration over a simple Function given by a Matrix
    uses AngIntegrate_trapezoid_fixedstep2 to integrate the funcion

    \param z_size The size of the zenith dimension, e.g. 1801
    \param a_size The size of the azimuth dimension, e.g. 3601
    \param stepsize The size of the steps, e.g. 0.1
    \param frequency Only for test purposes: frequency of integral computation

    \return The resulting integral
*/
Numeric test_xy_fixedstep_opt2(int z_size,
                               int a_size,
                               float stepsize,
                               int frequency) {
  std::cout << ">>>>>-----------test_xy_fixedstep_opt2---------------\n";
  Matrix Integrand(z_size, a_size);  // function to be integrated
  Vector za_grid(z_size);            // zenith (Theta) values
  Vector aa_grid(a_size);            // azimuth (Phi) values

  init_xy(stepsize, frequency, Integrand, za_grid, aa_grid);

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_fixedstep_opt2(
          Integrand, za_grid, aa_grid, stepsize);
  });

  double error = result / (4 * PI) - 1;

  std::cout.precision(15);
  std::cout << "stepsize is    : " << stepsize << '\n'
            << "z_size         : " << z_size << '\n'
            << "a_size         : " << a_size << '\n'
            << "1 is           : " << result / (4 * PI) << '\n'
            << "The result is  : " << result << '\n'
            << "The error is   : " << error * 100 << " %\n"
            << "Number of loops: " << frequency << '\n'
            << "elapsed time   : " << t.dt << '\n'
            << "----------------test_xy_fixedstep_opt2----------<<<<<\n";

  return result;
}

//! test_AngIntegrate_trapezoid_opti
/*! 
    Performs an integration over a simple Function given by a Matrix
    uses AngIntegrate_trapezoid_opti to integrate the funcion

    It uses the original arts function

    \param z_size The size of the zenith dimension, e.g. 1801
    \param a_size The size of the azimuth dimension, e.g. 3601
    \param stepsize The size of the steps, e.g. 0.1
    \param frequency Only for test purposes: frequency of integral computation

    \return The resulting integral
*/
Numeric test_AngIntegrate_trapezoid_opti(int z_size,
                                         int a_size,
                                         float stepsize,
                                         int frequency) {
  std::cout
      << ">>>>>-----------test_AngIntegrate_trapezoid_opti---------------\n";
  Matrix Integrand(z_size, a_size);  // function to be integrated
  Vector za_grid(z_size);            // zenith (Theta) values
  Vector aa_grid(a_size);            // azimuth (Phi) values

  Vector grid_stepsize(2);

  init_xy(stepsize, frequency, Integrand, za_grid, aa_grid);

  grid_stepsize[0] = za_grid[1] - za_grid[0];
  grid_stepsize[1] = aa_grid[1] - aa_grid[0];
  //grid_stepsize[0] = -1;
  //grid_stepsize[1] = -1;

  //  std::cout << za_grid << '\n';
  //  std::cout << grid_stepsize[0] << '\n';
  //  std::cout << grid_stepsize[1] << '\n';

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_opti(
          Integrand, za_grid, aa_grid, grid_stepsize);
  });

  double error = result / (4 * PI) - 1;

  std::cout.precision(15);
  std::cout
      << "stepsize is    : " << stepsize << '\n'
      << "z_size         : " << z_size << '\n'
      << "a_size         : " << a_size << '\n'
      << "1 is           : " << result / (4 * PI) << '\n'
      << "The result is  : " << result << '\n'
      << "The error is   : " << error * 100 << " %\n"
      << "Number of loops: " << frequency << '\n'
      << "elapsed time   : " << t.dt << '\n'
      << "----------------test_AngIntegrate_trapezoid_opti----------<<<<<\n";

  return result;
}

//! test_x
/*! 
    Performs an integration over a simple Function y=x given by a Vector
    uses AngIntegrate_trapezoid_original to integrate the funcion

    \param vsize The size of the Vetor who is integrated, e.g. 1801
    \param stepsize The size of the steps, e.g. 0.1
    \param frequency Only for test purposes: how often is the integral calculated

    \return The resulting integral
*/
Numeric test_x(int vsize, float stepsize, int frequency) {
  std::cout << ">>>>>-----------test_x---------------\n";
  Vector Integrand(vsize);  // function to be integrated
  Vector Theta(vsize);      // Theta values

  init_x(vsize, stepsize, frequency, Integrand, Theta);

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_original(Integrand, Theta);
  });

  // norming the result to Pi
  result = result / (2 * 180);

  double error = result / PI - 1;

  std::cout.precision(15);
  std::cout << "stepsize is    : " << stepsize << '\n'
            << "number of steps: " << vsize << '\n'
            << "1 is           : " << PI / PI << '\n'
            << "The result is  : " << result / PI << '\n'
            << "The error is   : " << error * 100 << "%\n"
            << "Number of loops: " << frequency << '\n'
            << "elapsed time   : " << t.dt << '\n'
            << "---------------test_x-----------<<<<<\n";

  return result;
}
//! test_x_fixedstep
/*! 
    Performs an integration over a simple Function y=x given by a Vector
    uses the other Integration function with constant Theta-spacing.

    \param vsize The size of the Vetor who is integrated, e.g. 1801
    \param frequency Only for test purposes: how often is the integral calculated

    \return The resulting integral
*/
Numeric test_x_fixedstep(int vsize, int frequency) {
  std::cout << ">>>>>-----------test_x_fixedstep---------------\n";
  std::cout.precision(12);
  Vector Integrand(vsize);  // function to be integrated
  Vector Theta(vsize);      // Theta values

  double stepsize;
  stepsize =
      180.0 /
      (vsize - 1);  // attention this only works with eaqually spaced intervals
  std::cout << "Neue berechnete Stepsize: " << stepsize << '\n';

  for (Size i = 0; i < Integrand.size(); i++) Integrand[i] = static_cast<Numeric>(i) * stepsize;

  //Theta is between 0 and 180
  for (Size i = 0; i < Theta.size(); i++) Theta[i] = static_cast<Numeric>(i) * stepsize;

  std::cout << "Integrand: von " << Integrand[0] << " bis "
            << Integrand[Integrand.size() - 1] << '\n'
            << "Theta: von " << Theta[0] << " bis " << Theta[Theta.size() - 1]
            << '\n';

  Numeric result = 0;

  Timing t("");
  t([&]() {
    for (int i = 0; i < frequency; i++)
      result = AngIntegrate_trapezoid_fixedstep(Integrand, Theta, stepsize);
  });

  // norming the result to Pi
  result = result / (2 * 180);

  double error = result / PI - 1;

  std::cout.precision(15);
  std::cout << "stepsize is    : " << stepsize << '\n'
            << "number of steps: " << vsize << '\n'
            << "1 is          : " << PI / PI << '\n'
            << "The result is  : " << result / PI << '\n'
            << "The error is   : " << error * 100 << "%\n"
            << "Number of loops: " << frequency << '\n'
            << "elapsed time   : " << t.dt << '\n'
            << "---------------test_x_fixedstep----------<<<<<\n";

  return result;
}
