#include <iostream>
#include <numbers>

#include "scattering/mie.h"
#include "scattering/maths.h"

using std::numbers::pi_v;

/** Tests Mie calculations.
 *
 * Test Mie calculations against results obtained using
 * the bhmie.f code.
 */
bool test_mie_calculation() {

    // Reference value calculated using the bhmie.f code from [1]

    double q_ext_ref = 1.9609;
    double q_sca_ref = 1.1723;
    double q_back_ref = 1.3609e-3;

    Eigen::ArrayX<double> s11_ref{19};
    s11_ref << 3.01299e+01,
        2.78231e+01,
        2.18554e+01,
        1.44846e+01,
        7.96539e+00,
        3.52364e+00,
        1.19086e+00,
        3.09333e-01,
        1.40275e-01,
        1.91035e-01,
        2.41253e-01,
        2.36120e-01,
        1.86923e-01,
        1.20391e-01,
        6.04771e-02,
        2.13912e-02,
        4.79433e-03,
        2.29062e-03,
        3.06208e-03;

    Eigen::ArrayX<double> s12_ref{19};
    s12_ref << 0.00000E+00,
        -2.63733E-01,
        -8.19629E-01,
        -1.18734E+00,
        -1.07794E+00,
        -6.03791E-01,
        -1.10250E-01,
        1.37221E-01,
        1.11125E-01,
        -4.63547E-02,
        -1.80532E-01,
        -2.21548E-01,
        -1.83775E-01,
        -1.15622E-01,
        -5.61279E-02,
        -2.01833E-02,
        -4.73674E-03,
        -5.12504E-04,
        0.00000E+00;

    Eigen::ArrayX<double> s33_ref{19};
    s33_ref << 3.01299e+01,
        2.78211e+01,
        2.18286e+01,
        1.43885e+01,
        7.78099e+00,
        3.29251e+00,
        9.86726e-01,
        1.74212e-01,
        5.39128e-02,
        9.96356e-02,
        1.09232e-01,
        7.20451e-02,
        2.94318e-02,
        5.58099e-03,
        -6.17044e-04,
        -2.27331e-05,
        -1.66440e-04,
        -1.98402e-03,
        -3.06208e-03;

    Eigen::ArrayX<double> s34_ref{19};
    s34_ref << 0.00000e+00,
        2.12135e-01,
        7.07786e-01,
        1.16811e+00,
        1.31973e+00,
        1.10040e+00,
        6.57542e-01,
        2.15656e-01,
        -6.64961e-02,
        -1.56264e-01,
        -1.16960e-01,
        -3.84481e-02,
        1.73463e-02,
        3.30813e-02,
        2.25111e-02,
        7.08654e-03,
        -7.21911e-04,
        -1.02369e-03,
        0.00000e+00;

    // Compare with MieSphere results.

    Eigen::ArrayX<double> theta = Eigen::ArrayX<double>::LinSpaced(19, 0.0, pi_v<double>);
    double lambda = 2.0 * pi_v<double>;
    double r = 3.0;
    double k = pi_v<double> * 2.0 / lambda;
    scattering::MieSphere<double> sphere(
        lambda,
        r,
        std::complex<double>{1.33, 0.1},
        theta
        );

    double delta = q_ext_ref - sphere.get_extinction_eff();
    if (abs(delta) > 1e-4) return false;

    delta = q_sca_ref - sphere.get_scattering_eff();
    if (abs(delta) > 1e-4) return false;

    delta = q_back_ref - sphere.get_backscattering_eff();
    if (abs(delta) > 1e-4) return false;


    Eigen::ArrayX<double> s11 = k * k * sphere.get_phase_function();
    delta = (s11 - s11_ref).maxCoeff();
    if (abs(delta) > 1e-4) return false;

    scattering::math::Matrix<double> phase_matrix = sphere.get_scattering_matrix_compact();

    delta = (k * k * phase_matrix(Eigen::all, 1).array() - s12_ref).maxCoeff();
    if (abs(delta) > 1e-4) return false;
    delta = (k * k * phase_matrix(Eigen::all, 3).array() - s33_ref).maxCoeff();
    if (abs(delta) > 1e-4) return false;
    delta = (k * k * phase_matrix(Eigen::all, 4).array() - s34_ref).maxCoeff();
    if (abs(delta) > 1e-4) return false;


    return true;

}


/*
 * Compares Mie calculation with results from ARTS SSDB.
 *
 * This agreement is only up to a few percent. This is likely due to
 * inaccuracies in the sphere diameter.
 */
bool test_mie_liquid() {
    double frequency = 88.8e9;
    double temperature = 230.0;
    Eigen::ArrayX<double> theta = Eigen::ArrayX<double>::LinSpaced(901, 0.0, pi_v<double>);

    auto sphere = scattering::MieSphere<double>::Liquid(
        frequency,
        temperature,
        145925e-8 / 2.0,
        theta
        );

    double c_ext_ref = 3.216e-06;
    double c_abs_ref = 2.223e-06;

    double c_ext = sphere.get_extinction_coeff();
    double c_sca = sphere.get_scattering_coeff();
    double c_abs = c_ext - c_sca;

    double delta = abs(c_ext - c_ext) / c_ext_ref;
    if (!(abs(delta) < 1e-2)) return false;
    delta = abs(c_abs - c_abs_ref) / c_abs_ref;
    if (!(abs(delta) < 1e-2)) return false;


    Eigen::ArrayX<double> pm = sphere.get_phase_function();

    Eigen::ArrayX<double> pm_ref{20};
    pm_ref << 2.70849618e-07, 2.70846645e-07, 2.70837725e-07, 2.70822860e-07,
        2.70802050e-07, 2.70775297e-07, 2.70742603e-07, 2.70703971e-07,
        2.70659403e-07, 2.70608902e-07, 2.70552473e-07, 2.70490119e-07,
        2.70421845e-07, 2.70347657e-07, 2.70267559e-07, 2.70181558e-07,
        2.70089660e-07, 2.69991871e-07, 2.69888200e-07, 2.69778653e-07;

    delta = (pm_ref - pm(Eigen::seq(0, 19, 1))).abs().maxCoeff() / (c_ext_ref - c_abs_ref);
    if (!(abs(delta) < 1e-3)) return false;

    return true;
}


/*
 * Compares Mie calculation with results from ARTS SSDB.
 *
 * This agreement is only up to a few percent. This is likely due to
 * inaccuracies in the sphere diameter.
 */
bool test_mie_ice() {
    double frequency = 88.8e9;
    double temperature = 230.0;
    Eigen::ArrayX<double> theta = Eigen::ArrayX<double>::LinSpaced(901, 0.0, pi_v<double>);

    auto sphere = scattering::MieSphere<double>::Ice(
        frequency,
        temperature,
        5096e-6 / 2.0,
        theta
        );

    double c_ext_ref = 3.91347217e-05;
    double c_abs_ref = 6.71771906e-07;

    double c_ext = sphere.get_extinction_coeff();
    double c_sca = sphere.get_scattering_coeff();
    double c_abs = c_ext - c_sca;

    double delta = abs(c_ext - c_ext_ref) / c_ext_ref;
    if (!(abs(delta) < 1e-2)) return false;
    delta = abs(c_abs - c_abs_ref) / c_abs_ref;
    if (!(abs(delta) < 1e-2)) return false;


    Eigen::ArrayX<double> pm = sphere.get_phase_function();

    Eigen::ArrayX<double> pm_ref{20};
    pm_ref << 3.90861299e-05, 3.90806709e-05, 3.90642974e-05, 3.90370203e-05,
        3.89988575e-05, 3.89498341e-05, 3.88899823e-05, 3.88193415e-05,
        3.87379582e-05, 3.86458857e-05, 3.85431846e-05, 3.84299223e-05,
        3.83061731e-05, 3.81720180e-05, 3.80275451e-05, 3.78728487e-05,
        3.77080302e-05, 3.75331972e-05, 3.73484638e-05, 3.71539504e-05;

    delta = (pm_ref - pm(Eigen::seq(0, 19, 1))).abs().maxCoeff() / (c_ext_ref - c_abs_ref);
    if (!(abs(delta) < 1e-3)) return false;

    return true;
}

int main(int /*nargs*/, char **/*argv*/) {
    bool passed = test_mie_calculation();
    std::cout << "Test MieSphere: ";
    if (!passed) {
        std::cout << "FAILED!" << std::endl;
        return 1;
    }
    std::cout << "PASSED." << std::endl;

    passed = test_mie_liquid();
    std::cout << "Test liquid sphere: ";
    if (!passed) {
        std::cout << "FAILED!" << std::endl;
        return 1;
    }
    std::cout << "PASSED." << std::endl;

    passed = test_mie_ice();
    std::cout << "Test ice sphere: ";
    if (!passed) {
        std::cout << "FAILED!" << std::endl;
        return 1;
    }
    std::cout << "PASSED." << std::endl;
    return 0;
}
