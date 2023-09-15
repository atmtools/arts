#include <arts_conversions.h>
#include <arts_constexpr_math.h>
#include <matpack.h>
#include <rtepack.h>
#include <Faddeeva.hh>
#include <valarray>

/**
 * @brief Contains the Rosenkranz 20xx absorption models
 * 
 * Currently only water vapour, oxygen and nitrogen are included
 * 
 */

namespace Absorption::PredefinedModel::PWR20xx{

void compute_h2o(PropmatVector& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric& p_pa,
                 const Numeric& t,
                 const Numeric& h2o_vmr,
                 const std::valarray<Numeric>& frequency_ghz,
                 const std::valarray<Numeric>& strength_296,
                 const std::valarray<Numeric>& B,
                 const std::valarray<Numeric>& w0_air,
                 const std::valarray<Numeric>& xw_air,
                 const std::valarray<Numeric>& w0_self,
                 const std::valarray<Numeric>& xw_self,
                 const std::valarray<Numeric>& d_air,
                 std::valarray<Numeric>& xd_air,
                 const std::valarray<Numeric>& d_self,
                 std::valarray<Numeric>& xd_self,
                 const std::valarray<Numeric>& a_air,
                 const std::valarray<Numeric>& a_self,
                 const std::valarray<Numeric>& w2_air,
                 std::valarray<Numeric>& x2_air,
                 const std::valarray<Numeric>& w2_self,
                 std::valarray<Numeric>& x2_self,
                 const std::valarray<Numeric>& d2_air,
                 const std::valarray<Numeric>& d2_self,
                 const Numeric &tref_lines,
                 const Numeric &tref_cont,
                 const Numeric &c_f,
                 const Numeric &xc_f,
                 const Numeric &c_s,
                 const Numeric &xc_s) noexcept {
    // Water vapur absorption routine - common to all versions
    // Based on abh2o.f
    using Math::pow2;
    using Constant::inv_pi, Constant::boltzmann_constant;
    using Conversion::pa2hpa, Conversion::hpa2bar, Conversion::hz2ghz;

    if (h2o_vmr <= 0){
        return;
    }

    // Replacement of unknown temperature exponents
    // done on line reading in abh2o.f
    for (std::size_t iline = 0; iline < xd_air.size(); iline++){
        if (xd_air[iline] <= 0) {
            xd_air[iline] = xw_air[iline];
        }
        if (xd_self[iline] <= 0) {
            xd_self[iline] = xw_self[iline];
        }
        if (x2_air[iline] <= 0) {
            x2_air[iline] = xw_air[iline];
        }
        if (x2_self[iline] <= 0) {
            x2_self[iline] = xw_self[iline];
        }
    }
    
    // Unit conversions. In abh2o.f, line parameters are converted from
    // GHz/bar to GHz/mbar, but here we will keep pressure in bar. However,
    // the continuum coefficients use hPa units
    const Numeric p_hpa = pa2hpa(p_pa);
    const Numeric pvap_hpa = h2o_vmr * p_hpa;
    const Numeric pdry_hpa = p_hpa - pvap_hpa;
    const Numeric pvap_bar = hpa2bar(pvap_hpa);
    const Numeric pdry_bar = hpa2bar(pdry_hpa);

    // Relative inverse temperatures for continua/lines
    const Numeric theta_cont = tref_cont / t;
    const Numeric theta_line = tref_lines / t;

    // Line width calculation
    const std::valarray<Numeric> w0 = w0_air * pdry_bar * std::pow(theta_line, xw_air) +
                    w0_self * pvap_bar * std::pow(theta_line, xw_self);
    const std::valarray<Numeric> w2 = w2_air * pdry_bar * std::pow(theta_line, x2_air) + 
                    w2_self * pvap_bar * std::pow(theta_line, x2_self);                

    // Speed-dependent parameter
    const std::valarray<Numeric> d2 = d2_air * pdry_bar + d2_self * pvap_bar;

    // Frequency shifts
    const Numeric log_theta_line = std::log(theta_line);
    const std::valarray<Numeric> shift_f = d_air * pdry_bar * 
                         (1.0 - a_air * log_theta_line) * std::pow(theta_line, xd_air);
    const std::valarray<Numeric> shift_s = d_self * pvap_bar *
                         (1.0 - a_self * log_theta_line) * std::pow(theta_line, xd_self);
    const std::valarray<Numeric> shift = shift_f + shift_s;

    // Line strengths
    const std::valarray<Numeric> strength = strength_296 * std::pow(theta_line,2.5) *
                          std::exp(B * (1.0 - theta_line));                         

    // Definition of local line contribution
    constexpr Numeric line_cutoff = 750.0;
    const std::valarray<Numeric> base = w0 / (pow2(line_cutoff) + pow2(w0));

    const Index nf = f_grid.nelem();
    for (Index iv = 0; iv < nf; ++iv) {
        const Numeric f = hz2ghz(f_grid[iv]);
        // Continuum absorbtion
        constexpr Numeric conv_cont = 1e-3; // Unit conversion km^-1 to m^-1
        const Numeric cont = ((c_f * pdry_hpa * std::pow(theta_cont, xc_f)
                             + c_s * pvap_hpa * std::pow(theta_cont, xc_s))
            * pvap_hpa * pow2(f) * conv_cont);

        // Line absorption
        const std::valarray<Numeric> df_1 = f - frequency_ghz - shift;
        const std::valarray<Numeric> df_2 = f + frequency_ghz + shift;

        Numeric line_sum = 0.0;
        for (std::size_t iline = 0; iline < frequency_ghz.size(); iline++){
            Numeric resonant = 0.0;
            if ((w2[iline] > 0) && 
                (std::abs(df_1[iline]) < (10.0*w0[iline]))){
                    // Speed-dependent resonant shape factor
                    const Complex denom = Complex(w2[iline], -d2[iline]);
                    const Complex xc = Complex(w0[iline]-1.5*w2[iline], df_1[iline]+1.5*d2[iline]) /
                                    denom;                    
                    const Complex xrt = std::sqrt(xc);
                    constexpr Numeric magic_number = 1.77245385090551603;
                    const Complex pxw = magic_number * xrt * Faddeeva::erfcx(xrt);
                    const Complex sd = 2.0 * (1.0 - pxw) / denom;
                    resonant += sd.real() - base[iline];
                }
            else if (std::abs(df_1[iline]) < line_cutoff) {
                // Lorentzian
                resonant += w0[iline] / (pow2(df_1[iline]) + pow2(w0[iline])) - base[iline];
            }
            if (std::abs(df_2[iline]) < line_cutoff) {
                // Lorentzian for negative frequency
                resonant += w0[iline] / (pow2(df_2[iline]) + pow2(w0[iline])) - base[iline];
            }
            line_sum += strength[iline] * resonant * pow2(f / frequency_ghz[iline]);
        }
        constexpr Numeric conv = 1e-13; // Conversion from Hz / GHz cm^2 / m^2 m^-1 to m^-1
        line_sum = conv * inv_pi * line_sum * p_pa * h2o_vmr / (boltzmann_constant * t);
        propmat_clearsky[iv].A() += line_sum + cont;
  }
}

void compute_h2o_2021(PropmatVector& propmat_clearsky,
                         const Vector& f_grid,
                         const Numeric& p_pa,
                         const Numeric& t,
                         const Numeric& h2o_vmr) noexcept {

    // Line parameters taken from h2o_sdlist.asc
    static const std::valarray<Numeric> frequency_ghz = {
        22.23508 , 183.310087, 321.22563 , 325.152888, 380.197353,
       439.150807, 443.018343, 448.001085, 470.888999, 474.689092,
       488.490108, 556.935985, 620.700807, 658.006072, 752.033113,
       916.171582
    };
    static const std::valarray<Numeric> strength_296 = {
       1.335e-14, 2.319e-12, 7.657e-14, 2.721e-12, 2.477e-11, 2.137e-12,
       4.440e-13, 2.588e-11, 8.196e-13, 3.268e-12, 6.628e-13, 1.570e-09,
       1.700e-11, 9.033e-13, 1.035e-09, 4.275e-11
    };
    static const std::valarray<Numeric> B = {
       2.172, 0.677, 6.262, 1.561, 1.062, 3.643, 5.116, 1.424, 3.645,
       2.411, 2.89 , 0.161, 2.423, 7.921, 0.402, 1.461
    };
    static const std::valarray<Numeric> w0_air = {
       2.74 , 3.033, 2.426, 2.847, 2.868, 2.055, 1.819, 2.612, 2.169,
       2.366, 2.616, 3.115, 2.468, 3.154, 3.114, 2.695
    };
    static const std::valarray<Numeric> xw_air = {
       0.76, 0.62, 0.73, 0.64, 0.54, 0.69, 0.7 , 0.7 , 0.73, 0.71, 0.75,
       0.75, 0.79, 0.73, 0.77, 0.79
    };
    static const std::valarray<Numeric> w0_self = {
       13.63, 15.01, 10.65, 13.95, 14.4 ,  9.06,  7.96, 13.01,  9.7 ,
       11.24, 13.58, 14.24, 11.94, 13.84, 13.58, 13.55
    };
    static const std::valarray<Numeric> xw_self = {
       1.2 , 0.82, 0.54, 0.74, 0.89, 0.52, 0.5 , 0.67, 0.65, 0.64, 0.72,
       1.  , 0.75, 1.  , 0.84, 0.48
    };
    static const std::valarray<Numeric> d_air = {
       -0.033 , -0.074, -0.143, -0.013, -0.074, 0.051, 0.140 , -0.116, 0.061, 
       -0.027,  -0.065,  0.187, 0.0, 0.176  , 0.162, 0.0
    };
    std::valarray<Numeric> xd_air = {
       2.6, 1.8, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. 
    };
    static const std::valarray<Numeric> d_self = {
        0.814,  0.136,  0.278,  1.325,  0.24 ,  0.165, -0.229, -0.615,
       -0.465, -0.72 , -0.36 , -1.693,  0.687, -1.496, -0.878,  0.521
    };
    std::valarray<Numeric> xd_self = {
       0.  , 0.98, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.92, 0.  , 0.  , 0.47        
    };
    static const std::valarray<Numeric> a_air = {
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    };
    static const std::valarray<Numeric> a_self = {
        0. , 12.6,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0. ,  0. ,  0. ,  0. ,  0. 
    };
    static const std::valarray<Numeric> w2_air = {
        0.435 , 0.407,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0.    ,  0.  ,  0. ,  0. ,  0. 
    };
    std::valarray<Numeric> x2_air = {
       0.   , 0.412, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   
    };
    static const std::valarray<Numeric> w2_self = {
       1.91, 1.46, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  , 0.  
    };
    std::valarray<Numeric> x2_self = {
       0.   , 0.571, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   
    };
    static const std::valarray<Numeric> d2_air = {
        0.   , -0.016,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0. 
    };
    static const std::valarray<Numeric> d2_self = {
       0.  , 0.16, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  , 0.  
    };

    constexpr Numeric tref_lines = 296.0;

    // Continuum parameters taken from h2o_sdlist.asc
    constexpr Numeric tref_cont = 300.0;
    constexpr Numeric c_f = 5.919e-10;
    constexpr Numeric xc_f = 3.0;
    constexpr Numeric c_s = 1.416e-8;
    constexpr Numeric xc_s = 7.5;

    compute_h2o(propmat_clearsky,
                f_grid,
                p_pa,
                t,
                h2o_vmr,
                frequency_ghz,
                strength_296,
                B,
                w0_air,
                xw_air,
                w0_self,
                xw_self,
                d_air,
                xd_air,
                d_self,
                xd_self,
                a_air,
                a_self,
                w2_air,
                x2_air,
                w2_self,
                x2_self,
                d2_air,
                d2_self,
                tref_lines,
                tref_cont,
                c_f,
                xc_f,
                c_s,
                xc_s);
}

void compute_h2o_2022(PropmatVector& propmat_clearsky,
                      const Vector& f_grid,
                      const Numeric& p_pa,
                      const Numeric& t,
                      const Numeric& h2o_vmr) noexcept {

    // Line parameters taken from h2o_sdlist.asc
    static const std::valarray<Numeric> frequency_ghz = {
        22.23508 , 183.310087, 321.22563 , 325.152888, 380.197353,
       439.150807, 443.018343, 448.001085, 470.888999, 474.689092,
       488.490108, 556.935985, 620.700807, 658.006072, 752.033113,
       859.965608, 916.171582, 970.315045, 987.926803, 1097.36487
    };
    static const std::valarray<Numeric> strength_296 = {
       1.334e-14, 2.319e-12, 7.654e-14, 2.720e-12, 2.476e-11, 2.136e-12,
       4.440e-13, 2.587e-11, 8.193e-13, 3.268e-12, 6.628e-13, 1.570e-09,
       1.700e-11, 9.027e-13, 1.035e-09, 5.705e-13, 4.272e-11, 4.806e-11,
       7.528e-10, 4.890e-9
    };
    static const std::valarray<Numeric> B = {
       2.172, 0.677, 6.262, 1.561, 1.062, 3.643, 5.116, 1.424, 3.645,
       2.411, 2.89 , 0.161, 2.423, 7.921, 0.402, 8.163, 1.461, 1.944,
       0.261, 0.754
    };
    static const std::valarray<Numeric> w0_air = {
       2.74 , 3.033, 2.426, 2.847, 2.868, 2.055, 1.819, 2.612, 2.169,
       2.366, 2.616, 3.115, 2.468, 3.154, 3.114, 3.121, 2.695, 2.574,
       2.976, 3.095
    };
    static const std::valarray<Numeric> xw_air = {
       0.76, 0.62, 0.73, 0.64, 0.54, 0.69, 0.7 , 0.7 , 0.73, 0.71, 0.75,
       0.75, 0.79, 0.73, 0.77, 0.76, 0.79, 0.70, 0.75, 0.75
    };
    static const std::valarray<Numeric> w0_self = {
       13.63, 15.01, 10.65, 13.95, 14.4 ,  9.06,  7.96, 13.01,  9.7 ,
       11.24, 13.58, 14.24, 11.94, 13.84, 13.58, 14.08, 13.55, 25.95,
       14.35, 15.27 
    };
    static const std::valarray<Numeric> xw_self = {
       1.2 , 0.82, 0.54, 0.74, 0.89, 0.52, 0.5 , 0.67, 0.65, 0.64, 0.72,
       1.  , 0.75, 1.  , 0.84, 0.76, 0.48, 0.7 , 0.75, 0.75
    };
    static const std::valarray<Numeric> d_air = {
       -0.033, -0.074, -0.143, -0.013, -0.074, 0.051, 0.140, -0.116,  0.061, 
       -0.027, -0.065,  0.187,  0.001,  0.176, 0.162, 0.005, -0.001, -0.003, 
       -0.002,  0.002
    };
    std::valarray<Numeric> xd_air = {
       2.6, 1.8, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0., 0., 0., 0., 0. 
    };
    static const std::valarray<Numeric> d_self = {
        0.814,  0.136,  0.278,  1.325,  0.24 ,  0.165, -0.229, -0.615,
       -0.465, -0.72 , -0.36 , -1.693,  0.687, -1.496, -0.878,  0.,
       0.521, 0., 0., 0.
    };
    std::valarray<Numeric> xd_self = {
       0.  , 0.98, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.92, 0.  , 0.  , 0., 0.47, 0., 0., 0.        
    };
    static const std::valarray<Numeric> a_air = {
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    };
    static const std::valarray<Numeric> a_self = {
        0. , 12.6,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.
    };
    static const std::valarray<Numeric> w2_air = {
        0.435 , 0.407,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0.    ,  0.  ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0.
    };
    std::valarray<Numeric> x2_air = {
       0.   , 0.412, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.
    };
    static const std::valarray<Numeric> w2_self = {
       1.91, 1.46, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.
    };
    std::valarray<Numeric> x2_self = {
       0.   , 0.571, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.  
    };
    static const std::valarray<Numeric> d2_air = {
        0.   , -0.016,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        0.   ,  0.   ,  0.   ,  0.   
    };
    static const std::valarray<Numeric> d2_self = {
       0.  , 0.16, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.
    };

    constexpr Numeric tref_lines = 296.0;

    // Continuum parameters taken from h2o_sdlist.asc
    constexpr Numeric tref_cont = 300.0;
    constexpr Numeric c_f = 5.9197e-10;
    constexpr Numeric xc_f = 3.0;
    constexpr Numeric c_s = 1.4162e-8;
    constexpr Numeric xc_s = 7.5;

    compute_h2o(propmat_clearsky,
                f_grid,
                p_pa,
                t,
                h2o_vmr,
                frequency_ghz,
                strength_296,
                B,
                w0_air,
                xw_air,
                w0_self,
                xw_self,
                d_air,
                xd_air,
                d_self,
                xd_self,
                a_air,
                a_self,
                w2_air,
                x2_air,
                w2_self,
                x2_self,
                d2_air,
                d2_self,
                tref_lines,
                tref_cont,
                c_f,
                xc_f,
                c_s,
                xc_s);
}

void compute_o2(PropmatVector& propmat_clearsky,
                const Vector& f_grid,
                const Numeric& p_pa,
                const Numeric& t,
                const Numeric& o2_vmr,
                const Numeric& h2o_vmr,
                const std::valarray<Numeric>& frequency_ghz,
                const std::valarray<Numeric>& strength_300,
                const std::valarray<Numeric>& be,
                const std::valarray<Numeric>& width_300,
                const std::valarray<Numeric>& y0,
                const std::valarray<Numeric>& y1,
                const std::valarray<Numeric>& g0,
                const std::valarray<Numeric>& g1,
                const std::valarray<Numeric>& dnu0,
                const std::valarray<Numeric>& dnu1
                ) noexcept {
    // Oxygen absorption routine
    // Based on o2abs_19.f
    using Math::pow2, Math::pow3;
    using Constant::inv_pi, Constant::boltzmann_constant;
    using Conversion::hz2ghz, Conversion::pa2bar;

    // Continuum line width
    constexpr Numeric cont_width_300 = 0.56;
    
    constexpr Numeric x = 0.754;

    // Relative inverse temperature
    constexpr Numeric t_ref = 300.0;
    const Numeric theta = t_ref / t;
    const Numeric theta_minus_1 = theta - 1.0;
    
    const Numeric b = std::pow(theta, x);

    const Numeric pvap_pa = h2o_vmr * p_pa;
    const Numeric pdry_pa = p_pa - pvap_pa;
    const Numeric pvap_bar = pa2bar(pvap_pa);
    const Numeric pdry_bar = pa2bar(pdry_pa);

    // Factor of 1.2 accounts for higher broadening efficiency
    // for water vapour compared to dry air
    const Numeric den = pdry_bar*b + 1.2*pvap_bar*theta;

    const Numeric df_cont = cont_width_300 * den;
    const Numeric pe2 = pow2(den);

    const std::valarray<Numeric> y = den * (y0 + y1 * theta_minus_1);
    const std::valarray<Numeric> delta_nu = pe2 * (dnu0 + dnu1 * theta_minus_1);
    const std::valarray<Numeric> g = 1.0 + pe2 * (g0 + g1 * theta_minus_1);
    const std::valarray<Numeric> width = width_300 * den;
    const std::valarray<Numeric> strength = strength_300 * std::exp(-be * theta_minus_1);

    const Index nf = f_grid.nelem();
    for (Index iv = 0; iv < nf; ++iv) {
        const Numeric f_ghz = hz2ghz(f_grid[iv]);
        const Numeric f2_ghz = pow2(f_ghz);
        const Numeric cont = 1.584e-17*f2_ghz*df_cont / (theta * (f2_ghz + pow2(df_cont)));
        const std::valarray<Numeric> df_1 = f_ghz - frequency_ghz - delta_nu;
        const std::valarray<Numeric> df_2 = f_ghz + frequency_ghz + delta_nu;
        const std::valarray<Numeric> den_1 = pow2(df_1) + pow2(width);
        const std::valarray<Numeric> den_2 = pow2(df_2) + pow2(width);
        const std::valarray<Numeric> sfac_1 = (width * g + df_1*y) / den_1;
        const std::valarray<Numeric> sfac_2 = (width * g - df_2*y) / den_2;
        const std::valarray<Numeric> lines = strength * (sfac_1 + sfac_2) * pow2((f_ghz / frequency_ghz));
        const Numeric sum = lines.sum() + cont;
        
        constexpr Numeric conv = 1e-13; // Conversion from Hz / GHz cm^2 / m^2 m^-1 to m^-1
        // Factor of 1.004 comes from Koshelev 2017 paper according to o2abs_19.f code
        const Numeric absorption = 1.004 * conv * o2_vmr * inv_pi / (boltzmann_constant * t_ref) * sum * pdry_pa * pow3(theta);
        if (absorption > 0){
            propmat_clearsky[iv].A() += absorption;
        }
    }

}

void compute_o2_2021(PropmatVector& propmat_clearsky,
                     const Vector& f_grid,
                     const Numeric& p_pa,
                     const Numeric& t,
                     const Numeric& o2_vmr,
                     const Numeric& h2o_vmr) noexcept {
    // Line parameters
    // Line frequency (GHz)
    static const std::valarray<Numeric> frequency_ghz = {
        118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,
        59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
        56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,
        55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,
        53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,
        52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310,
        50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630,
        487.2493, 566.8956, 715.3929, 731.1866,
        773.8395, 834.1455, 895.0710
    };
    // Line strength at 300K
    static const std::valarray<Numeric> strength_300 = {
        0.2906e-14,0.7957e-15,0.2444e-14,0.2194e-14,
        0.3301e-14,0.3243e-14,0.3664e-14,0.3834e-14,
        0.3588e-14,0.3947e-14,0.3179e-14,0.3661e-14,
        0.2590e-14,0.3111e-14,0.1954e-14,0.2443e-14,
        0.1373e-14,0.1784e-14,0.9013e-15,0.1217e-14,
        0.5545e-15,0.7766e-15,0.3201e-15,0.4651e-15,
        0.1738e-15,0.2619e-15,0.8880e-16,0.1387e-15,
        0.4272e-16,0.6923e-16,0.1939e-16,0.3255e-16,
        0.8301e-17,0.1445e-16,0.3356e-17,0.6049e-17,
        0.1280e-17,0.2394e-17,
        0.3287e-16,0.6463e-15,0.1334e-16,0.7049e-14,
        0.3011e-14,0.1797e-16,0.1826e-14,0.2193e-16,
        0.1153e-13,0.3974e-14,0.2512e-16
    };
    static const std::valarray<Numeric> be = {
    0.010, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.387, 0.621, 0.621,
    0.910, 0.910, 1.255, 1.255, 1.654, 1.654, 2.109, 2.109, 2.618, 2.618, 
    3.182, 3.182, 3.800, 3.800, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 
    6.819, 6.819, 7.709, 7.709, 8.653, 8.653, 9.651, 9.651,
    0.019, 0.048, 0.045, 0.044, 0.049, 0.084, 0.145, 0.136, 0.141, 0.145, 
    0.201
    };
    // Line width at 300K (GHz/bar)
    static const std::valarray<Numeric> width_300 = {
    1.685, 1.703, 1.513, 1.495, 1.433, 1.408,
    1.353, 1.353, 1.303, 1.319, 1.262, 1.265, 
    1.238, 1.217, 1.207, 1.207, 1.137, 1.137,
    1.101, 1.101, 1.037, 1.038, 0.996, 0.996,
    0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 
    0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 
    0.669, 0.669, 1.65,  1.64,  1.64,  1.64, 
    1.60,  1.60,  1.60,  1.60,  1.62,  1.47, 
    1.47        
    };
    // First order mixing coefficients (1/bar)
    static const std::valarray<Numeric> y0 = {
    -0.041, 0.277, -0.373, 0.560, -0.573, 0.618,
    -0.366, 0.278, -0.089, -0.021, 0.0599, -0.152,
    0.216, -0.293, 0.374, -0.436, 0.491, -0.542,
    0.571, -0.613, 0.636, -0.670, 0.690, -0.718,
    0.740, -0.763, 0.788, -0.807, 0.834, -0.849,
    0.876, -0.887, 0.915, -0.922, 0.950, -0.955,
    0.987, -0.988, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    static const std::valarray<Numeric> y1 = {
    0.000, 0.11, -0.009, 0.007, 0.049, -0.1,
    0.260, -0.346, 0.364, -0.422, 0.315, -0.341,
    0.483, -0.503, 0.598, -0.610, 0.630, -0.633,
    0.613, -0.611, 0.570, -0.564, 0.58, -0.57,
    0.61, -0.60, 0.64, -0.62, 0.65, -0.64,
    0.66, -0.64, 0.66, -0.64, 0.66, -0.64,
    0.65, -0.63, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    // Second order mixing coefficients (1/bar^2)
    static const std::valarray<Numeric> g0 = {
    -0.000695, -0.090, -0.103, -0.239, -0.172, -0.171,
    0.028, 0.150, 0.132, 0.170, 0.087, 0.069, 
    0.083, 0.068, 0.007, 0.016, -0.021, -0.066,
    -0.095, -0.116, -0.118, -0.140, -0.173, -0.186,
    -0.217, -0.227, -0.234, -0.242, -0.266, -0.272,
    -0.301, -0.304, -0.334, -0.333, -0.362, -0.358,
    -0.348, -0.344, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    static const std::valarray<Numeric> g1 = {
    0.000, -0.042, 0.004, 0.025, 0.083, 0.167,
    0.178, 0.223, 0.054, 0.003, 0.002, -0.044,
    -0.019, -0.054, -0.177, -0.208, -0.294, -0.334,
    -0.368, -0.386, -0.374, -0.384, -0.387, -0.389,
    -0.423, -0.422, -0.46, -0.46, -0.51, -0.50,
    -0.55, -0.53, -0.58, -0.56, -0.62, -0.59,
    -0.68, -0.65, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    // Pressure-induced frequency shift coefficients (GHz/bar^2)
    static const std::valarray<Numeric> dnu0 = {
    -0.00028, 0.00596, -0.01950, 0.032, -0.0475, 0.0541,
    -0.0232, 0.0155, 0.0007, -0.0086, -0.0026, -0.0013,
    -0.0004, -0.002, 0.005, -0.007, 0.007, -0.008,
    0.006, -0.007, 0.006, -0.006, 0.005, -0.0049,
    0.0040, -0.0041, 0.0036, -0.0037, 0.0033, -0.0034,
    0.0032, -0.0032, 0.0030, -0.0030, 0.0028, -0.0029,
    0.0029, -0.0029, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    static const std::valarray<Numeric> dnu1 = {
    -0.00037, 0.0086, -0.013, 0.019, -0.026, 0.027,
    0.005, -0.014, 0.012, -0.018, -0.015, 0.015,
    0.003, -0.004, 0.012, -0.013, 0.012, -0.012,
    0.009, -0.009, 0.002, -0.002, 0.0005, -0.0005,
    0.002, -0.002, 0.002, -0.002, 0.002, -0.002,
    0.002, -0.002, 0.002, -0.002, 0.001, -0.001,
    0.0004, -0.0004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };

    compute_o2(propmat_clearsky,
               f_grid,
               p_pa,
               t,
               o2_vmr,
               h2o_vmr,
               frequency_ghz,
               strength_300,
               be,
               width_300,
               y0,
               y1,
               g0,
               g1,
               dnu0,
               dnu1
            );
}

void compute_o2_2022(PropmatVector& propmat_clearsky,
                     const Vector& f_grid,
                     const Numeric& p_pa,
                     const Numeric& t,
                     const Numeric& o2_vmr,
                     const Numeric& h2o_vmr) noexcept {
    // Line parameters
    // Line frequency (GHz)
    static const std::valarray<Numeric> frequency_ghz = {
        118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,
        59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
        56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,
        55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,
        53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,
        52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310,
        50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630,
        487.2493, 566.8956, 715.3929, 731.1866,
        773.8395, 834.1455, 895.0710
    };
    // Line strength at 300K
    static const std::valarray<Numeric> strength_300 = {
        0.2906e-14,0.7957e-15,0.2444e-14,0.2194e-14,
        0.3301e-14,0.3243e-14,0.3664e-14,0.3834e-14,
        0.3588e-14,0.3947e-14,0.3179e-14,0.3661e-14,
        0.2590e-14,0.3111e-14,0.1954e-14,0.2443e-14,
        0.1373e-14,0.1784e-14,0.9013e-15,0.1217e-14,
        0.5545e-15,0.7766e-15,0.3201e-15,0.4651e-15,
        0.1738e-15,0.2619e-15,0.8880e-16,0.1387e-15,
        0.4272e-16,0.6923e-16,0.1939e-16,0.3255e-16,
        0.8301e-17,0.1445e-16,0.3356e-17,0.6049e-17,
        0.1280e-17,0.2394e-17,
        0.8299e-17,0.6492e-15,0.1354e-16,0.7082e-14,
        0.3026e-14,0.1823e-16,0.1835e-14,0.2226e-16,
        0.1158e-13,0.3992e-14,0.2550e-16
    };
    static const std::valarray<Numeric> be = {
    0.010, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.387, 0.621, 0.621,
    0.910, 0.910, 1.255, 1.255, 1.654, 1.654, 2.109, 2.109, 2.618, 2.618, 
    3.182, 3.182, 3.800, 3.800, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 
    6.819, 6.819, 7.709, 7.709, 8.653, 8.653, 9.651, 9.651,
    0.019, 0.048, 0.045, 0.044, 0.049, 0.084, 0.145, 0.136, 0.140, 0.145, 
    0.201
    };
    // Line width at 300K (GHz/bar)
    static const std::valarray<Numeric> width_300 = {
    1.685, 1.703, 1.513, 1.495, 1.433, 1.408,
    1.353, 1.353, 1.303, 1.319, 1.262, 1.265, 
    1.238, 1.217, 1.207, 1.207, 1.137, 1.137,
    1.101, 1.101, 1.037, 1.038, 0.996, 0.996,
    0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 
    0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 
    0.669, 0.669, 1.65,  1.64,  1.64,  1.64, 
    1.60,  1.60,  1.60,  1.60,  1.62,  1.47, 
    1.47        
    };
    // First order mixing coefficients (1/bar)
    static const std::valarray<Numeric> y0 = {
    -0.041, 0.277, -0.373, 0.560, -0.573, 0.618,
    -0.366, 0.278, -0.089, -0.021, 0.0599, -0.152,
    0.216, -0.293, 0.374, -0.436, 0.491, -0.542,
    0.571, -0.613, 0.636, -0.670, 0.690, -0.718,
    0.740, -0.763, 0.788, -0.807, 0.834, -0.849,
    0.876, -0.887, 0.915, -0.922, 0.950, -0.955,
    0.987, -0.988, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    static const std::valarray<Numeric> y1 = {
    0.000, 0.11, -0.009, 0.007, 0.049, -0.1,
    0.260, -0.346, 0.364, -0.422, 0.315, -0.341,
    0.483, -0.503, 0.598, -0.610, 0.630, -0.633,
    0.613, -0.611, 0.570, -0.564, 0.58, -0.57,
    0.61, -0.60, 0.64, -0.62, 0.65, -0.64,
    0.66, -0.64, 0.66, -0.64, 0.66, -0.64,
    0.65, -0.63, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    // Second order mixing coefficients (1/bar^2)
    static const std::valarray<Numeric> g0 = {
    -0.000695, -0.090, -0.103, -0.239, -0.172, -0.171,
    0.028, 0.150, 0.132, 0.170, 0.087, 0.069, 
    0.083, 0.068, 0.007, 0.016, -0.021, -0.066,
    -0.095, -0.116, -0.118, -0.140, -0.173, -0.186,
    -0.217, -0.227, -0.234, -0.242, -0.266, -0.272,
    -0.301, -0.304, -0.334, -0.333, -0.362, -0.358,
    -0.348, -0.344, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    static const std::valarray<Numeric> g1 = {
    0.000, -0.042, 0.004, 0.025, 0.083, 0.167,
    0.178, 0.223, 0.054, 0.003, 0.002, -0.044,
    -0.019, -0.054, -0.177, -0.208, -0.294, -0.334,
    -0.368, -0.386, -0.374, -0.384, -0.387, -0.389,
    -0.423, -0.422, -0.46, -0.46, -0.51, -0.50,
    -0.55, -0.53, -0.58, -0.56, -0.62, -0.59,
    -0.68, -0.65, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    // Pressure-induced frequency shift coefficients (GHz/bar^2)
    static const std::valarray<Numeric> dnu0 = {
    -0.00028, 0.00596, -0.01950, 0.032, -0.0475, 0.0541,
    -0.0232, 0.0155, 0.0007, -0.0086, -0.0026, -0.0013,
    -0.0004, -0.002, 0.005, -0.007, 0.007, -0.008,
    0.006, -0.007, 0.006, -0.006, 0.005, -0.0049,
    0.0040, -0.0041, 0.0036, -0.0037, 0.0033, -0.0034,
    0.0032, -0.0032, 0.0030, -0.0030, 0.0028, -0.0029,
    0.0029, -0.0029, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };
    static const std::valarray<Numeric> dnu1 = {
    -0.00037, 0.0086, -0.013, 0.019, -0.026, 0.027,
    0.005, -0.014, 0.012, -0.018, -0.015, 0.015,
    0.003, -0.004, 0.012, -0.013, 0.012, -0.012,
    0.009, -0.009, 0.002, -0.002, 0.0005, -0.0005,
    0.002, -0.002, 0.002, -0.002, 0.002, -0.002,
    0.002, -0.002, 0.002, -0.002, 0.001, -0.001,
    0.0004, -0.0004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    };

    compute_o2(propmat_clearsky,
               f_grid,
               p_pa,
               t,
               o2_vmr,
               h2o_vmr,
               frequency_ghz,
               strength_300,
               be,
               width_300,
               y0,
               y1,
               g0,
               g1,
               dnu0,
               dnu1
            );
}

void compute_n2(PropmatVector& propmat_clearsky,
                const Vector& f_grid,
                const Numeric& p_pa,
                const Numeric& t,
                const Numeric& n2_vmr,
                const Numeric& h2o_vmr) noexcept {
    // Dry air continuum absorption routine based on absn2.f
    // Note that in spite of the name, this is for air and
    // not pure Nitrogen, i.e. it includes O2-N2 and O2-O2 
    // collision induced absorption. As such, it should only be 
    // used with N2 (dry) vmr = 0.781 and  the N2/O2 rato applicable to 
    // Earth. However, the N2-vmr dependence is included as
    // the absorption must have a non-zero Jacobian for
    // the calculation of partial absorption in 
    // propmat_clearsky_fieldCalc. 

    using Math::pow2;
    using Conversion::pa2hpa, Conversion::hz2ghz;

    const Numeric theta = 300.0 / t;
    const Numeric pdry_pa = p_pa * (1.0 - h2o_vmr);
    const Numeric pdry_hpa = pa2hpa(pdry_pa);
    // absn2.f includes the N2 vmr in the continuum coefficient
    // Here, we remove the (assumed) value to permit absorption
    // to vary with N2 VMR - ARTS requires this for 
    // propmat_clearsky_fieldCalc to work properly
    constexpr Numeric assumed_n2_vmr = 0.781;
    constexpr Numeric continuum_coefficient = 9.95e-14; // Np / km / (hPa GHz)^2
    const Numeric cont = (n2_vmr / assumed_n2_vmr) * continuum_coefficient * pow2(pdry_hpa) * std::pow(theta, 3.22);
    const Index nf = f_grid.nelem();
    for (Index iv = 0; iv < nf; ++iv) {
        const Numeric f_ghz = hz2ghz(f_grid[iv]);
        const Numeric frequency_dependence = 0.5 + 0.5/(1.0 + pow2(f_ghz / 450.0));
        const Numeric continuum = cont * frequency_dependence * pow2(f_ghz) / 1000.0;
        propmat_clearsky[iv].A() += continuum;
    }
}

} // namespace Absorption::PredefinedModel::PWR20xx