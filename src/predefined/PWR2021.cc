#include <constants.h>
#include <propagationmatrix.h>
#include <Faddeeva.hh>
#include <valarray>

/**
 * @brief Contains the Rosenkranz 2021 absorption models
 * 
 * Currently only water vapour is included
 * 
 */

namespace Absorption::PredefinedModel::PWR2021{

void compute_h2o(PropagationMatrix& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric& p_pa,
                 const Numeric& t,
                 const Numeric& h2o_vmr) noexcept {
    
    using Constant::pow2, Constant::inv_pi, Constant::boltzmann_constant;
    using Conversion::pa2hpa, Conversion::hpa2bar, Conversion::hz2ghz;

    if (h2o_vmr <= 0){
        return;
    }

    // Line parameters taken from h2o_sdlist.asc
    const std::valarray<Numeric> frequency_ghz = {
        22.23508 , 183.310087, 321.22563 , 325.152888, 380.197353,
       439.150807, 443.018343, 448.001085, 470.888999, 474.689092,
       488.490108, 556.935985, 620.700807, 658.006072, 752.033113,
       916.171582
    };
    const std::valarray<Numeric> strength_296 = {
       1.335e-14, 2.319e-12, 7.657e-14, 2.721e-12, 2.477e-11, 2.137e-12,
       4.440e-13, 2.588e-11, 8.196e-13, 3.268e-12, 6.628e-13, 1.570e-09,
       1.700e-11, 9.033e-13, 1.035e-09, 4.275e-11
    };
    const std::valarray<Numeric> B = {
       2.172, 0.677, 6.262, 1.561, 1.062, 3.643, 5.116, 1.424, 3.645,
       2.411, 2.89 , 0.161, 2.423, 7.921, 0.402, 1.461
    };
    const std::valarray<Numeric> w0_air = {
       2.74 , 3.033, 2.426, 2.847, 2.868, 2.055, 1.819, 2.612, 2.169,
       2.366, 2.616, 3.115, 2.468, 3.154, 3.114, 2.695
    };
    const std::valarray<Numeric> xw_air = {
       0.76, 0.62, 0.73, 0.64, 0.54, 0.69, 0.7 , 0.7 , 0.73, 0.71, 0.75,
       0.75, 0.79, 0.73, 0.77, 0.79
    };
    const std::valarray<Numeric> w0_self = {
       13.63, 15.01, 10.65, 13.95, 14.4 ,  9.06,  7.96, 13.01,  9.7 ,
       11.24, 13.58, 14.24, 11.94, 13.84, 13.58, 13.55
    };
    const std::valarray<Numeric> xw_self = {
       1.2 , 0.82, 0.54, 0.74, 0.89, 0.52, 0.5 , 0.67, 0.65, 0.64, 0.72,
       1.  , 0.75, 1.  , 0.84, 0.48
    };
    const std::valarray<Numeric> d_air = {
       1.2 , 0.82, 0.54, 0.74, 0.89, 0.52, 0.5 , 0.67, 0.65, 0.64, 0.72,
       1.  , 0.75, 1.  , 0.84, 0.48
    };
    std::valarray<Numeric> xd_air = {
       2.6, 1.8, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. 
    };
    const std::valarray<Numeric> d_self = {
        0.814,  0.136,  0.278,  1.325,  0.24 ,  0.165, -0.229, -0.615,
       -0.465, -0.72 , -0.36 , -1.693,  0.687, -1.496, -0.878,  0.521
    };
    std::valarray<Numeric> xd_self = {
       0.  , 0.98, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.92, 0.  , 0.  , 0.47        
    };
    const std::valarray<Numeric> a_air = {
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    };
    const std::valarray<Numeric> a_self = {
        0. , 12.6,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0. ,  0. ,  0. ,  0. ,  0. 
    };
    const std::valarray<Numeric> w2_air = {
        0.435 , 0.407,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0.    ,  0.  ,  0. ,  0. ,  0. 
    };
    std::valarray<Numeric> x2_air = {
       0.   , 0.412, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   
    };
    const std::valarray<Numeric> w2_self = {
       1.91, 1.46, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  , 0.  
    };
    std::valarray<Numeric> x2_self = {
       0.   , 0.571, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   
    };
    const std::valarray<Numeric> d2_air = {
        0.   , -0.016,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0. 
    };
    const std::valarray<Numeric> d2_self = {
       0.  , 0.16, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  , 0.  
    };

    constexpr Numeric tref_lines = 296.0;

    // Continuum parameters taken from h2o_sdlist.asc
    constexpr Numeric tref_cont = 300.0;
    constexpr Numeric c_f = 5.919E-10;
    constexpr Numeric xc_f = 3.0;
    constexpr Numeric c_s = 1.416E-8;
    constexpr Numeric xc_s = 7.5;

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
    const auto w0 = w0_air * pdry_bar * std::pow(theta_line, xw_air) +
                    w0_self * pvap_bar * std::pow(theta_line, xw_self);
    const auto w2 = w2_air * pdry_bar * std::pow(theta_line, x2_air) + 
                    w2_self * pvap_bar * std::pow(theta_line, x2_self);                

    // Speed-dependent parameter
    const auto d2 = d2_air * pdry_bar + d2_self * pvap_bar;

    // Frequency shifts
    const Numeric log_theta_line = std::log(theta_line);
    const auto shift_f = d_air * pdry_bar * 
                         (1.0 - a_air * log_theta_line) * std::pow(theta_line, xd_air);
    const auto shift_s = d_self * pvap_bar *
                         (1.0 - a_self * log_theta_line) * std::pow(theta_line, xd_self);
    const auto shift = shift_f + shift_s;

    // Line strengths
    const auto strength = strength_296 * std::pow(theta_line,2.5) *
                          std::exp(B * (1.0 - theta_line));                         

    // Definition of local line contribution
    constexpr Numeric line_cutoff = 750.0;
    const auto base = w0 / (pow2(line_cutoff) + pow2(w0));

    const Index nf = f_grid.nelem();
    for (Index iv = 0; iv < nf; ++iv) {
        const Numeric f = hz2ghz(f_grid[iv]);
        // Continuum absorbtion
        constexpr Numeric conv_cont = 1e-3; // Unit conversion km^-1 to m^-1
        const Numeric cont = ((c_f * pdry_hpa * std::pow(theta_cont, xc_f)
                             + c_s * pvap_hpa * std::pow(theta_cont, xc_s))
            * pvap_hpa * pow2(f) * conv_cont);

        // Line absorption
        const auto df_1 = f - frequency_ghz - shift;
        const auto df_2 = f + frequency_ghz + shift;

        Numeric line_sum = 0.0;
        for (std::size_t iline = 0; iline < frequency_ghz.size(); iline++){
            Numeric resonant = 0.0;
            if ((w2[iline] > 0) && 
                (std::abs(df_1[iline]) < (10.0*w0[iline]))){
                    // Speed-dependent resonant shape factor
                    const auto denom = std::complex<Numeric> (w2[iline], -d2[iline]);
                    const auto xc = std::complex<Numeric> (w0[iline]-1.5*w2[iline], df_1[iline]+1.5*d2[iline]) /
                                    denom;
                    
                    const auto xrt = std::sqrt(xc);
                    constexpr Numeric magic_number = 1.77245385090551603;
                    const auto pxw = magic_number * xrt * Faddeeva::erfcx(xrt);
                    const auto sd = 2.0 * (1.0 - pxw) / denom;
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
        propmat_clearsky.Kjj()[iv] += line_sum + cont;
  }

}

} // namespace Absorption::PredefinedModel::PWR2021