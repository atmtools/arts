#include <constants.h>
#include <propagationmatrix.h>

/**
 * @brief Contains the Rosenkranz 2021 absorption models
 * 
 * Currently only water vapour is included
 * 
 */

namespace Absorption::PredefinedModel::PWR2021{

constexpr Index num_h2o_lines = 16;

void compute_h2o(PropagationMatrix& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric& p_pa,
                 const Numeric& t,
                 const Numeric& h2o_vmr) noexcept {
    
    using Constant::pow2;
    using Conversion::pa2hpa, Conversion::hpa2bar;

    if (h2o_vmr <= 0){
        return
    }

    // Line parameters taken from h2o_sdlist.asc
    constexpr std::array<Numeric, num_h2o_lines> frequency_ghz {
        22.23508 , 183.310087, 321.22563 , 325.152888, 380.197353,
       439.150807, 443.018343, 448.001085, 470.888999, 474.689092,
       488.490108, 556.935985, 620.700807, 658.006072, 752.033113,
       916.171582
    };
    constexpr std::array<Numeric, num_h2o_lines> strength_296 {
       1.335e-14, 2.319e-12, 7.657e-14, 2.721e-12, 2.477e-11, 2.137e-12,
       4.440e-13, 2.588e-11, 8.196e-13, 3.268e-12, 6.628e-13, 1.570e-09,
       1.700e-11, 9.033e-13, 1.035e-09, 4.275e-11
    };
    constexpr std::array<Numeric, num_h2o_lines> B{
       2.172, 0.677, 6.262, 1.561, 1.062, 3.643, 5.116, 1.424, 3.645,
       2.411, 2.89 , 0.161, 2.423, 7.921, 0.402, 1.461
    };
    std::array<Numeric, num_h2o_lines> w0_air{
       2.74 , 3.033, 2.426, 2.847, 2.868, 2.055, 1.819, 2.612, 2.169,
       2.366, 2.616, 3.115, 2.468, 3.154, 3.114, 2.695
    };
    constexpr std::array<Numeric, num_h2o_lines> xw_air{
       0.76, 0.62, 0.73, 0.64, 0.54, 0.69, 0.7 , 0.7 , 0.73, 0.71, 0.75,
       0.75, 0.79, 0.73, 0.77, 0.79
    };
    std::array<Numeric, num_h2o_lines> w0_self{
       13.63, 15.01, 10.65, 13.95, 14.4 ,  9.06,  7.96, 13.01,  9.7 ,
       11.24, 13.58, 14.24, 11.94, 13.84, 13.58, 13.55
    };
    constexpr std::array<Numeric, num_h2o_lines> xw_self{
       1.2 , 0.82, 0.54, 0.74, 0.89, 0.52, 0.5 , 0.67, 0.65, 0.64, 0.72,
       1.  , 0.75, 1.  , 0.84, 0.48
    };
    std::array<Numeric, num_h2o_lines> d_air{
       1.2 , 0.82, 0.54, 0.74, 0.89, 0.52, 0.5 , 0.67, 0.65, 0.64, 0.72,
       1.  , 0.75, 1.  , 0.84, 0.48
    };
    std::array<Numeric, num_h2o_lines> xd_air{
       2.6, 1.8, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. 
    };
    std::array<Numeric, num_h2o_lines> d_self{
        0.814,  0.136,  0.278,  1.325,  0.24 ,  0.165, -0.229, -0.615,
       -0.465, -0.72 , -0.36 , -1.693,  0.687, -1.496, -0.878,  0.521
    };
    std::array<Numeric, num_h2o_lines> xd_self{
       0.  , 0.98, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.92, 0.  , 0.  , 0.47        
    };
    constexpr std::array<Numeric, num_h2o_lines> a_air {
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    };
    constexpr std::array<Numeric, num_h2o_lines> a_self {
        0. , 12.6,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0. ,  0. ,  0. ,  0. ,  0. 
    };
    std::array<Numeric, num_h2o_lines> w2_air {
        0.435 , 0.407,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0.    ,  0.  ,  0. ,  0. ,  0. 
    };
    std::array<Numeric, num_h2o_lines> x2_air {
       0.   , 0.412, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   
    };
    std::array<Numeric, num_h2o_lines> w2_self {
       1.91, 1.46, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  , 0.  
    };
    std::array<Numeric, num_h2o_lines> x2_self {
       0.   , 0.571, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   
    };
    std::array<Numeric, num_h2o_lines> d2_air {
        0.   , -0.016,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0. 
    };
    std::array<Numeric, num_h2o_lines> d2_self {
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
    constexpr auto replace_x [] (auto &a, auto &b){
        return a <= 0 ? b : a;
    };
    std::transform(xd_air.begin(), xd_air.end(), xw_air.begin(), xd_air.begin(), replace_x);
    std::transform(xd_self.begin(), xd_self.end(), xw_self.begin(), xd_self.begin(), replace_x);
    std::transform(x2_air.begin(), x2_air.end(), xw_air.begin(), x2_air.begin(), replace_x);
    std::transform(x2_self.begin(), x2_self.end(), xw_self.begin(), x2_self.begin(), replace_x);
    
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
    const auto width = [pdry_bar, theta_line](auto &w, auto &x) {
        return w > 0 ? w * pdry_bar * std::pow(theta_line, x) : 0.0;
    };
    std::transform(w0_air.begin(), w0_air.end(), xw_air.begin(), w0_air.begin(), width);
    std::transform(w0_self.begin(), w0_self.end(), xw_self.begin(), w0_self.begin(), width);
    std::array<Numeric, num_h2o_lines> w0;
    std::transform(w0_air.begin(), w0_air.end(), w0_self.begin(), w0.begin(), std::plus<Numeric>);

    std::transform(w2_air.begin(), w2_air.end(), x2_air.begin(), w2_air.begin(), width);
    std::transform(w2_self.begin(), w2_self.end(), x2_self.begin(), w2_self.begin(), width);
    std::array<Numeric, num_h2o_lines> w2;
    std::transform(w2_air.begin(), w2_air.end(), w2_self.begin(), w2.begin(), std::plus<Numeric>);

    // Speed-dependent parameter
    const auto delta = [pdry_bar, pvap_bar](auto &a, auto &s){
        return pdry_bar * a + pvap_bar * s;
    };
    std::array<Numeric, num_h2o_lines> d2;
    std::transform(d2_air.begin(), d2_air.end(), d2_self.begin(), d2.begin(), delta);

    // Frequency shifts
    const Numeric log_theta_line = std::log(theta_line);
    const 

    const Index nf = f_grid.nelem();
    for (Index iv = 0; iv < nf; iv++) {
        const Numeric f = hz2ghz(f_grid[iv]);

        // Continuum absorbtion
        const Numeric cont = ((c_f * pdry_hpa * std::pow(theta_cont, xc_f)
                             + c_s * pvap_hpa * std::pow(theta_cont, xc_s))
            * pvap_hpa * pow2(f));


    if (const Numeric a = sum_lines(f, c, ga, g0, y0, dv0, f0); a > 0)
      propmat_clearsky.Kjj()[iv] += conv * oxygen_vmr * pow2(f) * a;
  }

}

} // namespace Absorption::PredefinedModel::PWR2021