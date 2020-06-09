#ifndef LINEMIXING_HITRAN_H
#define LINEMIXING_HITRAN_H

#include "complex.h"
#include "constants.h"
#include "matpackIV.h"
#include "mystring.h"

namespace lm_hitran_2017 {
namespace parameters {
  static constexpr Index nSigmx=5000001;  // Max number of spectral points
  static constexpr Index iFile=3;  // Unit Number for File access
  static constexpr Index nBmx=7000;  // Max Number of Bands
  static constexpr Index nLmx=700;  // Max Number of Lines per Band
  static constexpr Index nIsotp=10;  // Number of CO2 Isotopomers
  static constexpr Index Nlifmax=9;  // Max number of l values
  static constexpr Index Jmax=130;  // Max number of j values
  static constexpr Index Nwmax=100000;  // Max Number of W coupling
  
  static constexpr Numeric Ct=1.4387686e0;  // Constant
  static constexpr Numeric T0=296;  // Constant
  static constexpr Numeric CtGamD=1.1325e-08;  // Constant
  static constexpr Numeric aMolAtm=7.33889e+21;  // Constant
  static constexpr Numeric Pi=3.141592654e0;  // Constant
  
  static constexpr auto aMass = stdarrayify(44.e-3,45.e-3,46.e-3,45.e-3,47.e-3,46.e-3,48.e-3,47.e-3,46.e-3,49.e-3);  // Constant
  
  static constexpr auto Qcoef = stdarrayify(
    stdarrayify(-.13617e+01, .94899e+00,-.69259e-03, .25974e-05),  // O(16)-C(12)-O(16)
    stdarrayify(-.20631e+01, .18873e+01,-.13669e-02, .54032e-05),  // O(16)-C(13)-O(16)
    stdarrayify(-.29175e+01, .20114e+01,-.14786e-02, .55941e-05),  // O(16)-C(12)-O(18)
    stdarrayify(-.16558e+02, .11733e+02,-.85844e-02, .32379e-04),  // O(16)-C(12)-O(17)
    stdarrayify(-.44685e+01, .40330e+01,-.29590e-02, .11770e-04),  // O(16)-C(13)-O(18)
    stdarrayify(-.26263e+02, .23350e+02,-.17032e-01, .67532e-04),  // O(16)-C(13)-O(17)
    stdarrayify(-.14811e+01, .10667e+01,-.78758e-03, .30133e-05),  // O(18)-C(12)-O(18)
    stdarrayify(-.17600e+02, .12445e+02,-.91837e-02, .34915e-04)); // O(17)-C(12)-O(18)
};

struct CommonBlock {
struct Bands {
  Index nBand;
  std::array<Index, parameters::nBmx> Isot;
  std::array<Index, parameters::nBmx> nLines;
  std::array<Index, parameters::nBmx> li;
  std::array<Index, parameters::nBmx> lf;
  std::array<String, parameters::nBmx> BandFile;
} Bands;

struct LineSg {
  Eigen::MatrixXd Sig;
  LineSg() : Sig(parameters::nLmx, parameters::nBmx) {};
} LineSg;

struct DipoRigid {
  Eigen::MatrixXd Dipo0;
  DipoRigid() : Dipo0(parameters::nLmx, parameters::nBmx) {};
} DipoRigid;

struct Energy {
  Eigen::MatrixXd  E;
  Energy() : E(parameters::nLmx, parameters::nBmx) {};
} Energy;

struct GamVT0AIR {
  Eigen::MatrixXd HWVT0AIR;
  GamVT0AIR() : HWVT0AIR(parameters::nLmx, parameters::nBmx) {};
} GamVT0AIR;

struct GamSDVT0AIR {
  Eigen::MatrixXd HWSDVT0AIR;
  Eigen::MatrixXd rHWT0AIR;
  GamSDVT0AIR() : HWSDVT0AIR(parameters::nLmx, parameters::nBmx), rHWT0AIR(parameters::nLmx, parameters::nBmx) {};
} GamSDVT0AIR;

struct DTGAMAIR {
  Eigen::MatrixXd BHWAIR;
  DTGAMAIR() : BHWAIR(parameters::nLmx, parameters::nBmx) {};
} DTGAMAIR;

struct GamVT0CO2 {
  Eigen::MatrixXd HWVT0SELF;
  GamVT0CO2() : HWVT0SELF(parameters::nLmx, parameters::nBmx) {};
} GamVT0CO2;

struct GamSDVT0CO2 {
  Eigen::MatrixXd HWSDVT0SELF;
  Eigen::MatrixXd rHWT0SELF;
  GamSDVT0CO2() : HWSDVT0SELF(parameters::nLmx, parameters::nBmx), rHWT0SELF(parameters::nLmx, parameters::nBmx) {};
} GamSDVT0CO2;

struct DTGAMCO2 {
  Eigen::MatrixXd BHWSELF;
  DTGAMCO2() : BHWSELF(parameters::nLmx, parameters::nBmx) {};
} DTGAMCO2;

struct GamVT0H2O {
  Eigen::MatrixXd HWVT0H2O;
  GamVT0H2O() : HWVT0H2O(parameters::nLmx, parameters::nBmx) {};
} GamVT0H2O;

struct GamSDVT0H2O {
  Eigen::MatrixXd HWSDVT0H2O;
  Eigen::MatrixXd rHWT0H2O;
  GamSDVT0H2O() : HWSDVT0H2O(parameters::nLmx, parameters::nBmx), rHWT0H2O(parameters::nLmx, parameters::nBmx) {};
} GamSDVT0H2O;

struct DTGAMH2O {
  Eigen::MatrixXd BHWH2O;
  DTGAMH2O() : BHWH2O(parameters::nLmx, parameters::nBmx) {};
} DTGAMH2O;

struct GamT {
  std::array<Numeric, parameters::nLmx> HWT;
  std::array<Numeric, parameters::nLmx> HWSDV2T;
} GamT;

struct SHIFT {
  std::array<Numeric, parameters::nLmx> shft;
} SHIFT;

struct SHIFT0 {
  Eigen::MatrixXd shft0;
  SHIFT0() : shft0(parameters::nLmx, parameters::nBmx) {};
} SHIFT0;

struct PopuT {
  std::array<Numeric, parameters::nLmx> PopuT;
} PopuT;

struct PopTrf {
  Eigen::MatrixXd PopuT0;
  PopTrf() : PopuT0(parameters::nLmx, parameters::nBmx) {};
} PopTrf;

struct DipoTcm {
  Eigen::MatrixXd DipoT;
  DipoTcm() : DipoT(parameters::nLmx, parameters::nBmx) {};
} DipoTcm;

struct Jiln {
  Eigen::MatrixXi Ji;
  Jiln() : Ji(parameters::nLmx, parameters::nBmx) {};
} Jiln;

struct Jfln {
  Eigen::MatrixXi Jf;
  Jfln() : Jf(parameters::nLmx, parameters::nBmx) {};
} Jfln;

struct FicLSR {
  std::array<Numeric, parameters::nLmx> SSR;
} FicLSR;

struct FicLSI {
  std::array<Numeric, parameters::nLmx> SSI;
} FicLSI;

struct FicLPR {
  std::array<Numeric, parameters::nLmx> AlphR;
} FicLPR;

struct FicLPI {
  std::array<Numeric, parameters::nLmx> AlphI;
} FicLPI;

struct Zss {
  std::array<Complex, parameters::nLmx> ZS;
} Zss;

struct Zaa {
  std::array<Complex, parameters::nLmx> ZA;
} Zaa;

struct Wmatrix {
  Eigen::MatrixXd W;
  Wmatrix() : W(parameters::nLmx, parameters::nLmx) {};
} Wmatrix;

struct Wfittedp {
  Tensor4 W0pp;
  Tensor4 W0pq;
  Tensor4 W0pr;
  Wfittedp() :
  W0pp(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  W0pq(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  W0pr(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax) {}
} Wfittedp;

struct Wfittedq {
  Tensor4 W0qp;
  Tensor4 W0qq;
  Tensor4 W0qr;
  Wfittedq() :
  W0qp(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  W0qq(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  W0qr(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax) {}
} Wfittedq;

struct Wfittedr {
  Tensor4 W0rp;
  Tensor4 W0rq;
  Tensor4 W0rr;
  Wfittedr() :
  W0rp(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  W0rq(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  W0rr(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax) {}
} Wfittedr;

struct Bfittedp {
  Tensor4 B0pp;
  Tensor4 B0pq;
  Tensor4 B0pr;
  Bfittedp() :
  B0pp(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  B0pq(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  B0pr(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax) {}
} Bfittedp;

struct Bfittedq {
  Tensor4 B0qp;
  Tensor4 B0qq;
  Tensor4 B0qr;
  Bfittedq() :
  B0qp(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  B0qq(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  B0qr(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax) {}
} Bfittedq;

struct Bfittedr {
  Tensor4 B0rp;
  Tensor4 B0rq;
  Tensor4 B0rr;
  Bfittedr() :
  B0rp(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  B0rq(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax),
  B0rr(parameters::Nlifmax, parameters::Nlifmax, parameters::Jmax, parameters::Jmax) {}
} Bfittedr;

struct DiagnR {
  Eigen::MatrixXd OpR;
  DiagnR() : OpR(parameters::nLmx, parameters::nLmx) {};
} DiagnR;

struct DiagnI {
  Eigen::MatrixXd OpI;
  DiagnI() : OpI(parameters::nLmx, parameters::nLmx) {};
} DiagnI;

struct YLT {
  std::array<Numeric, parameters::nLmx> YT;
} YLT;
};  // CommonBlock

void compute(Vector& absorption,
             const Numeric p,
             const Numeric t,
             const Numeric xco2,
             const Numeric xh2o,
             const Numeric sigmin,
             const Numeric sigmax,
             const Numeric stotmax,
             const Numeric dsig);

};  // lm_hitran_2017

#endif  // LINEMIXING_HITRAN_H
