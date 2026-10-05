// Tests of the RT3 wrapper against references from outside RT3:
//   - Evans' benchmark outputs of the original RT3 program, kept verbatim in
//     3rdparty/polradtran/runmietest (the Mie case of Evans and
//     Stephens, 1991) and runtesta, and copied into the tables below;
//   - closed forms of the transfer equation for gas-only layers;
//   - single scattering by an optically thin Rayleigh layer, built from the
//     vector geometry of the scattered field.
#include <arts_constants.h>
#include <physics_funcs.h>
#include <rt3.h>
#include <rtepack.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
constexpr Numeric pi = Constant::pi;

Index size(const Vector& v) { return static_cast<Index>(v.size()); }

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

//! One line of an rt3.f output file: Z, PHI [deg], MU, then I, Q, U, V.
//! MU = -2 / +2 are the up / down fluxes; MU < 0 is upwelling radiance.
struct table_row {
  Numeric                z, phi_deg, mu;
  std::array<Numeric, 4> iquv;
};

// 3rdparty/polradtran/runmietest, mietest.out.check
// (W m-2 sr-1 um-1 and W m-2 um-1)
constexpr table_row mietest_check[] = {
    {1.000, .0, -2.00000, {.365133E+00, .489285E-01, .000000E+00, .000000E+00}},
    {1.000, .0, 2.00000, {.628319E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, -.98940, {.611250E-01, -.254753E-01, .000000E+00, .000000E+00}},
    {1.000, .0, -.94458, {.759512E-01, -.301313E-01, .000000E+00, .000000E+00}},
    {1.000, .0, -.86563, {.101371E+00, -.337486E-01, .000000E+00, .000000E+00}},
    {1.000, .0, -.75540, {.142650E+00, -.352959E-01, .000000E+00, .000000E+00}},
    {1.000, .0, -.61788, {.208195E+00, -.334187E-01, .000000E+00, .000000E+00}},
    {1.000, .0, -.45802, {.312367E+00, -.266875E-01, .000000E+00, .000000E+00}},
    {1.000, .0, -.28160, {.484328E+00, -.143500E-01, .000000E+00, .000000E+00}},
    {1.000, .0, -.09501, {.807490E+00, .153739E-03, .000000E+00, .000000E+00}},
    {1.000, .0, .09501, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, .28160, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, .45802, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, .61788, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, .75540, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, .86563, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, .94458, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, .0, .98940, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, -.98940, {.553862E-01, .219558E-01, .181689E-02, .204991E-04}},
    {1.000, 90.0, -.94458, {.585238E-01, .237238E-01, .429232E-02, .441628E-04}},
    {1.000, 90.0, -.86563, {.645430E-01, .271399E-01, .708604E-02, .611213E-04}},
    {1.000, 90.0, -.75540, {.741562E-01, .326561E-01, .104047E-01, .674601E-04}},
    {1.000, 90.0, -.61788, {.885759E-01, .410670E-01, .145388E-01, .595037E-04}},
    {1.000, 90.0, -.45802, {.109780E+00, .537914E-01, .199176E-01, .334282E-04}},
    {1.000, 90.0, -.28160, {.140755E+00, .735436E-01, .272017E-01, -.161735E-04}},
    {1.000, 90.0, -.09501, {.184828E+00, .106090E+00, .375727E-01, -.100117E-03}},
    {1.000, 90.0, .09501, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, .28160, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, .45802, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, .61788, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, .75540, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, .86563, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, .94458, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 90.0, .98940, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, -.98940, {.514661E-01, -.175194E-01, .353333E-08, -.241748E-11}},
    {1.000, 180.0, -.94458, {.511342E-01, -.122200E-01, .334104E-08, -.702084E-11}},
    {1.000, 180.0, -.86563, {.544662E-01, -.681647E-02, .323669E-08, -.125962E-10}},
    {1.000, 180.0, -.75540, {.612331E-01, -.129271E-02, .316845E-08, -.180448E-10}},
    {1.000, 180.0, -.61788, {.716076E-01, .436260E-02, .306566E-08, -.218610E-10}},
    {1.000, 180.0, -.45802, {.861062E-01, .100360E-01, .283014E-08, -.217545E-10}},
    {1.000, 180.0, -.28160, {.105191E+00, .151938E-01, .232724E-08, -.133365E-10}},
    {1.000, 180.0, -.09501, {.127362E+00, .171386E-01, .139210E-08, .139267E-10}},
    {1.000, 180.0, .09501, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, .28160, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, .45802, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, .61788, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, .75540, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, .86563, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, .94458, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {1.000, 180.0, .98940, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -2.00000, {.274800E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, 2.00000, {.273970E+00, .317170E-01, .000000E+00, .000000E+00}},
    {.000, .0, -.98940, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.94458, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.86563, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.75540, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.61788, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.45802, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.28160, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.09501, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, .09501, {.934986E-01, .567385E-02, .000000E+00, .000000E+00}},
    {.000, .0, .28160, {.161225E+00, .969435E-02, .000000E+00, .000000E+00}},
    {.000, .0, .45802, {.184881E+00, .668130E-02, .000000E+00, .000000E+00}},
    {.000, .0, .61788, {.171583E+00, -.814586E-03, .000000E+00, .000000E+00}},
    {.000, .0, .75540, {.145162E+00, -.920166E-02, .000000E+00, .000000E+00}},
    {.000, .0, .86563, {.117483E+00, -.162311E-01, .000000E+00, .000000E+00}},
    {.000, .0, .94458, {.936345E-01, -.209165E-01, .000000E+00, .000000E+00}},
    {.000, .0, .98940, {.753632E-01, -.231434E-01, .000000E+00, .000000E+00}},
    {.000, 90.0, -.98940, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.94458, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.86563, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.75540, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.61788, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.45802, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.28160, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.09501, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, .09501, {.561758E-01, .155300E-01, .949841E-02, .140749E-04}},
    {.000, 90.0, .28160, {.764336E-01, .267527E-01, .127407E-01, -.974867E-05}},
    {.000, 90.0, .45802, {.815112E-01, .307685E-01, .128460E-01, -.168168E-04}},
    {.000, 90.0, .61788, {.788949E-01, .301894E-01, .111981E-01, -.116161E-04}},
    {.000, 90.0, .75540, {.744413E-01, .282078E-01, .896557E-02, -.330049E-05}},
    {.000, 90.0, .86563, {.702974E-01, .261542E-01, .657858E-02, .281815E-05}},
    {.000, 90.0, .94458, {.672162E-01, .245569E-01, .418008E-02, .474342E-05}},
    {.000, 90.0, .98940, {.654567E-01, .236249E-01, .181520E-02, .287773E-05}},
    {.000, 180.0, -.98940, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.94458, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.86563, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.75540, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.61788, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.45802, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.28160, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.09501, {.872073E-02, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, .09501, {.488879E-01, .841337E-03, -.132059E-08, -.213737E-10}},
    {.000, 180.0, .28160, {.600969E-01, .584603E-03, -.215189E-08, -.129592E-10}},
    {.000, 180.0, .45802, {.597792E-01, -.317858E-02, -.295357E-08, -.429956E-11}},
    {.000, 180.0, .61788, {.555423E-01, -.804907E-02, -.352726E-08, .214458E-12}},
    {.000, 180.0, .75540, {.518829E-01, -.128188E-01, -.387878E-08, .153983E-11}},
    {.000, 180.0, .86563, {.504972E-01, -.170131E-01, -.405942E-08, .118486E-11}},
    {.000, 180.0, .94458, {.522332E-01, -.203676E-01, -.412493E-08, .381463E-12}},
    {.000, 180.0, .98940, {.577439E-01, -.226181E-01, -.412231E-08, -.802475E-13}},
};

// 3rdparty/polradtran/runtesta, testa.out.check
constexpr table_row testa_check[] = {
    {15.000, .0, -2.00000, {.173477E+01, .103900E+00, .000000E+00, .000000E+00}},
    {15.000, .0, 2.00000, {.500000E+01, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, .0, -.96029, {.415400E+00, -.188827E+00, .000000E+00, .000000E+00}},
    {15.000, .0, -.79667, {.451745E+00, -.219329E+00, .000000E+00, .000000E+00}},
    {15.000, .0, -.52553, {.579206E+00, -.203552E+00, .000000E+00, .000000E+00}},
    {15.000, .0, -.18343, {.799029E+00, -.137430E+00, .000000E+00, .000000E+00}},
    {15.000, .0, .18343, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, .0, .52553, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, .0, .79667, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, .0, .96029, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 90.0, -.96029, {.443776E+00, .145556E+00, .508178E-01, .100812E-03}},
    {15.000, 90.0, -.79667, {.465895E+00, .155671E+00, .123671E+00, .172961E-03}},
    {15.000, 90.0, -.52553, {.506812E+00, .180072E+00, .215800E+00, .119817E-03}},
    {15.000, 90.0, -.18343, {.572843E+00, .230011E+00, .347868E+00, .250606E-04}},
    {15.000, 90.0, .18343, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 90.0, .52553, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 90.0, .79667, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 90.0, .96029, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 180.0, -.96029, {.496458E+00, -.882229E-01, .203017E-07, -.120093E-10}},
    {15.000, 180.0, -.79667, {.608726E+00, -.173126E-01, .124362E-07, -.290620E-10}},
    {15.000, 180.0, -.52553, {.750658E+00, .267501E-01, .488537E-09, -.297141E-10}},
    {15.000, 180.0, -.18343, {.917189E+00, -.951650E-02, -.209440E-07, -.424919E-11}},
    {15.000, 180.0, .18343, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 180.0, .52553, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 180.0, .79667, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {15.000, 180.0, .96029, {.000000E+00, .000000E+00, .000000E+00, .000000E+00}},
    {5.000, .0, -2.00000, {.124860E+01, .450068E-01, .000000E+00, .000000E+00}},
    {5.000, .0, 2.00000, {.241529E+01, .820379E-01, .000000E+00, .000000E+00}},
    {5.000, .0, -.96029, {.316952E+00, -.571097E-01, .000000E+00, .000000E+00}},
    {5.000, .0, -.79667, {.378777E+00, -.746844E-01, .000000E+00, .000000E+00}},
    {5.000, .0, -.52553, {.528770E+00, -.771516E-01, .000000E+00, .000000E+00}},
    {5.000, .0, -.18343, {.836870E+00, -.524060E-01, .000000E+00, .000000E+00}},
    {5.000, .0, .18343, {.640948E+00, .851421E-02, .000000E+00, .000000E+00}},
    {5.000, .0, .52553, {.556200E+00, .221156E-01, .000000E+00, .000000E+00}},
    {5.000, .0, .79667, {.420970E+00, -.145090E-01, .000000E+00, .000000E+00}},
    {5.000, .0, .96029, {.319709E+00, -.692199E-01, .000000E+00, .000000E+00}},
    {5.000, 90.0, -.96029, {.310467E+00, .406900E-01, .158511E-01, .199327E-03}},
    {5.000, 90.0, -.79667, {.336355E+00, .480606E-01, .404127E-01, .389393E-03}},
    {5.000, 90.0, -.52553, {.386372E+00, .658319E-01, .752772E-01, .381569E-03}},
    {5.000, 90.0, -.18343, {.460334E+00, .100630E+00, .127413E+00, .373166E-04}},
    {5.000, 90.0, .18343, {.438575E+00, .148283E+00, .202021E+00, .428581E-04}},
    {5.000, 90.0, .52553, {.362782E+00, .134683E+00, .157325E+00, .241136E-04}},
    {5.000, 90.0, .79667, {.301335E+00, .118329E+00, .922815E-01, .131270E-04}},
    {5.000, 90.0, .96029, {.274216E+00, .111656E+00, .382001E-01, .529811E-05}},
    {5.000, 180.0, -.96029, {.313108E+00, -.204385E-01, .540970E-08, -.240508E-10}},
    {5.000, 180.0, -.79667, {.346697E+00, .166653E-02, .343395E-08, -.676079E-10}},
    {5.000, 180.0, -.52553, {.398128E+00, .151647E-01, .481160E-09, -.106247E-09}},
    {5.000, 180.0, -.18343, {.453772E+00, .758518E-02, -.501378E-08, -.967635E-10}},
    {5.000, 180.0, .18343, {.566832E+00, -.656012E-01, -.231479E-07, -.374676E-11}},
    {5.000, 180.0, .52553, {.390841E+00, -.143244E+00, -.278118E-07, -.210807E-11}},
    {5.000, 180.0, .79667, {.273935E+00, -.161544E+00, -.256513E-07, -.114760E-11}},
    {5.000, 180.0, .96029, {.246343E+00, -.142586E+00, -.223435E-07, -.463174E-12}},
    {.000, .0, -2.00000, {.262221E+00, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, 2.00000, {.510073E+00, -.474961E-02, .000000E+00, .000000E+00}},
    {.000, .0, -.96029, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.79667, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.52553, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, -.18343, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, .0, .18343, {.124286E+00, -.409200E-02, .000000E+00, .000000E+00}},
    {.000, .0, .52553, {.151508E+00, -.428797E-02, .000000E+00, .000000E+00}},
    {.000, .0, .79667, {.177364E+00, -.560435E-02, .000000E+00, .000000E+00}},
    {.000, .0, .96029, {.191069E+00, -.970279E-02, .000000E+00, .000000E+00}},
    {.000, 90.0, -.96029, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.79667, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.52553, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, -.18343, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 90.0, .18343, {.121056E+00, -.102552E-02, .196468E-02, -.126150E-04}},
    {.000, 90.0, .52553, {.146406E+00, .221397E-02, .336247E-02, -.516552E-04}},
    {.000, 90.0, .79667, {.170769E+00, .742230E-02, .428931E-02, -.834948E-04}},
    {.000, 90.0, .96029, {.187139E+00, .119689E-01, .276823E-02, -.572090E-04}},
    {.000, 180.0, -.96029, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.79667, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.52553, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, -.18343, {.825162E-01, .000000E+00, .000000E+00, .000000E+00}},
    {.000, 180.0, .18343, {.122545E+00, -.628812E-02, -.457769E-09, -.108266E-10}},
    {.000, 180.0, .52553, {.146764E+00, -.923285E-02, -.103255E-08, -.957088E-11}},
    {.000, 180.0, .79667, {.168759E+00, -.129317E-01, -.181443E-08, -.255211E-11}},
    {.000, 180.0, .96029, {.184590E+00, -.147828E-01, -.235030E-08, .252312E-11}},
};

Matrix legendre(std::initializer_list<std::array<Numeric, 6>> rows) {
  Matrix m(static_cast<Index>(rows.size()), 6);
  Index  l = 0;
  for (const auto& r : rows) {
    for (Index k = 0; k < 6; k++) m[l, k] = r[k];
    l++;
  }
  return m;
}

//! mietest.sca of runmietest and runtesta (columns F11, F12, F33, F34, F22, F44)
Matrix mie_legendre() {
  return legendre({
      {1.00000000, -.32071711, .71206342, -.01882245, 1.00000000, .71206342},
      {1.45529318, -.20350675, 1.76014119, -.04725108, 1.45529318, 1.76014119},
      {1.05402631, .24638948, 1.06682431, .00894436, 1.05402631, 1.06682431},
      {.39758994, .18605748, .39651104, .04505815, .39758994, .39651104},
      {.11659302, .07124848, .09576412, .00958275, .11659302, .09576412},
      {.02387477, .01700757, .01765088, .00215761, .02387477, .01765088},
      {.00395010, .00302534, .00261549, .00029195, .00395010, .00261549},
      {.00053888, .00043592, .00032713, .00003502, .00053888, .00032713},
      {.00006372, .00005326, .00003583, .00000337, .00006372, .00003583},
      {.00000667, .00000572, .00000351, .00000029, .00000667, .00000351},
      {.00000063, .00000055, .00000031, .00000002, .00000063, .00000031},
      {.00000006, .00000005, .00000003, .00000000, .00000006, .00000003},
  });
}

//! rayleigh.sca of runtesta
Matrix rayleigh_legendre() {
  return legendre({{1.0, -0.5, 0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 1.5, 0.0, 0.0, 1.5}, {0.5, 0.5, 0.0, 0.0, 0.5, 0.0}});
}

//! Henyey-Greenstein F11 = (2 l + 1) g^l to degree nleg; F22 = F11, the rest 0
Matrix henyey_greenstein(Numeric g, Index nleg) {
  Matrix m(nleg + 1, 6, 0.0);
  for (Index l = 0; l <= nleg; l++) m[l, 0] = m[l, 4] = static_cast<Numeric>(2 * l + 1) * std::pow(g, l);
  return m;
}

//! rt3.f USER_INPUT turns the solar zenith angle into DIRECT_MU with a truncated pi / 180
Numeric rt3_direct_mu(Numeric theta_deg) { return std::abs(std::cos(0.017453292 * theta_deg)); }

//! The frequency of a wavelength [um], and RT3's per-micrometre to per-Hz factor lambda[um] / f
struct spectral {
  Numeric frequency, per_um_to_per_hz;
};

spectral at_wavelength(Numeric lambda_um) {
  const Numeric f = Constant::c / (lambda_um * 1e-6);
  return {.frequency = f, .per_um_to_per_hz = lambda_um / f};
}

//! Half a unit in the 6th significant digit of a value printed by rt3.f
//! (E13.6); 0 for a printed 0, which is an exact REAL*4 zero
Numeric print_half_unit(Numeric v) {
  if (v == 0.0) return 0.0;
  return 0.5 * std::pow(10.0, std::floor(std::log10(std::abs(v))) + 1.0 - 6.0);
}

/** Bound on the REAL*4 error of rt3.f's OUTPUT_FILE for one radiance entry.
 *
 * OUTPUT_FILE sums sum_m c_m t(m phi), t = cos for I, Q and sin for U, V,
 * in single precision: PHI = 3.1415927 * PHID / 180 is REAL*4, M*PHI is
 * rounded to REAL*4, the cosine is a REAL*4 function and the running sum
 * OUT is REAL*4.  Per term this gives at most
 *   |c_m| (m |phi_f - phi| + 2^-24 m phi + 2^-23),
 * and the running sum one rounding of at most 2^-24 sum_m |c_m| per term.
 */
Numeric single_precision_bound(const Tensor4& c, Index level, Index stream, Index stokes, Numeric phi_deg) {
  constexpr Numeric eps32 = 0x1p-24;
  const Numeric     phi   = phi_deg * pi / 180.0;
  const float       phif  = static_cast<float>(pi) * static_cast<float>(phi_deg) / 180.0f;
  const Numeric     dphi  = std::abs(static_cast<Numeric>(phif) - phi);
  Numeric           terms = 0.0, sum = 0.0;
  for (Index m = 0; m < c.extent(1); m++) {
    const Numeric a   = std::abs(c[level, m, stream, stokes]);
    const auto    mm  = static_cast<Numeric>(m);
    terms            += a * (mm * dphi + eps32 * mm * phi + 2.0 * eps32);
    sum              += a;
  }
  return terms + static_cast<Numeric>(c.extent(1)) * eps32 * sum;
}

/** OUTPUT_FILE's own REAL*4 evaluation of one entry, as the original
 *  program computes it: PHI = PI * PHID / 180 with a REAL*4 PI, COS(M*PHI)
 *  or SIN(M*PHI) in REAL*4, the product with the REAL*8 coefficient in
 *  REAL*8, and the running sum OUT in REAL*4.  Diagnostic only: the result
 *  depends on the platform's single-precision cosine. */
float output_file_real4(const Tensor4& c, Index level, Index stream, Index stokes, Numeric phi_deg, Numeric scale) {
  const float pi32 = static_cast<float>(pi);
  const float phi  = pi32 * static_cast<float>(phi_deg) / 180.0f;
  float       out  = 0.0f;
  for (Index m = 0; m < c.extent(1); m++) {
    const float arg = static_cast<float>(m) * phi;
    const float t   = stokes < 2 ? std::cos(arg) : std::sin(arg);
    out = static_cast<float>(static_cast<Numeric>(out) + static_cast<Numeric>(t) * c[level, m, stream, stokes] * scale);
  }
  return out;
}

//! One table entry evaluated from a result
struct table_entry {
  Numeric got{}, ref{}, bound{};
  float   real4{};
  Index   stokes{};
  bool    flux{};
};

/** Evaluates a result at every entry of an rt3.f output table.
 *
 * height gives the Z of each level; to_per_um converts the wrapper's
 * per-Hz output back to RT3's per-micrometre units.  bound is the REAL*4
 * error bound of OUTPUT_FILE (fluxes: one REAL*4 rounding, SNGL).
 */
std::vector<table_entry> evaluate_table(const rt3::result&         r,
                                        const Vector&              height,
                                        std::span<const table_row> rows,
                                        Numeric                    to_per_um) {
  std::vector<table_entry> out;
  for (const auto& row : rows) {
    Index level = -1;
    for (Index l = 0; l < size(height); l++)
      if (std::abs(height[l] - row.z) < 5e-4) level = l;
    require(level >= 0, std::format("no level at Z = {}", row.z));

    const bool flux = std::abs(row.mu) > 1.5;
    const bool upw  = row.mu < 0.0;
    Index      i    = -1;
    if (not flux) {
      for (Index k = 0; k < size(r.mu); k++)
        if (std::abs(r.mu[k] - std::abs(row.mu)) < 5.000001e-6) i = k;
      require(i >= 0, std::format("no stream at printed mu = {}", row.mu));
    }
    const Tensor4 rad = flux ? Tensor4{} : rt3::azimuth_radiance(upw ? r.up : r.down, Vector{row.phi_deg * pi / 180.0});

    for (Index s = 0; s < r.up.extent(3); s++) {
      table_entry e{.ref = row.iquv[s], .stokes = s, .flux = flux};
      if (flux) {
        e.got   = (upw ? r.up_flux : r.down_flux)[level, s] * to_per_um;
        e.bound = 0x1p-24 * std::abs(e.got);
        e.real4 = static_cast<float>(e.got);
      } else {
        e.got   = rad[level, 0, i, s] * to_per_um;
        e.bound = single_precision_bound(upw ? r.up : r.down, level, i, s, row.phi_deg) * to_per_um;
        e.real4 = output_file_real4(upw ? r.up : r.down, level, i, s, row.phi_deg, to_per_um);
      }
      out.push_back(e);
    }
  }
  return out;
}

struct table_deviation {
  Numeric worst{};        //!< max |RT3 - table| / tolerance
  Numeric worst_iq{};     //!< max |RT3 - table| / print half-unit, I and Q
  Numeric worst_uv{};     //!< the same for U and V where the REAL*4 bound is below the print half-unit
  Numeric real4{};        //!< max |REAL*4 emulation - table| / print half-unit, entries that are not noise
  Numeric real4_noise{};  //!< max relative |REAL*4 emulation - table| of the noise entries
  Index   noise{};        //!< entries whose double-precision value is below their REAL*4 bound (zero by symmetry)
  Index   entries{};
};

/** The tolerance of each entry is half a unit in its last printed digit
 *  plus the REAL*4 error bound of OUTPUT_FILE. */
table_deviation compare_table(const std::vector<table_entry>& entries) {
  table_deviation d;
  for (const auto& e : entries) {
    const Numeric half = print_half_unit(e.ref);
    const Numeric dev  = std::abs(e.got - e.ref);
    const Numeric tol  = half + e.bound;
    d.worst            = std::max(d.worst, dev / tol);
    if (e.stokes < 2 and half > 0) d.worst_iq = std::max(d.worst_iq, dev / half);
    if (e.stokes >= 2 and half > e.bound) d.worst_uv = std::max(d.worst_uv, dev / half);
    const Numeric dev4 = std::abs(static_cast<Numeric>(e.real4) - e.ref);
    if (e.ref != 0.0 and std::abs(e.got) <= e.bound) {
      d.noise++;
      d.real4_noise = std::max(d.real4_noise, dev4 / std::abs(e.ref));
    } else {
      d.real4 = std::max(d.real4, half > 0.0 ? dev4 / half : (dev4 > 0.0 ? 1e99 : 0.0));
    }
    d.entries++;
    if (dev > tol)
      std::cout << std::format(
          "    Stokes {}: RT3 {:.6e}, table {:.6e}, tolerance {:.2e}\n", e.stokes, e.got, e.ref, tol);
  }
  return d;
}

void report_table(std::string_view name, const table_deviation& d) {
  std::cout << std::format(
      "{:<64} {} entries, max |dev| / tolerance {:.2f}; in print half-units: I, Q {:.2f}, U, V {:.2f}\n",
      name,
      d.entries,
      d.worst,
      d.worst_iq,
      d.worst_uv);
  std::cout << std::format(
      "    OUTPUT_FILE's REAL*4 sum emulated: max |dev| {:.2f} print half-units; the {} entries that are "
      "REAL*4 noise (zero by symmetry) to {:.1e} relative\n",
      d.real4,
      d.noise,
      d.real4_noise);
  require(d.worst <= 1.0, std::format("{}: deviation exceeds the tolerance", name));
}

/** (a) RT3's quadratures against their defining exactness, as for RT4. */
void test_quadrature() {
  Numeric worst = 0.0;
  for (Index n : {1, 2, 5, 8, 16}) {
    const auto D = rt3::get_quadrature(n, rt3::quadrature_type::double_gauss);
    const auto G = rt3::get_quadrature(n, rt3::quadrature_type::gauss);
    const auto L = rt3::get_quadrature(n, rt3::quadrature_type::lobatto);
    for (const auto* q : {&D, &G, &L}) {
      require(size(q->mu) == n and size(q->weights) == n, "quadrature size");
      for (Index i = 0; i < n; i++)
        require(q->mu[i] > 0 and q->mu[i] <= 1 and (i == 0 or q->mu[i] > q->mu[i - 1]),
                "quadrature nodes must be ascending in (0, 1]");
    }
    const auto moment = [](const rt3::quadrature& q, Index k) {
      Numeric s = 0.0;
      for (Index i = 0; i < size(q.mu); i++) s += q.weights[i] * std::pow(q.mu[i], k);
      return s;
    };
    for (Index k = 0; k <= 2 * n - 1; k++)
      worst = std::max(worst, std::abs(moment(D, k) - 1.0 / static_cast<Numeric>(k + 1)));
    for (Index k = 0; 2 * k <= 4 * n - 1; k++)
      worst = std::max(worst, std::abs(moment(G, 2 * k) - 1.0 / static_cast<Numeric>(2 * k + 1)));
    for (Index k = 0; 2 * k <= 4 * n - 3; k++)
      worst = std::max(worst, std::abs(moment(L, 2 * k) - 1.0 / static_cast<Numeric>(2 * k + 1)));
    require(L.mu[n - 1] == 1.0, "Lobatto must include mu = 1");
  }
  std::cout << std::format("{:<64} max moment error {:9.3e}\n", "(a) quadrature exactness G/D/L", worst);
  require(worst < 1e-13, "quadrature moments are not exact");
}

/** (b) Evans' Mie benchmark (runmietest): tau = 1, omega = 0.99, the Mie
 *  series of the paper, Lambertian albedo 0.1, solar flux 0.2 pi on the
 *  horizontal at mu0 = 0.2, G quadrature with 8 nodes, aziorder 8, I Q U V
 *  at the top (Z = 1) and bottom (Z = 0) at azimuths 0, 90 and 180 deg.
 *  The table is the output of the original program (per micrometre). */
void test_mietest() {
  const auto   sp = at_wavelength(0.951);
  rt3::problem p;
  p.nstokes                = 4;
  p.nmu                    = 8;
  p.quad                   = rt3::quadrature_type::gauss;
  p.aziorder               = 8;
  p.delta_m                = false;
  p.max_delta_tau          = 1e-6;  // rt3.f's MAX_DELTA_TAU
  p.direct_flux            = 0.628318531 * sp.per_um_to_per_hz;
  p.direct_mu              = rt3_direct_mu(78.46304097);
  p.thermal                = false;
  p.frequency              = sp.frequency;
  p.height                 = Vector{1.0, 0.0};
  p.temperature            = Vector{0.0, 0.0};
  p.gas_extinction         = Vector{0.0};
  p.scattering_sets        = {{.extinction = 1.0, .scattering = 0.99, .legendre = mie_legendre()}};
  p.layer_scattering_index = {0};
  p.sky_temperature        = 0.0;
  p.surface_temperature    = 0.0;
  p.ground                 = rt3::lambertian_surface{.albedo = 0.1};
  const auto r             = rt3::solve(p);

  Numeric uv0 = 0.0;
  for (const auto* t : {&r.up, &r.down})
    for (Index l = 0; l < t->extent(0); l++)
      for (Index i = 0; i < t->extent(2); i++)
        for (Index s = 2; s < 4; s++) uv0 = std::max(uv0, std::abs((*t)[l, 0, i, s]));
  require(uv0 == 0.0, "the m = 0 coefficients of U and V must be 0");
  report_table("(b) runmietest (Evans and Stephens 1991 Mie case) vs its table",
               compare_table(evaluate_table(r, p.height, mietest_check, 1.0 / sp.per_um_to_per_hz)));
}

//! Exact-SI radiation constants in RT3's units, 2 h c^2 [W m-2 sr-1 um^4] and h c / k [um K]
constexpr Numeric planck_c1 = 2.0 * Constant::h * Constant::c * Constant::c * 1e24;
constexpr Numeric planck_c2 = Constant::h * Constant::c / Constant::k * 1e6;

//! Evans' original 5-digit Planck function [W m-2 sr-1 um-1]
Numeric planck_5digit(Numeric lambda_um, Numeric t) {
  return 1.1911e8 / std::pow(lambda_um, 5) / (std::exp(1.4388e4 / (lambda_um * t)) - 1.0);
}

//! The temperature at which the exact Planck function equals planck_5digit(t)
Numeric exact_temperature_of_5digit(Numeric lambda_um, Numeric t) {
  const Numeric b = planck_5digit(lambda_um, t);
  return planck_c2 / (lambda_um * std::log1p(planck_c1 / (std::pow(lambda_um, 5) * b)));
}

/** (c) Evans' general benchmark (runtesta): a Rayleigh layer over a Mie
 *  layer with gas absorption, solar and thermal sources at 3 um, G
 *  quadrature with 4 nodes, aziorder 4, Lambertian albedo 0.25 at 300 K,
 *  output at all three levels.
 *
 *  The table was made with Evans' 5-digit radiation constants, which the
 *  ARTS3 RT3 replaces by exact ones; at 3 um and 300 to 200 K that raises
 *  B by 2.1e-4 to 3.4e-4 relative.  RT3 evaluates the Planck function only at
 *  the interface, surface and sky temperatures, so the test reproduces the
 *  5-digit values exactly by passing the temperatures T' at which the exact
 *  Planck function equals the 5-digit one at T.  The run with the true
 *  temperatures is reported, not checked. */
void test_testa() {
  const Numeric lambda = 3.0;
  const auto    sp     = at_wavelength(lambda);
  for (bool five_digit : {true, false}) {
    const auto   t = [&](Numeric x) { return five_digit ? exact_temperature_of_5digit(lambda, x) : x; };
    rt3::problem p;
    p.nstokes                = 4;
    p.nmu                    = 4;
    p.quad                   = rt3::quadrature_type::gauss;
    p.aziorder               = 4;
    p.delta_m                = false;
    p.max_delta_tau          = 1e-6;
    p.direct_flux            = 5.0 * sp.per_um_to_per_hz;
    p.direct_mu              = rt3_direct_mu(60.0);
    p.thermal                = true;
    p.frequency              = sp.frequency;
    p.height                 = Vector{15.0, 5.0, 0.0};
    p.temperature            = Vector{t(200.0), t(270.0), t(300.0)};
    p.gas_extinction         = Vector{0.02, 0.05};
    p.scattering_sets        = {{.extinction = 0.05, .scattering = 0.05, .legendre = rayleigh_legendre()},
                                {.extinction = 1.0, .scattering = 0.99, .legendre = mie_legendre()}};
    p.layer_scattering_index = {0, 1};
    p.sky_temperature        = 0.0;
    p.surface_temperature    = t(300.0);
    p.ground                 = rt3::lambertian_surface{.albedo = 0.25};
    const auto r             = rt3::solve(p);
    const auto entries       = evaluate_table(r, p.height, testa_check, 1.0 / sp.per_um_to_per_hz);
    if (five_digit) {
      report_table("(c) runtesta, 5-digit Planck values via T', vs its table", compare_table(entries));
    } else {
      Numeric worst = 0.0;
      for (const auto& e : entries)
        if (e.stokes == 0 and e.ref != 0.0) worst = std::max(worst, std::abs(e.got - e.ref) / e.ref);
      std::cout << std::format(
          "    with the true temperatures (exact Planck constants) I differs from the table "
          "by up to {:.1e} relative\n",
          worst);
    }
  }
}

/** Exact transfer through a gas layer of optical path x along a stream,
 *  with the Planck function linear in optical depth from b_start to b_end:
 *    I = I0 e^-x + b_end (1 - e^-x) - (b_end - b_start) (1 - (1 + x) e^-x) / x,
 *  and polarization only attenuated. */
rtepack::stokvec gas_path(const rtepack::stokvec& v0, Numeric x, Numeric b_start, Numeric b_end) {
  const Numeric    ex = std::exp(-x);
  rtepack::stokvec v{};
  v[0] = v0[0] * ex + b_end * (-std::expm1(-x)) - (b_end - b_start) * (1.0 - (1.0 + x) * ex) / x;
  for (Index s = 1; s < 4; s++) v[s] = v0[s] * ex;
  return v;
}

//! |r_v|^2, |r_h|^2 and r_v r_h* from Snell's law, index 1 above n
struct fresnel_coefficients {
  Numeric rv2, rh2;
  Complex rvrh;
};

fresnel_coefficients fresnel(Complex n, Numeric mu) {
  const Complex cos_t = std::sqrt(1.0 - (1.0 - mu * mu) / (n * n));
  const Complex rv    = (n * mu - cos_t) / (n * mu + cos_t);
  const Complex rh    = (mu - n * cos_t) / (mu + n * cos_t);
  return {std::norm(rv), std::norm(rh), rv * std::conj(rh)};
}

/** (d) Gas-only atmosphere, thermal source, every level and stream, for
 *  nstokes 1 to 4: over a Lambertian surface (double_gauss, whose discrete
 *  Lambertian operator 2 A mu_j w_j conserves energy) and over a Fresnel
 *  surface (gauss with two extra angles).  RT3 solves gas-only layers
 *  analytically, so only round-off remains.  The m > 0 modes and U, V must
 *  vanish exactly, and the fluxes must be 2 pi sum_i w_i mu_i I_i of the
 *  closed form. */
void test_gas_only() {
  const Vector  height{3000.0, 2000.0, 1000.0, 0.0};
  const Vector  temperature{220.0, 240.0, 265.0, 285.0};
  const Vector  gas{1e-4, 3e-4, 5e-4};  // tau = 0.1, 0.3, 0.5
  const Numeric sky = Constant::cosmic_microwave_background_temperature, tsurf = 290.0;
  const Index   nlay = 3;

  for (Numeric frequency : {50e9, 30e12}) {
    for (bool lambert : {true, false}) {
      for (Index ns = 1; ns <= 4; ns++) {
        constexpr Numeric A = 0.3;
        const Complex     n{3.0, 0.2};
        rt3::problem      p;
        p.nstokes                = ns;
        p.nmu                    = 8;
        p.quad                   = lambert ? rt3::quadrature_type::double_gauss : rt3::quadrature_type::gauss;
        p.extra_mu               = lambert ? Vector{} : Vector{0.45, 1.0};
        p.aziorder               = 2;
        p.thermal                = true;
        p.frequency              = frequency;
        p.height                 = height;
        p.temperature            = temperature;
        p.gas_extinction         = gas;
        p.layer_scattering_index = ArrayOfIndex(nlay, -1);
        p.sky_temperature        = sky;
        p.surface_temperature    = tsurf;
        if (lambert)
          p.ground = rt3::lambertian_surface{.albedo = A};
        else
          p.ground = rt3::fresnel_surface{.refractive_index = n};
        const auto r   = rt3::solve(p);
        const auto nmu = size(r.mu);

        const auto              B  = [&](Numeric t) { return planck(frequency, t); };
        const auto              dz = [&](Index l) { return std::abs(height[l] - height[l + 1]); };
        rtepack::stokvec_matrix dn(nlay + 1, nmu), up(nlay + 1, nmu);
        for (Index i = 0; i < nmu; i++) {
          dn[0, i] = {B(sky), 0, 0, 0};
          for (Index l = 0; l < nlay; l++)
            dn[l + 1, i] = gas_path(dn[l, i], gas[l] * dz(l) / r.mu[i], B(temperature[l]), B(temperature[l + 1]));
        }
        Numeric flux_down = 0.0;
        for (Index j = 0; j < nmu; j++) flux_down += r.weights[j] * r.mu[j] * dn[nlay, j][0];
        for (Index i = 0; i < nmu; i++) {
          if (lambert) {
            up[nlay, i] = {(1 - A) * B(tsurf) + 2 * A * flux_down, 0, 0, 0};
          } else {
            const auto    fc = fresnel(n, r.mu[i]);
            const Numeric r1 = 0.5 * (fc.rv2 + fc.rh2), r2 = 0.5 * (fc.rv2 - fc.rh2);
            up[nlay, i] = {(1 - r1) * B(tsurf) + r1 * dn[nlay, i][0], r2 * (dn[nlay, i][0] - B(tsurf)), 0, 0};
          }
          for (Index l = nlay - 1; l >= 0; l--)
            up[l, i] = gas_path(up[l + 1, i], gas[l] * dz(l) / r.mu[i], B(temperature[l + 1]), B(temperature[l]));
        }

        Numeric dev = 0.0, dev_flux = 0.0;
        bool    zero = true;
        for (Index l = 0; l <= nlay; l++) {
          Numeric fu = 0.0, fd = 0.0;
          for (Index i = 0; i < nmu; i++) {
            for (Index s = 0; s < ns; s++) {
              dev = std::max({dev,
                              std::abs(r.up[l, 0, i, s] - up[l, i][s]) / up[l, i][0],
                              std::abs(r.down[l, 0, i, s] - dn[l, i][s]) / dn[l, i][0]});
              for (Index m = 1; m <= p.aziorder; m++) zero = zero and r.up[l, m, i, s] == 0 and r.down[l, m, i, s] == 0;
              if (s >= 2) zero = zero and r.up[l, 0, i, s] == 0 and r.down[l, 0, i, s] == 0;
            }
            fu += 2 * pi * r.weights[i] * r.mu[i] * up[l, i][0];
            fd += 2 * pi * r.weights[i] * r.mu[i] * dn[l, i][0];
          }
          dev_flux = std::max({dev_flux, std::abs(r.up_flux[l, 0] - fu) / fu, std::abs(r.down_flux[l, 0] - fd) / fd});
        }
        if (not lambert and ns > 1) {
          // Oblique emission from a warm dielectric is vertically polarized
          for (Index i = 0; i < nmu; i++)
            if (r.mu[i] < 0.9) zero = zero and r.up[nlay, 0, i, 1] > 0.0;
        }
        require(zero, "gas-only thermal: m > 0, U and V must vanish and Fresnel emission must have Q > 0");
        if (ns == 4 or ns == 1)
          std::cout << std::format("{:<64} max rel dev {:9.3e}, flux {:9.3e}\n",
                                   std::format("(d) gas-only thermal, {:.0e} Hz, {}, nstokes {}",
                                               frequency,
                                               lambert ? "Lambertian A = 0.3 (D)" : "Fresnel 3+0.2i (G + 2 extra)",
                                               ns),
                                   dev,
                                   dev_flux);
        require(dev < 1e-12 and dev_flux < 1e-12, "gas-only thermal closed form");
      }
    }
  }
}

//! Propagation direction in a right-handed frame with z up; mu_z > 0 upward
Vector3 direction(Numeric mu_z, Numeric phi) {
  const Numeric s = std::sqrt(1.0 - mu_z * mu_z);
  return {s * std::cos(phi), s * std::sin(phi), mu_z};
}

//! Meridional basis of k: h = k x z / |k x z|, v = h x k
std::pair<Vector3, Vector3> meridional(const Vector3& k) {
  Vector3 h  = cross(k, Vector3{0.0, 0.0, 1.0});
  h         /= std::sqrt(dot(h, h));
  return {cross(h, k), h};
}

/** Rayleigh (dipole) scattering of unpolarized light from k_in into k_out,
 *  as the Stokes vector of the scattered field in the meridional basis of
 *  k_out: unpolarized light is the incoherent sum of two orthogonal linear
 *  polarizations e (the v and h of k_in), the dipole radiates the part of
 *  e perpendicular to k_out, and E_v = v . e, E_h = h . e, so
 *    I = sum (E_v^2 + E_h^2), Q = sum (E_v^2 - E_h^2), U = sum 2 E_v E_h,
 *    V = 0,
 *  scaled by 3/4 so that I = 3/4 (1 + cos^2 Theta) is the phase function
 *  normalised to 1 over 4 pi. */
rtepack::stokvec rayleigh_column(const Vector3& k_out, const Vector3& k_in) {
  const auto [vi, hi] = meridional(k_in);
  const auto [vo, ho] = meridional(k_out);
  rtepack::stokvec z{};
  for (const auto& e : {vi, hi}) {
    const Numeric ev = dot(vo, e), eh = dot(ho, e);
    z[0] += 0.75 * (ev * ev + eh * eh);
    z[1] += 0.75 * (ev * ev - eh * eh);
    z[2] += 0.75 * 2.0 * ev * eh;
  }
  return z;
}

/** (e) Single scattering by an optically thin, conservative Rayleigh layer
 *  (tau = 1e-5 and 1e-6) over a black surface, lit by the direct beam (flux F on
 *  the horizontal, mu0 = 0.6, propagating toward phi = 0).  The reference
 *  is the exact single-scattering solution
 *    up at the top:     (F / 4 pi) Z (1 - exp(-tau (1/mu0 + 1/mu))) / (mu0 + mu)
 *    down at the bottom: (F / 4 pi) Z (exp(-tau/mu) - exp(-tau/mu0)) / (mu - mu0)
 *  with the Stokes column Z of rayleigh_column(), evaluated at 8 azimuths
 *  for the 8 gauss streams and the extra angles 0.5 and 0.77.  It fixes
 *  RT3's azimuth origin and sense, its Q and its U sign independently of
 *  RT3's rotation formulas.  The neglected multiple scattering is of
 *  relative order tau, and RT3's first-order initial sublayer adds
 *  max_delta_tau / mu_min = 5e-8; the tolerance is 10 tau, and the
 *  deviation must shrink with tau, so that it is the remainder and not a
 *  convention error (which would be of order 1 in Q or U). */
void test_single_scattering() {
  const Numeric mu0 = 0.6, F = 2.5;
  const Vector  phis{0.0, 30.0, 60.0, 90.0, 135.0, 180.0, 250.0, 315.0};
  Vector        phi(phis.size());
  for (Index k = 0; k < size(phis); k++) phi[k] = phis[k] * pi / 180.0;

  std::vector<Numeric> devs;
  for (Numeric tau : {1e-5, 1e-6}) {
    for (Index ns : {1, 2, 3, 4}) {
      rt3::problem p;
      p.nstokes                = ns;
      p.nmu                    = 8;
      p.quad                   = rt3::quadrature_type::gauss;
      p.extra_mu               = Vector{0.5, 0.77};
      p.aziorder               = 4;
      p.max_delta_tau          = 1e-9;
      p.direct_flux            = F;
      p.direct_mu              = mu0;
      p.thermal                = false;
      p.frequency              = Constant::c / 0.5e-6;
      p.height                 = Vector{1.0, 0.0};
      p.temperature            = Vector{0.0, 0.0};
      p.gas_extinction         = Vector{0.0};
      p.scattering_sets        = {{.extinction = tau, .scattering = tau, .legendre = rayleigh_legendre()}};
      p.layer_scattering_index = {0};
      p.ground                 = rt3::lambertian_surface{.albedo = 0.0};
      const auto r             = rt3::solve(p);
      const auto up            = rt3::azimuth_radiance(r.up, phi);
      const auto dn            = rt3::azimuth_radiance(r.down, phi);

      const Vector3 k0    = direction(-mu0, 0.0);
      Numeric       scale = 0.0, dev = 0.0, u_scale = 0.0;
      for (int pass = 0; pass < 2; pass++) {
        for (Index k = 0; k < size(phi); k++) {
          for (Index i = 0; i < size(r.mu); i++) {
            const Numeric mu = r.mu[i];
            const auto    zu = rayleigh_column(direction(mu, phi[k]), k0);
            const auto    zd = rayleigh_column(direction(-mu, phi[k]), k0);
            const Numeric fu = F / (4 * pi) * -std::expm1(-tau * (1 / mu0 + 1 / mu)) / (mu0 + mu);
            const Numeric fd = F / (4 * pi) * std::exp(-tau / mu0) * std::expm1(tau / mu0 - tau / mu) / (mu - mu0);
            for (Index s = 0; s < ns; s++) {
              const Numeric ru = fu * zu[s], rd = fd * zd[s];
              if (pass == 0) {
                scale = std::max({scale, std::abs(ru), std::abs(rd)});
                if (s == 2) u_scale = std::max(u_scale, std::abs(ru));
              } else {
                dev = std::max({dev, std::abs(up[0, k, i, s] - ru), std::abs(dn[1, k, i, s] - rd)});
              }
            }
          }
        }
      }
      dev /= scale;
      if (ns >= 3) require(u_scale > 0.1 * scale, "the single-scattering U must be large enough to pin its sign");
      std::cout << std::format(
          "{:<64} max |dev| / max |I| {:9.3e}\n",
          std::format("(e) thin Rayleigh layer, tau {:.0e}, nstokes {}, vs vector geometry", tau, ns),
          dev);
      if (ns == 4) devs.push_back(dev);
      require(dev < 10.0 * tau, "single scattering: deviation exceeds 10 tau");
    }
  }
  require(devs[1] < 0.2 * devs[0], "single scattering: the remainder must shrink with tau");
}

/** (f) Delta-M with f = 0 must be the identity.  For a series of degree
 *  < M = 2 nmu_total RT3 sets f = 0 and still raises the degree to M - 1,
 *  reading coefficients the set does not have; READ_SCAT_FILE read stale
 *  entries there (with gfortran, uninitialised stack memory for the first
 *  file), GET_SCAT_SET zeroes them.  A
 *  regression test of that ARTS3 change, with the runtesta atmosphere.
 *
 *  The two runs differ in round-off only (the series is summed at 64
 *  instead of 32 azimuths, and (2 l + 1) (c / (2 l + 1)) is not always c).
 *  RT3's initial sublayer stores T = 1 - O(max_delta_tau), so 1 - T, which
 *  carries the physics, has a relative precision of eps / max_delta_tau;
 *  that is RT3's round-off floor, and the tolerance is 10 times it. */
void test_delta_m_identity() {
  const auto   sp = at_wavelength(3.0);
  rt3::problem p;
  p.nstokes                = 4;
  p.nmu                    = 8;
  p.aziorder               = 6;
  p.direct_flux            = 5.0 * sp.per_um_to_per_hz;
  p.direct_mu              = 0.5;
  p.thermal                = true;
  p.frequency              = sp.frequency;
  p.height                 = Vector{15.0, 5.0, 0.0};
  p.temperature            = Vector{200.0, 270.0, 300.0};
  p.gas_extinction         = Vector{0.02, 0.05};
  p.scattering_sets        = {{.extinction = 0.05, .scattering = 0.05, .legendre = rayleigh_legendre()},
                              {.extinction = 1.0, .scattering = 0.99, .legendre = mie_legendre()}};
  p.layer_scattering_index = {0, 1};
  p.surface_temperature    = 300.0;
  p.ground                 = rt3::lambertian_surface{.albedo = 0.25};
  const auto plain         = rt3::solve(p);
  p.delta_m                = true;
  const auto scaled        = rt3::solve(p);

  Numeric scale = 0.0, dev = 0.0;
  for (Index x = 0; x < static_cast<Index>(plain.up.size()); x++) {
    scale = std::max({scale, std::abs(plain.up.data_handle()[x]), std::abs(plain.down.data_handle()[x])});
    dev   = std::max({dev,
                      std::abs(plain.up.data_handle()[x] - scaled.up.data_handle()[x]),
                      std::abs(plain.down.data_handle()[x] - scaled.down.data_handle()[x])});
  }
  const Numeric tol = 10.0 * std::numeric_limits<Numeric>::epsilon() / p.max_delta_tau;
  std::cout << std::format("{:<64} max |dev| / max |I| {:9.3e} (tolerance {:.1e})\n",
                           "(f) delta-M with f = 0 (degree 11 < M = 16)",
                           dev / scale,
                           tol);
  require(dev < tol * scale, "delta-M with f = 0 must not change the result");
}

/** (g) The direct beam: down_flux - 2 pi sum_i w_i mu_i c_0 must be
 *  F exp(-tau / mu0) at every level, with tau from the gas and the
 *  (delta-M scaled) particle extinction, (1 - omega f) k for a
 *  Henyey-Greenstein set with f = g^M. */
void test_direct_beam() {
  const Numeric g = 0.7, F = 3.0, mu0 = 0.35;
  for (bool dm : {false, true}) {
    rt3::problem p;
    p.nstokes                = 2;
    p.nmu                    = 8;
    p.aziorder               = 3;
    p.delta_m                = dm;
    p.direct_flux            = F;
    p.direct_mu              = mu0;
    p.thermal                = false;
    p.frequency              = Constant::c / 0.5e-6;
    p.height                 = Vector{3.0, 2.0, 1.0, 0.0};
    p.temperature            = Vector(4, 0.0);
    p.gas_extinction         = Vector{0.05, 0.0, 0.1};
    p.scattering_sets        = {{.extinction = 0.4, .scattering = 0.3, .legendre = henyey_greenstein(g, 29)}};
    p.layer_scattering_index = {-1, 0, 0};
    p.ground                 = rt3::lambertian_surface{.albedo = 0.2};
    const auto r             = rt3::solve(p);

    const Numeric f  = dm ? std::pow(g, 2 * p.nmu) : 0.0;
    const Numeric kp = (1.0 - 0.3 / 0.4 * f) * 0.4;
    const Vector  ext{0.05, kp, kp + 0.1};
    Numeric       tau = 0.0, dev = 0.0;
    for (Index l = 0; l <= 3; l++) {
      if (l > 0) tau += ext[l - 1];
      Numeric diffuse = 0.0;
      for (Index i = 0; i < size(r.mu); i++) diffuse += 2 * pi * r.weights[i] * r.mu[i] * r.down[l, 0, i, 0];
      const Numeric direct = F * std::exp(-tau / mu0);
      dev                  = std::max(dev, std::abs(r.down_flux[l, 0] - diffuse - direct) / direct);
    }
    std::cout << std::format("{:<64} max rel dev {:9.3e}\n",
                             std::format("(g) direct beam F exp(-tau / mu0){}", dm ? ", delta-M scaled tau" : ""),
                             dev);
    require(dev < 1e-13, "direct beam attenuation");
  }
}

void expect_throw(std::string_view what, std::string_view match, const std::function<void()>& f) {
  try {
    f();
  } catch (const std::exception& e) {
    const std::string_view msg{e.what()};
    require(msg.find(match) != std::string_view::npos,
            std::format("{} threw, but not with \"{}\": {}", what, match, msg));
    std::cout << std::format("    {:<50} throws: {}\n", what, msg.substr(msg.find(match), 90));
    return;
  }
  throw std::runtime_error(std::format("{} did not throw", what));
}

/** (h) Error paths; each must throw before the Fortran code (which would
 *  STOP the process or overrun a buffer) is called. */
void test_errors() {
  const auto good = [] {
    rt3::problem p;
    p.nstokes                = 4;
    p.nmu                    = 4;
    p.aziorder               = 2;
    p.direct_flux            = 1.0;
    p.direct_mu              = 0.5;
    p.thermal                = true;
    p.frequency              = Constant::c / 3e-6;
    p.height                 = Vector{2.0, 1.0, 0.0};
    p.temperature            = Vector{250.0, 260.0, 270.0};
    p.gas_extinction         = Vector{0.1, 0.0};
    p.scattering_sets        = {{.extinction = 1.0, .scattering = 0.9, .legendre = rayleigh_legendre()}};
    p.layer_scattering_index = {-1, 0};
    p.surface_temperature    = 280.0;
    p.ground                 = rt3::lambertian_surface{.albedo = 0.1};
    return p;
  };
  rt3::solve(good());
  const auto isotropic = legendre({{1.0, 0.0, 0.0, 0.0, 1.0, 0.0}});
  const auto layers    = [](rt3::problem& p, Index nlay) {
    p.height = Vector(nlay + 1, 0.0);
    for (Index i = 0; i <= nlay; i++) p.height[i] = static_cast<Numeric>(nlay - i);
    p.temperature            = Vector(nlay + 1, 250.0);
    p.gas_extinction         = Vector(nlay, 0.0);
    p.layer_scattering_index = ArrayOfIndex(nlay, -1);
  };

  expect_throw("nstokes = 0", "nstokes 1 to 4", [&] {
    auto p    = good();
    p.nstokes = 0;
    rt3::solve(p);
  });
  expect_throw("nstokes = 5", "nstokes 1 to 4", [&] {
    auto p    = good();
    p.nstokes = 5;
    rt3::solve(p);
  });
  expect_throw("nmu = 0", "at least one quadrature node", [&] {
    auto p = good();
    p.nmu  = 0;
    rt3::solve(p);
  });
  expect_throw("extra_mu with double_gauss", "only to the gauss", [&] {
    auto p     = good();
    p.quad     = rt3::quadrature_type::double_gauss;
    p.extra_mu = Vector{0.5};
    rt3::solve(p);
  });
  expect_throw("extra_mu = 0", "extra_mu values must be in (0, 1]", [&] {
    auto p     = good();
    p.extra_mu = Vector{0.0};
    rt3::solve(p);
  });
  expect_throw("nstokes * nmu_total = 4 * (16 + 1) > 64", "<= 64", [&] {
    auto p     = good();
    p.nmu      = 16;
    p.extra_mu = Vector{0.5};
    rt3::solve(p);
  });
  expect_throw("aziorder = -1", "aziorder must be >= 0", [&] {
    auto p     = good();
    p.aziorder = -1;
    rt3::solve(p);
  });
  expect_throw("2 aziorder + 1 = 513 > 512 with a beam", "2 * aziorder + 1 <= 512", [&] {
    auto p            = good();
    p.nstokes         = 1;
    p.nmu             = 1;
    p.aziorder        = 256;
    p.scattering_sets = {{.extinction = 1.0, .scattering = 0.9, .legendre = isotropic}};
    rt3::solve(p);
  });
  expect_throw("2 aziorder + 1 = 1025 > 1024 without a beam", "2 * aziorder + 1 <= 1024", [&] {
    auto p            = good();
    p.nstokes         = 1;
    p.nmu             = 1;
    p.aziorder        = 512;
    p.direct_flux     = 0.0;
    p.scattering_sets = {{.extinction = 1.0, .scattering = 0.9, .legendre = isotropic}};
    rt3::solve(p);
  });
  expect_throw("no layer", "at least 2 interfaces", [&] {
    auto p = good();
    layers(p, 0);
    rt3::solve(p);
  });
  expect_throw("201 layers", "at most 200 layers", [&] {
    auto p = good();
    layers(p, 201);
    rt3::solve(p);
  });
  expect_throw("(nlay + 1) * 64^2 > 101 * 4096", "(nlay + 1) * (nstokes", [&] {
    auto p = good();
    p.nmu  = 16;
    layers(p, 101);
    rt3::solve(p);
  });
  expect_throw("201 scattering sets", "at most 200 scattering sets", [&] {
    auto p            = good();
    p.scattering_sets = std::vector<rt3::scattering_set>(201, p.scattering_sets[0]);
    rt3::solve(p);
  });
  expect_throw(
      "scattering-matrix buffer (190 sets, N = 64, aziorder 16)", "scattering_sets.size() * (aziorder + 1)", [&] {
        auto p            = good();
        p.nmu             = 16;
        p.aziorder        = 16;
        p.direct_flux     = 0.0;
        p.scattering_sets = std::vector<rt3::scattering_set>(190, p.scattering_sets[0]);
        rt3::solve(p);
      });
  expect_throw("direct-beam buffer (100 layers, N = 64, aziorder 32)", "With a direct beam RT3 requires", [&] {
    auto p     = good();
    p.nmu      = 16;
    p.aziorder = 32;
    layers(p, 100);
    rt3::solve(p);
  });
  expect_throw("max_delta_tau = 0", "max_delta_tau must be positive", [&] {
    auto p          = good();
    p.max_delta_tau = 0.0;
    rt3::solve(p);
  });
  expect_throw("frequency = 0", "frequency must be positive", [&] {
    auto p      = good();
    p.frequency = 0.0;
    rt3::solve(p);
  });
  expect_throw("direct_flux < 0", "direct_flux must be >= 0", [&] {
    auto p        = good();
    p.direct_flux = -1.0;
    rt3::solve(p);
  });
  expect_throw("direct_mu = 0", "direct_mu must be in (0, 1]", [&] {
    auto p      = good();
    p.direct_mu = 0.0;
    rt3::solve(p);
  });
  expect_throw("direct beam over a Fresnel surface", "only over a Lambertian surface", [&] {
    auto p   = good();
    p.ground = rt3::fresnel_surface{.refractive_index = Complex{1.5, 0.0}};
    rt3::solve(p);
  });
  expect_throw("temperature with nlay values", "temperature needs nlay + 1", [&] {
    auto p        = good();
    p.temperature = Vector{250.0, 260.0};
    rt3::solve(p);
  });
  expect_throw("temperature 0 K with thermal", "temperature values must be positive", [&] {
    auto p           = good();
    p.temperature[0] = 0.0;
    rt3::solve(p);
  });
  expect_throw("negative gas extinction", "gas_extinction values must be non-negative", [&] {
    auto p              = good();
    p.gas_extinction[0] = -1e-3;
    rt3::solve(p);
  });
  expect_throw("layer_scattering_index out of range", "layer_scattering_index values must be <", [&] {
    auto p                      = good();
    p.layer_scattering_index[0] = 1;
    rt3::solve(p);
  });
  expect_throw("legendre with 5 columns", "must have legendre [nleg + 1, 6]", [&] {
    auto p                        = good();
    p.scattering_sets[0].legendre = Matrix(3, 5, 0.0);
    rt3::solve(p);
  });
  expect_throw("legendre of degree 1024", "at most 1024 Legendre coefficients", [&] {
    auto p                                 = good();
    p.scattering_sets[0].legendre          = henyey_greenstein(0.0, 1024);
    p.scattering_sets[0].legendre[1024, 0] = 1e-30;
    rt3::solve(p);
  });
  expect_throw("legendre[0, 0] = 0.99", "must be normalised", [&] {
    auto p                              = good();
    p.scattering_sets[0].legendre[0, 0] = 0.99;
    rt3::solve(p);
  });
  expect_throw("delta-M with f = 1", "divides by 1 - f", [&] {
    auto p                        = good();
    p.delta_m                     = true;
    p.scattering_sets[0].legendre = henyey_greenstein(1.0, 8);
    rt3::solve(p);
  });
  expect_throw("Mie series of degree 11 with gauss nmu = 2", "silently truncate the Legendre series", [&] {
    auto p                        = good();
    p.nmu                         = 2;
    p.scattering_sets[0].legendre = mie_legendre();
    rt3::solve(p);
  });
  expect_throw("delta-M with double_gauss", "delta-M scaled Legendre series", [&] {
    auto p                        = good();
    p.quad                        = rt3::quadrature_type::double_gauss;
    p.delta_m                     = true;
    p.scattering_sets[0].legendre = henyey_greenstein(0.8, 40);
    rt3::solve(p);
  });
  expect_throw("degree 252 with aziorder > 0 (FFT1DR)", "summed to degree 252", [&] {
    auto p                        = good();
    p.nstokes                     = 1;
    p.nmu                         = 64;
    p.aziorder                    = 1;
    p.direct_flux                 = 0.0;
    p.scattering_sets[0].legendre = henyey_greenstein(0.5, 252);
    rt3::solve(p);
  });
  {
    // The boundary of the last limit runs (RT3 would STOP above it)
    auto p                        = good();
    p.nstokes                     = 1;
    p.nmu                         = 64;
    p.aziorder                    = 1;
    p.direct_flux                 = 0.0;
    p.scattering_sets[0].legendre = henyey_greenstein(0.5, 251);
    rt3::solve(p);
    std::cout << "    degree 251 with aziorder 1 and nmu 64 runs\n";
  }
}
}  // namespace

int main() try {
  if (not rt3::available()) throw std::runtime_error("rt3-test needs ENABLE_RT3=ON");
  test_quadrature();
  test_mietest();
  test_testa();
  test_gas_only();
  test_single_scattering();
  test_delta_m_identity();
  test_direct_beam();
  test_errors();
  std::cout << "All RT3 tests passed\n";
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::cerr << "rt3-test failed: " << e.what() << '\n';
  return EXIT_FAILURE;
}
