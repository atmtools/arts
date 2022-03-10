#include "absorptionlines.h"
#include "auto_md.h"
#include "gui/plot.h"
#include "linemixing.h"
#include "lineshape.h"
#include "messages.h"
#include "physics_funcs.h"
#include "predefined_absorption_models.h"

Vector VectorNLinSpaceConst(Numeric f0, Numeric f1, Index n) {
  Vector x;
  VectorNLinSpace(x, n, f0, f1, Verbosity{});
  return x;
}

Vector VectorNLogSpaceConst(Numeric f0, Numeric f1, Index n) {
  Vector x;
  VectorNLogSpace(x, n, f0, f1, Verbosity{});
  return x;
}
  
int main() {
  std::stringstream ss(
    "50.47422e9 1.154500755758e-22 3.9951775840448995e-20 73 75 5.16e-10 0.00153987482219061 -0.0540571637755708 7500.616827041697 0.754 0 0 6225.511966444609 0.754 0 0 36 37 37 37\n"
    "50.98776e9 3.066876845339999e-22 3.5819124286696463e-20 69 71 5.317e-10 0.00170486349206349 -0.0571515949553766 7950.6538366642 0.754 0 0 6675.548976067111 0.754 0 0 34 35 35 35\n"
    "51.50335e9 7.68368069854e-22 3.1909272496583036e-20 65 67 5.48e-10 0.00190028520499109 -0.0606208643829584 8400.690846286701 0.754 0 0 7200.5921539600295 0.754 0 0 32 33 33 33\n"
    "52.02142e9 1.8161427105639998e-21 2.822271112222729e-20 61 63 5.647e-10 0.00213411290322581 -0.0645374497084855 8850.727855909203 0.754 0 0 7650.629163582532 0.754 0 0 30 31 31 31\n"
    "52.54244e9 4.047198183e-21 2.476004007026811e-20 57 59 5.82e-10 0.00241711264367816 -0.0689937941931067 9465.778435726623 0.754 0 0 8280.680977054035 0.754 0 0 28 29 29 29\n"
    "53.06695e9 8.50211410888e-21 2.1521819518427926e-20 53 55 5.994e-10 0.00276410582010582 -0.0741096869563658 10500.863557858376 0.754 0 0 9300.764865531704 0.754 0 0 26 27 27 27\n"
    "53.59572e9 1.68033672709e-20 1.850857587485014e-20 49 51 6.173e-10 0.00319595076923077 -0.0800431763685832 10980.903034789046 0.754 0 0 9863.311127559833 0.754 0 0 24 25 25 25\n"
    "54.12997e9 3.1238374123599996e-20 1.572079581876167e-20 45 47 6.36e-10 0.00374275362318841 -0.0870071595435633 11400.937577103381 0.754 0 0 10365.852454971626 0.754 0 0 22 23 23 23\n"
    "54.67115e9 5.4562227356e-20 1.3158930273364604e-20 41 43 6.549e-10 0.00444926406926407 -0.0952954158056413 11558.450530471257 0.754 0 0 10650.87589439921 0.754 0 0 20 21 21 21\n"
    "55.22137e9 8.93681317298e-20 1.0823392419390375e-20 37 39 6.748e-10 0.00538433684210526 -0.105325006369713 11895.978287688133 0.754 0 0 10890.895632864545 0.754 0 0 18 19 19 19\n"
    "55.78382e9 1.3706511179760002e-19 8.714557695099746e-21 33 35 6.953e-10 0.00665837908496732 -0.117708429714795 12181.001727115718 0.754 0 0 11265.92647421663 0.754 0 0 16 17 17 17\n"
    "56.26479e9 8.172342405080002e-20 4.140349031321195e-23 5 3 5.844e-10 1.0011 0.970122237409002 16786.38045891932 0.754 0 0 17326.424870466322 0.754 0 0 2 1 1 1\n"
    "56.36339e9 1.962741222526e-19 6.832765782728652e-21 29 31 7.167e-10 0.00845753333333333 -0.133383069062231 12548.53195164076 0.754 0 0 11880.977054034049 0.754 0 0 14 15 15 15\n"
    "56.96818e9 2.615389403592e-19 5.1783166355965434e-21 25 27 7.394e-10 0.0111158241758242 -0.153858929609147 12698.544288181594 0.754 0 0 12278.509745867259 0.754 0 0 12 13 13 13\n"
    "57.61249e9 3.2257668480799997e-19 3.751464518768915e-21 21 23 7.637e-10 0.0152824242424242 -0.181732800635926 13073.57512953368 0.754 0 0 12668.541820873426 0.754 0 0 10 11 11 11\n"
    "58.32388e9 3.65447006302e-19 2.5524259548405997e-21 17 19 7.902e-10 0.0223600888888889 -0.221873181488793 13223.587466074512 0.754 0 0 13051.073279052554 0.754 0 0 8 9 9 9\n"
    "58.44659e9 2.251141567122e-19 3.2553079588197106e-22 9 7 7.682e-10 0.166946666666667 0.495924531107892 15083.740439180856 0.754 0 0 15061.238588699729 0.754 0 0 4 3 3 3\n"
    "59.16421e9 3.7444078004199996e-19 1.5813499272484117e-21 13 15 8.21e-10 0.0358654285714286 -0.284561797673629 13756.131260794473 0.754 0 0 13651.12262521589 0.754 0 0 6 7 7 7\n"
    "59.59098e9 3.32170043464e-19 8.387568847982832e-22 13 11 8.339e-10 0.0668482666666667 0.332459702435534 14176.165803108808 0.754 0 0 14251.171971379226 0.754 0 0 6 5 5 5\n"
    "60.30605e9 3.38165892624e-19 8.382821242463063e-22 9 11 8.59e-10 0.0668482666666667 -0.395594131107892 14401.184307920059 0.754 0 0 14281.174438687392 0.754 0 0 4 5 5 5\n"
    "60.43478e9 3.91828742606e-19 1.5805076742189629e-21 17 15 8.773e-10 0.0358654285714286 0.249910292599905 13621.120157907722 0.754 0 0 13516.11152232914 0.754 0 0 8 7 7 7\n"
    "61.15057e9 4.0202168617800005e-19 2.5505527364284057e-21 21 19 9.123e-10 0.0223600888888889 0.200164509726835 13253.58993338268 0.754 0 0 13006.069578090304 0.754 0 0 10 9 9 9\n"
    "61.80016e9 3.7144285546199997e-19 3.7486894539525465e-21 25 23 9.431e-10 0.0152824242424242 0.166924057814275 13103.577596841846 0.754 0 0 12676.04243770047 0.754 0 0 12 11 11 11\n"
    "62.41122e9 3.14182495984e-19 5.174709249979848e-21 29 27 9.728e-10 0.0111158241758242 0.143148250014611 12758.549222797928 0.754 0 0 12181.001727115718 0.754 0 0 14 13 13 13\n"
    "62.48626e9 2.508063703628e-19 3.228550533568117e-22 5 7 9.147e-10 0.166946666666667 -0.636228904075668 15503.774981495191 0.754 0 0 15241.25339254873 0.754 0 0 2 3 3 3\n"
    "62.998e9 2.454101061188e-19 6.82836977811976e-21 33 31 1.001e-09 0.00845753333333333 0.125300606185383 12511.02886750555 0.754 0 0 11693.461633358007 0.754 0 0 16 15 15 15\n"
    "63.56852e9 1.781066992978e-19 8.709400881740195e-21 37 35 1.029e-09 0.00665837908496732 0.111410713972052 12143.498642980509 0.754 0 0 11280.927707870713 0.754 0 0 18 17 17 17\n"
    "64.12778e9 1.206364850992e-19 1.081749068884676e-20 41 39 1.057e-09 0.00538433684210526 0.100294006281832 11828.472736244757 0.754 0 0 10658.376511226252 0.754 0 0 20 19 19 19\n"
    "64.6789e9 7.641709754420001e-20 1.3152299517203517e-20 45 43 1.085e-09 0.00444926406926407 0.091195618041587 11565.951147298298 0.754 0 0 10650.87589439921 0.754 0 0 22 21 21 21\n"
    "65.22412e9 4.53885781412e-20 1.5713445969212232e-20 49 47 1.112e-09 0.00374275362318841 0.0836117897019165 11130.915371329878 0.754 0 0 10125.832716506291 0.754 0 0 24 23 23 23\n"
    "65.76474e9 2.532646685184e-20 1.8500512891249828e-20 53 51 1.141e-09 0.00319595076923077 0.0771934932241721 10950.900567480878 0.754 0 0 9825.808043424624 0.754 0 0 26 25 25 25\n"
    "66.30207e9 1.328380381398e-20 2.151305134656003e-20 57 55 1.169e-09 0.00276410582010582 0.071691272025619 10493.362941031335 0.754 0 0 9300.764865531704 0.754 0 0 28 27 27 27\n"
    "66.83678e9 6.5594589810399995e-21 2.4750568696578466e-20 61 59 1.198e-09 0.00241711264367816 0.0669221034719264 9465.778435726623 0.754 0 0 8280.680977054035 0.754 0 0 30 29 29 29\n"
    "67.36954e9 3.0488892978600002e-21 2.8212542506053354e-20 65 63 1.227e-09 0.00213411290322581 0.0627486674132615 8850.727855909203 0.754 0 0 7650.629163582532 0.754 0 0 32 31 31 31\n"
    "67.90087e9 1.3370743626800001e-21 3.1898406637924825e-20 69 67 1.256e-09 0.00190028520499109 0.0590659210057967 8400.690846286701 0.754 0 0 7200.5921539600295 0.754 0 0 34 33 33 33\n"
    "68.431e9 5.528172925519999e-22 3.580756714489144e-20 73 71 1.286e-09 0.00170486349206349 0.0557921427545498 7950.6538366642 0.754 0 0 6675.548976067111 0.754 0 0 36 35 35 35\n"
    "68.96031e9 2.158205905142e-22 3.993938637784364e-20 77 75 1.316e-09 0.00153987482219061 0.0528628014412826 7500.616827041697 0.754 0 0 6225.511966444609 0.754 0 0 38 37 37 37\n"
    "118.75034e9 2.9847337118480003e-19 0.0 1 3 4.48e-09 1.0011 0 16726.375524302985 0.754 0 0 16838.88477670861 0.754 0 0 0 1 1 1\n");
  
  const SpeciesTag sp("O2-66");
  const QuantumIdentifier qid("O2-66 S 1 1 Lambda 0 0 ElecStateLabel X X v1 0 0");
  const LineShape::Model metamodel = LineShape::MetaData2ModelShape("G0 T1 T1");
  AbsorptionLines band(false, false, 38, Absorption::CutoffType::None,
                       Absorption::MirroringType::None, Absorption::PopulationType::ByMakarovFullRelmat,
                       Absorption::NormalizationType::None, LineShape::Type::VP, 296, 750e9, -1, qid, 
                       {Species::Species::Oxygen, Species::Species::Nitrogen}, QuantumNumberLocalState{"J 0 0", "N 0 0"}, metamodel);
  ss >> band;
  Index wigner_initialized;
  Wigner6Init(wigner_initialized, 20000000, 250, Verbosity{});
  
  // Initializing values
  constexpr Index nfreq = 100000;
  Vector f_grid = VectorNLinSpaceConst(40e9, 130e9, nfreq);
  const Numeric P=Conversion::torr2pa(755.0);
  const Numeric T=296;
  const Numeric H=50e-6;
  const Vector VMR = {0.2098, 1 - 0.2098};
  PropagationMatrix mpm_abs(nfreq, 1);
  ArrayOfPropagationMatrix dmpm_abs(0);
  ArrayOfMatrix dxsec, dsource, dphase;
  ArrayOfArrayOfSpeciesTag specs(2, ArrayOfSpeciesTag(1));
  specs[0][0] = SpeciesTag("O2");
  specs[1][0] = SpeciesTag("N2");
  Matrix VMRmat(2, 1); VMRmat(0, 0) = VMR[0]; VMRmat(1, 0) = VMR[1];
  Matrix xsec2(nfreq, 1, 0), xsec3(nfreq, 1, 0), source, phase;
  
  // Retrieval Quantities
  RetrievalQuantity rq;
  rq.Target(Jacobian::Target(Jacobian::Atm::Temperature));
  rq.Target().perturbation = 0.1;
  
  Absorption::LineMixing::ErrorCorrectedSuddenData ecs_data{QuantumIdentifier("O2-66")};
  ecs_data[Species::Species::Oxygen].scaling = LineShapeModelParameters(LineShapeTemperatureModel::T0, 1.0, 0, 0, 0);
  ecs_data[Species::Species::Oxygen].collisional_distance = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::angstrom2meter(0.61), 0, 0, 0);
  ecs_data[Species::Species::Oxygen].lambda = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.39, 0, 0, 0);
  ecs_data[Species::Species::Oxygen].beta = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.567, 0, 0, 0);
  ecs_data[Species::Species::Oxygen].mass = 31.989830;
  ecs_data[Species::Species::Nitrogen].scaling = LineShapeModelParameters(LineShapeTemperatureModel::T0, 1.0, 0, 0, 0);
  ecs_data[Species::Species::Nitrogen].collisional_distance = LineShapeModelParameters(LineShapeTemperatureModel::T0, Conversion::angstrom2meter(0.61), 0, 0, 0);
  ecs_data[Species::Species::Nitrogen].lambda = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.39, 0, 0, 0);
  ecs_data[Species::Species::Nitrogen].beta = LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.567, 0, 0, 0);
  ecs_data[Species::Species::Nitrogen].mass = 28.006148;
  
  
  // Line Mixing full calculations
  const auto [abs, dabs, error] = Absorption::LineMixing::ecs_absorption(T, 0, P, 1, VMR, ecs_data, f_grid,  Zeeman::Polarization::None, band, {rq});
  
  // Line Mixing full calculations
  ComplexVector absdT = Absorption::LineMixing::ecs_absorption(T+0.1, 0, P, 1, VMR, ecs_data, f_grid, Zeeman::Polarization::None, band).abs;
  
  // Line Mixing full calculations with Zeeman (ignoring polarization...)
  auto [absZ, dabsZ, errorZ] = Absorption::LineMixing::ecs_absorption(T, H, P, 1, VMR, ecs_data, f_grid,  Zeeman::Polarization::Pi, band);
  absZ += Absorption::LineMixing::ecs_absorption(T, H, P, 1, VMR, ecs_data, f_grid, 
                                                        Zeeman::Polarization::SigmaMinus, band).abs;
  absZ += Absorption::LineMixing::ecs_absorption(T, H, P, 1, VMR, ecs_data, f_grid, 
                                                        Zeeman::Polarization::SigmaPlus, band).abs;
  
  // Line Mixing reimplementation of MPM19
  Absorption::PredefinedModel::VMRS vmrs_predef;
  vmrs_predef.O2 = 1;
  Absorption::PredefinedModel::compute(mpm_abs, dmpm_abs, Species::Isotopologues[Species::find_species_index("O2", "MPM2020")], f_grid, P, T, vmrs_predef, ArrayOfRetrievalQuantity(0));
  
  // Line by line calculations
  band.normalization = Absorption::NormalizationType::SFS;
  band.population = Absorption::PopulationType::LTE;
  LineShape::ComputeData com_lte(f_grid, {rq}, false);
  LineShape::ComputeData sparse_com_lte(Vector(0), {rq}, false);
  LineShape::compute(com_lte, sparse_com_lte, band, {rq}, {},
                     band.BroadeningSpeciesVMR(VMR, specs), {}, 1.0, 1.0, P, T, 0, 0,
                     Zeeman::Polarization::None, Options::LblSpeedup::None, false);
  
  band.population = Absorption::PopulationType::ByMakarovFullRelmat;
//   auto data = Absorption::LineMixing::ecs_eigenvalue_adaptation_test(band,
//     VectorNLinSpaceConst(200, 350, 76), ecs_data,
//     VectorNLogSpaceConst(1, 1'000'000'000'000, 101));
//   WriteXML("ascii", data, "prestemp.xml", 0, "", "", "", Verbosity());
  Verbosity v;
  Absorption::LineMixing::ecs_eigenvalue_adaptation(band,
                                                    VectorNLinSpaceConst(200, 350, 76),
                                                    ecs_data,
                                                    Conversion::atm2pa(1), 2, false, false, v);
  LineShape::ComputeData com(f_grid, {rq}, false);
  LineShape::ComputeData sparse_com(Vector(0), {rq}, false);
  LineShape::compute(com, sparse_com, band, {rq}, {},
                     band.BroadeningSpeciesVMR(VMR, specs), {}, 1.0, 1.0, P, T, 0, 0,
                     Zeeman::Polarization::None, Options::LblSpeedup::None, false);
  
  // Plot it all
  f_grid /= 1e9;  // Rescale for easier axis
  ARTSGUI::PlotConfig::Legend = {"Full", "MPM2020", "LBL", "FullZeeman", "Adapted"};
  ARTSGUI::PlotConfig::X = "Frequency [GHz]";
  ARTSGUI::PlotConfig::Y = "Absorption [1/m]";
  ARTSGUI::plot(f_grid, abs.real(), f_grid, mpm_abs.Kjj(), f_grid, com_lte.F.real(), f_grid, absZ.real(), f_grid, com.F.real());
}
