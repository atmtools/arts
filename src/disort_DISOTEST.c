/* DISOTEST.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Common Block Declarations */

struct dochek_1_ {
    real tstfir[360]	/* was [5][8][9] */, tstfdn[360]	/* was [5][8][
	    9] */, tstfup[360]	/* was [5][8][9] */, tstdfd[360]	/* 
	    was [5][8][9] */, tstuu[10800]	/* was [5][10][3][8][9] */;
};

#define dochek_1 (*(struct dochek_1_ *) &dochek_)

/* Initialized data */

struct {
    real e_1[2];
    integer fill_2[3];
    real e_3[2];
    integer fill_4[3];
    real e_5[2];
    integer fill_6[3];
    real e_7[2];
    integer fill_8[3];
    real e_9[2];
    integer fill_10[3];
    real e_11[2];
    integer fill_12[13];
    real e_13[2];
    integer fill_14[3];
    real e_15[2];
    integer fill_16[3];
    real e_17[2];
    integer fill_18[3];
    real e_19[2];
    integer fill_20[23];
    real e_21[2];
    integer fill_22[3];
    real e_23[2];
    integer fill_24[3];
    real e_25[2];
    integer fill_26[3];
    real e_27[2];
    integer fill_28[3];
    real e_29[2];
    integer fill_30[3];
    real e_31[2];
    integer fill_32[13];
    real e_33[3];
    integer fill_34[2];
    real e_35[3];
    integer fill_36[2];
    real e_37[3];
    integer fill_38[27];
    real e_39[3];
    integer fill_40[2];
    real e_41[3];
    integer fill_42[32];
    real e_43[2];
    integer fill_44[3];
    real e_45[3];
    integer fill_46[2];
    real e_47[3];
    integer fill_48[2];
    real e_49[3];
    integer fill_50[2];
    real e_51[3];
    integer fill_52[2];
    real e_53[3];
    integer fill_54[2];
    real e_55[3];
    integer fill_56[2];
    real e_57[3];
    integer fill_58[2];
    real e_59[3];
    integer fill_60[2];
    real e_61[3];
    integer fill_62[2];
    real e_63[3];
    integer fill_64[27];
    real e_65[3];
    integer fill_66[2];
    real e_67[3];
    integer fill_68[2];
    real e_69[3];
    integer fill_70[27];
    real e_71[15];
    integer fill_72[25];
    real e_73[2];
    integer fill_74[3];
    real e_75[2];
    integer fill_76[3];
    real e_77[2];
    integer fill_78[3];
    real e_79[2];
    integer fill_80[3];
    real e_81[2];
    integer fill_82[3];
    real e_83[2];
    integer fill_84[13];
    real e_85[2];
    integer fill_86[3];
    real e_87[2];
    integer fill_88[3];
    real e_89[2];
    integer fill_90[3];
    real e_91[2];
    integer fill_92[23];
    real e_93[2];
    integer fill_94[3];
    real e_95[2];
    integer fill_96[3];
    real e_97[2];
    integer fill_98[3];
    real e_99[2];
    integer fill_100[3];
    real e_101[2];
    integer fill_102[3];
    real e_103[2];
    integer fill_104[13];
    real e_105[3];
    integer fill_106[2];
    real e_107[3];
    integer fill_108[2];
    real e_109[3];
    integer fill_110[27];
    real e_111[3];
    integer fill_112[2];
    real e_113[3];
    integer fill_114[32];
    real e_115[2];
    integer fill_116[3];
    real e_117[3];
    integer fill_118[2];
    real e_119[3];
    integer fill_120[2];
    real e_121[3];
    integer fill_122[2];
    real e_123[3];
    integer fill_124[2];
    real e_125[3];
    integer fill_126[2];
    real e_127[3];
    integer fill_128[2];
    real e_129[3];
    integer fill_130[2];
    real e_131[3];
    integer fill_132[2];
    real e_133[3];
    integer fill_134[2];
    real e_135[3];
    integer fill_136[27];
    real e_137[3];
    integer fill_138[2];
    real e_139[3];
    integer fill_140[2];
    real e_141[3];
    integer fill_142[27];
    real e_143[15];
    integer fill_144[25];
    real e_145[2];
    integer fill_146[3];
    real e_147[2];
    integer fill_148[3];
    real e_149[2];
    integer fill_150[3];
    real e_151[2];
    integer fill_152[3];
    real e_153[2];
    integer fill_154[3];
    real e_155[2];
    integer fill_156[13];
    real e_157[2];
    integer fill_158[3];
    real e_159[2];
    integer fill_160[3];
    real e_161[2];
    integer fill_162[3];
    real e_163[2];
    integer fill_164[23];
    real e_165[2];
    integer fill_166[3];
    real e_167[2];
    integer fill_168[3];
    real e_169[2];
    integer fill_170[3];
    real e_171[2];
    integer fill_172[3];
    real e_173[2];
    integer fill_174[3];
    real e_175[2];
    integer fill_176[13];
    real e_177[3];
    integer fill_178[2];
    real e_179[3];
    integer fill_180[2];
    real e_181[3];
    integer fill_182[27];
    real e_183[3];
    integer fill_184[2];
    real e_185[3];
    integer fill_186[32];
    real e_187[2];
    integer fill_188[3];
    real e_189[3];
    integer fill_190[2];
    real e_191[3];
    integer fill_192[2];
    real e_193[3];
    integer fill_194[2];
    real e_195[3];
    integer fill_196[2];
    real e_197[3];
    integer fill_198[2];
    real e_199[3];
    integer fill_200[2];
    real e_201[3];
    integer fill_202[2];
    real e_203[3];
    integer fill_204[2];
    real e_205[3];
    integer fill_206[2];
    real e_207[3];
    integer fill_208[27];
    real e_209[3];
    integer fill_210[2];
    real e_211[3];
    integer fill_212[2];
    real e_213[3];
    integer fill_214[27];
    real e_215[15];
    integer fill_216[25];
    real e_217[2];
    integer fill_218[3];
    real e_219[2];
    integer fill_220[3];
    real e_221[2];
    integer fill_222[3];
    real e_223[2];
    integer fill_224[3];
    real e_225[2];
    integer fill_226[3];
    real e_227[2];
    integer fill_228[13];
    real e_229[2];
    integer fill_230[3];
    real e_231[2];
    integer fill_232[3];
    real e_233[2];
    integer fill_234[3];
    real e_235[2];
    integer fill_236[23];
    real e_237[2];
    integer fill_238[3];
    real e_239[2];
    integer fill_240[3];
    real e_241[2];
    integer fill_242[3];
    real e_243[2];
    integer fill_244[3];
    real e_245[2];
    integer fill_246[3];
    real e_247[2];
    integer fill_248[13];
    real e_249[3];
    integer fill_250[2];
    real e_251[3];
    integer fill_252[2];
    real e_253[3];
    integer fill_254[27];
    real e_255[3];
    integer fill_256[2];
    real e_257[3];
    integer fill_258[32];
    real e_259[2];
    integer fill_260[3];
    real e_261[3];
    integer fill_262[2];
    real e_263[3];
    integer fill_264[2];
    real e_265[3];
    integer fill_266[2];
    real e_267[3];
    integer fill_268[2];
    real e_269[3];
    integer fill_270[2];
    real e_271[3];
    integer fill_272[2];
    real e_273[3];
    integer fill_274[2];
    real e_275[3];
    integer fill_276[2];
    real e_277[3];
    integer fill_278[2];
    real e_279[3];
    integer fill_280[27];
    real e_281[3];
    integer fill_282[2];
    real e_283[3];
    integer fill_284[2];
    real e_285[3];
    integer fill_286[27];
    real e_287[15];
    integer fill_288[25];
    real e_289[2];
    integer fill_290[3];
    real e_291[2];
    integer fill_292[3];
    real e_293[2];
    integer fill_294[3];
    real e_295[2];
    integer fill_296[3];
    real e_297[2];
    integer fill_298[3];
    real e_299[2];
    integer fill_300[123];
    real e_301[2];
    integer fill_302[3];
    real e_303[2];
    integer fill_304[3];
    real e_305[2];
    integer fill_306[3];
    real e_307[2];
    integer fill_308[3];
    real e_309[2];
    integer fill_310[3];
    real e_311[2];
    integer fill_312[123];
    real e_313[2];
    integer fill_314[3];
    real e_315[2];
    integer fill_316[3];
    real e_317[2];
    integer fill_318[3];
    real e_319[2];
    integer fill_320[3];
    real e_321[2];
    integer fill_322[3];
    real e_323[2];
    integer fill_324[123];
    real e_325[2];
    integer fill_326[3];
    real e_327[2];
    integer fill_328[3];
    real e_329[2];
    integer fill_330[3];
    real e_331[2];
    integer fill_332[3];
    real e_333[2];
    integer fill_334[3];
    real e_335[2];
    integer fill_336[123];
    real e_337[2];
    integer fill_338[3];
    real e_339[2];
    integer fill_340[3];
    real e_341[2];
    integer fill_342[3];
    real e_343[2];
    integer fill_344[3];
    real e_345[2];
    integer fill_346[3];
    real e_347[2];
    integer fill_348[123];
    real e_349[2];
    integer fill_350[3];
    real e_351[2];
    integer fill_352[3];
    real e_353[2];
    integer fill_354[3];
    real e_355[2];
    integer fill_356[3];
    real e_357[2];
    integer fill_358[3];
    real e_359[2];
    integer fill_360[423];
    real e_361[2];
    integer fill_362[3];
    real e_363[2];
    integer fill_364[3];
    real e_365[2];
    integer fill_366[3];
    real e_367[2];
    integer fill_368[3];
    real e_369[2];
    integer fill_370[3];
    real e_371[2];
    integer fill_372[123];
    real e_373[2];
    integer fill_374[3];
    real e_375[2];
    integer fill_376[3];
    real e_377[2];
    integer fill_378[3];
    real e_379[2];
    integer fill_380[3];
    real e_381[2];
    integer fill_382[3];
    real e_383[2];
    integer fill_384[123];
    real e_385[2];
    integer fill_386[3];
    real e_387[2];
    integer fill_388[3];
    real e_389[2];
    integer fill_390[3];
    real e_391[2];
    integer fill_392[3];
    real e_393[2];
    integer fill_394[3];
    real e_395[2];
    integer fill_396[123];
    real e_397[2];
    integer fill_398[3];
    real e_399[2];
    integer fill_400[3];
    real e_401[2];
    integer fill_402[3];
    real e_403[2];
    integer fill_404[3];
    real e_405[2];
    integer fill_406[3];
    real e_407[2];
    integer fill_408[723];
    real e_409[2];
    integer fill_410[3];
    real e_411[2];
    integer fill_412[3];
    real e_413[2];
    integer fill_414[3];
    real e_415[2];
    integer fill_416[3];
    real e_417[2];
    integer fill_418[3];
    real e_419[2];
    integer fill_420[123];
    real e_421[2];
    integer fill_422[3];
    real e_423[2];
    integer fill_424[3];
    real e_425[2];
    integer fill_426[3];
    real e_427[2];
    integer fill_428[3];
    real e_429[2];
    integer fill_430[3];
    real e_431[2];
    integer fill_432[123];
    real e_433[2];
    integer fill_434[3];
    real e_435[2];
    integer fill_436[3];
    real e_437[2];
    integer fill_438[3];
    real e_439[2];
    integer fill_440[3];
    real e_441[2];
    integer fill_442[3];
    real e_443[2];
    integer fill_444[123];
    real e_445[2];
    integer fill_446[3];
    real e_447[2];
    integer fill_448[3];
    real e_449[2];
    integer fill_450[3];
    real e_451[2];
    integer fill_452[3];
    real e_453[2];
    integer fill_454[3];
    real e_455[2];
    integer fill_456[123];
    real e_457[2];
    integer fill_458[3];
    real e_459[2];
    integer fill_460[3];
    real e_461[2];
    integer fill_462[3];
    real e_463[2];
    integer fill_464[3];
    real e_465[2];
    integer fill_466[3];
    real e_467[2];
    integer fill_468[123];
    real e_469[2];
    integer fill_470[3];
    real e_471[2];
    integer fill_472[3];
    real e_473[2];
    integer fill_474[3];
    real e_475[2];
    integer fill_476[3];
    real e_477[2];
    integer fill_478[3];
    real e_479[2];
    integer fill_480[423];
    real e_481[3];
    integer fill_482[2];
    real e_483[3];
    integer fill_484[2];
    real e_485[3];
    integer fill_486[2];
    real e_487[3];
    integer fill_488[2];
    real e_489[3];
    integer fill_490[2];
    real e_491[3];
    integer fill_492[122];
    real e_493[3];
    integer fill_494[2];
    real e_495[3];
    integer fill_496[2];
    real e_497[3];
    integer fill_498[2];
    real e_499[3];
    integer fill_500[2];
    real e_501[3];
    integer fill_502[2];
    real e_503[3];
    integer fill_504[122];
    real e_505[3];
    integer fill_506[2];
    real e_507[3];
    integer fill_508[2];
    real e_509[3];
    integer fill_510[2];
    real e_511[3];
    integer fill_512[2];
    real e_513[3];
    integer fill_514[2];
    real e_515[3];
    integer fill_516[22];
    real e_517[3];
    integer fill_518[2];
    real e_519[3];
    integer fill_520[2];
    real e_521[3];
    integer fill_522[2];
    real e_523[3];
    integer fill_524[2];
    real e_525[3];
    integer fill_526[2];
    real e_527[3];
    integer fill_528[22];
    real e_529[3];
    integer fill_530[2];
    real e_531[3];
    integer fill_532[2];
    real e_533[3];
    integer fill_534[2];
    real e_535[3];
    integer fill_536[2];
    real e_537[3];
    integer fill_538[2];
    real e_539[3];
    integer fill_540[772];
    real e_541[3];
    integer fill_542[2];
    real e_543[3];
    integer fill_544[2];
    real e_545[3];
    integer fill_546[2];
    real e_547[3];
    integer fill_548[2];
    real e_549[3];
    integer fill_550[2];
    real e_551[3];
    integer fill_552[122];
    real e_553[3];
    integer fill_554[2];
    real e_555[3];
    integer fill_556[2];
    real e_557[3];
    integer fill_558[2];
    real e_559[3];
    integer fill_560[2];
    real e_561[3];
    integer fill_562[2];
    real e_563[3];
    integer fill_564[1022];
    real e_565[2];
    integer fill_566[3];
    real e_567[2];
    integer fill_568[3];
    real e_569[2];
    integer fill_570[3];
    real e_571[2];
    integer fill_572[133];
    real e_573[3];
    integer fill_574[2];
    real e_575[3];
    integer fill_576[2];
    real e_577[3];
    integer fill_578[2];
    real e_579[3];
    integer fill_580[132];
    real e_581[3];
    integer fill_582[2];
    real e_583[3];
    integer fill_584[2];
    real e_585[3];
    integer fill_586[2];
    real e_587[3];
    integer fill_588[132];
    real e_589[3];
    integer fill_590[2];
    real e_591[3];
    integer fill_592[2];
    real e_593[3];
    integer fill_594[2];
    real e_595[3];
    integer fill_596[132];
    real e_597[3];
    integer fill_598[2];
    real e_599[3];
    integer fill_600[2];
    real e_601[3];
    integer fill_602[2];
    real e_603[3];
    integer fill_604[132];
    real e_605[3];
    integer fill_606[2];
    real e_607[3];
    integer fill_608[2];
    real e_609[3];
    integer fill_610[2];
    real e_611[3];
    integer fill_612[132];
    real e_613[3];
    integer fill_614[2];
    real e_615[3];
    integer fill_616[2];
    real e_617[3];
    integer fill_618[2];
    real e_619[3];
    integer fill_620[132];
    real e_621[3];
    integer fill_622[2];
    real e_623[3];
    integer fill_624[2];
    real e_625[3];
    integer fill_626[2];
    real e_627[3];
    integer fill_628[132];
    real e_629[3];
    integer fill_630[2];
    real e_631[3];
    integer fill_632[2];
    real e_633[3];
    integer fill_634[2];
    real e_635[3];
    integer fill_636[32];
    real e_637[3];
    integer fill_638[2];
    real e_639[3];
    integer fill_640[2];
    real e_641[3];
    integer fill_642[2];
    real e_643[3];
    integer fill_644[82];
    real e_645[3];
    integer fill_646[2];
    real e_647[3];
    integer fill_648[2];
    real e_649[3];
    integer fill_650[2];
    real e_651[3];
    integer fill_652[32];
    real e_653[3];
    integer fill_654[2];
    real e_655[3];
    integer fill_656[2];
    real e_657[3];
    integer fill_658[2];
    real e_659[3];
    integer fill_660[82];
    real e_661[3];
    integer fill_662[2];
    real e_663[3];
    integer fill_664[2];
    real e_665[3];
    integer fill_666[2];
    real e_667[3];
    integer fill_668[32];
    real e_669[3];
    integer fill_670[2];
    real e_671[3];
    integer fill_672[2];
    real e_673[3];
    integer fill_674[2];
    real e_675[3];
    integer fill_676[832];
    real e_677[3];
    integer fill_678[2];
    real e_679[3];
    integer fill_680[2];
    real e_681[3];
    integer fill_682[2];
    real e_683[3];
    integer fill_684[132];
    real e_685[3];
    integer fill_686[2];
    real e_687[3];
    integer fill_688[2];
    real e_689[3];
    integer fill_690[2];
    real e_691[3];
    integer fill_692[132];
    real e_693[3];
    integer fill_694[2];
    real e_695[3];
    integer fill_696[2];
    real e_697[3];
    integer fill_698[2];
    real e_699[3];
    integer fill_700[882];
    real e_701[20];
    integer fill_702[130];
    real e_703[20];
    integer fill_704[130];
    real e_705[20];
    integer fill_706[30];
    real e_707[20];
    integer fill_708[30];
    real e_709[20];
    integer fill_710[780];
    } dochek_ = { (float)3.14159, (float)2.29844, {0}, (float)3.14159, (float)
	    2.29844, {0}, (float)0., (float)0., {0}, (float)3.14159, (float)
	    0., {0}, (float)3.14159, (float)0., {0}, (float)0., (float)0., {0}
	    , (float).252716, (float).0210311, {0}, (float).252716, (float)
	    .0210311, {0}, (float).252716, (float)2.56077e-28, {0}, (float)
	    .252716, (float)2.56077e-28, {0}, (float)3.14159, (float)
	    1.42628e-4, {0}, (float)3.14159, (float)1.42628e-4, {0}, (float)
	    0., (float)0., {0}, (float)3.14159, (float)5.67011e-35, {0}, (
	    float)3.14159, (float)5.67011e-35, {0}, (float)0., (float)0., {0},
	     (float)3.14159, (float)1.90547, (float)1.15573, {0}, (float)
	    3.14159, (float)1.90547, (float)1.15573, {0}, (float)1.5708, (
	    float).577864, (float).212584, {0}, (float)3.14159, (float)
	    3.97856e-14, (float)5.03852e-28, {0}, (float).128058, (float)
	    8.67322e-6, (float)4.47729e-21, {0}, (float)100., (float)100., {0}
	    , (float)100., (float)36.7879, (float)13.5335, {0}, (float)100., (
	    float)36.7879, (float)13.5335, {0}, (float)100., (float)36.7879, (
	    float)13.5335, {0}, (float)100., (float)36.7879, (float)13.5335, {
	    0}, (float)100., (float)36.7879, (float)13.5335, {0}, (float)100.,
	     (float)36.7879, (float)13.5335, {0}, (float)100., (float)13.5335,
	     (float)2.06115e-7, {0}, (float)100., (float)36.7879, (float)
	    13.5335, {0}, (float)100., (float)36.7879, (float)13.5335, {0}, (
	    float)100., (float)36.7879, (float)13.5335, {0}, (float)0., (
	    float)0., (float)0., {0}, (float)0., (float)0., (float)0., {0}, (
	    float)0., (float)0., (float)0., {0}, (float)0., (float)0., (float)
	    0., (float)0., (float)0., (float)0., (float)0., (float)0., (float)
	    0., (float)0., (float)1.5708, (float).192354, (float).023555, (
	    float)9.65131e-6, (float)9.03133e-19, {0}, (float)0., (float)
	    .0794108, {0}, (float)0., (float).420233, {0}, (float)3.14159, (
	    float)3.04897, {0}, (float)0., (float)0., {0}, (float)0., (float)
	    .0676954, {0}, (float)3.14159, (float).00460048, {0}, (float)0., (
	    float).0441791, {0}, (float)0., (float).106123, {0}, (float)0., (
	    float)2.51683e-4, {0}, (float)0., (float).0268008, {0}, (float)0.,
	     (float).0470303, {0}, (float)0., (float)1.31455, {0}, (float)
	    3.14159, (float)2.4966, {0}, (float)0., (float)1.99813e-5, {0}, (
	    float)0., (float).591973, {0}, (float)3.14159, (float)1.01167, {0}
	    , (float)0., (float)1.17401, (float)1.81264, {0}, (float)0., (
	    float)1.01517, (float)1.51554, {0}, (float)0., (float).702764, (
	    float).803294, {0}, (float)0., (float)2.24768, (float).479851, {0}
	    , (float)1.74767, (float).233975, (float)6.38347e-5, {0}, (float)
	    0., (float)0., {0}, (float)0., (float)0., (float)0., {0}, (float)
	    0., (float)0., (float)0., {0}, (float)0., (float)0., (float)0., {
	    0}, (float)0., (float)0., (float)0., {0}, (float)321.497, (float)
	    142.493, (float)70.5306, {0}, (float)321.497, (float)304.775, (
	    float)363.632, {0}, (float)321.497, (float)255.455, (float)
	    443.444, {0}, (float)319.83, (float)354.099, (float)301.334, {0}, 
	    (float)319.83, (float)350.555, (float)292.063, {0}, (float)319.83,
	     (float)353.25, (float)298.583, {0}, (float)1., (float).722235, (
	    float).513132, {0}, (float)1., (float).795332, (float).650417, {0}
	    , (float)1., (float).486157, (float).159984, {0}, (float)1., (
	    float).355151, (float).144265, (float).00671445, (float)
	    6.16968e-7, (float)1., (float).452357, (float).236473, (float)
	    .0276475, (float)7.41854e-5, (float)6.09217, (float)4.97279, (
	    float)4.46616, (float)4.22731, (float)4.73767, {0}, (float)
	    .0799451, (float)0., {0}, (float).422922, (float)0., {0}, (float)
	    .0906556, (float)0., {0}, (float).259686, (float)0., {0}, (float)
	    3.0739, (float)0., {0}, (float)2.49618, (float)0., {0}, (float)
	    .0535063, (float)0., {0}, (float).125561, (float)0., {0}, (float)
	    .062473, (float)0., {0}, (float).225915, (float)0., {0}, (float)
	    .169951, (float)0., {0}, (float)1.8269, (float)0., {0}, (float)
	    .583122, (float)0., {0}, (float).170062, (float)0., {0}, (float)
	    2.54962, (float)0., {0}, (float)1.6875, (float)0., {0}, (float)
	    .173223, (float).111113, (float)0., {0}, (float).123666, (float)
	    .0788691, (float)0., {0}, (float).225487, (float).123848, (float)
	    0., {0}, (float)2.66174, (float)1.76783, (float)0., {0}, (float)
	    .270485, (float).0374253, (float)1.02904e-5, {0}, (float)0., (
	    float)0., {0}, (float)0., (float)0., (float)0., {0}, (float)
	    1.4845, (float)2.99914, (float)6.76676, {0}, (float).619105, (
	    float)1.33245, (float)3.4779, {0}, (float)82.8836, (float)165.215,
	     (float)359.895, {0}, (float)85.2994, (float)170.343, (float)
	    372.842, {0}, (float)336.115, (float)413.977, (float)443.05, {0}, 
	    (float)237.35, (float)261.13, (float)456.205, {0}, (float)429.572,
	     (float)447.018, (float)594.576, {0}, (float)312.563, (float)
	    268.126, (float)305.596, {0}, (float)407.04, (float)411.058, (
	    float)530.504, {0}, (float).0929634, (float).0278952, (float)0., {
	    0}, (float).225136, (float).126349, (float)0., {0}, (float)
	    .378578, (float).243397, (float)0., {0}, (float).227973, (float)
	    .0875098, (float).0361819, (float).00219291, (float)0., (float)
	    .100079, (float).0452015, (float).0241941, (float).00416017, (
	    float)0., (float)4.68414, (float)4.24381, (float)4.16941, (float)
	    4.30667, (float)5.11524, {0}, (float)25.4067, (float)18.6531, {0},
	     (float)0., (float)0., {0}, (float).066687, (float).0588936, {0}, 
	    (float)25.7766, (float)0., {0}, (float)0., (float)0., {0}, (float)
	    .114239, (float)7.93633e-5, {0}, (float)1.6657, (float).189848, {
	    0}, (float)0., (float)0., {0}, (float)1.67462, (float)1.75464e-4, 
	    {0}, (float)0., (float)0., {0}, (float)25.9221, (float).0743638, {
	    0}, (float)0., (float)0., {0}, (float).0806945, (float).0435243, {
	    0}, (float)25.9222, (float)1.91941e-5, {0}, (float)0., (float)0., 
	    {0}, (float).100459, (float).0171385, {0}, (float)0., (float)0., (
	    float)0., {0}, (float).343724, (float).35239, (float).31945, {0}, 
	    (float).385003, (float).337317, (float).216403, {0}, (float)0., (
	    float)0., (float)0., {0}, (float).310129, (float).0452671, (float)
	    1.25022e-5, {0}, (float)200., (float)200., {0}, (float)200., (
	    float)73.5759, (float)27.0671, {0}, (float)202.01, (float)77.9962,
	     (float)40.6006, {0}, (float)200.899, (float)75.7612, (float)
	    36.5186, {0}, (float)310.552, (float)311.018, (float)679.533, {0},
	     (float)957.001, (float)529.253, (float)807.923, {0}, (float)
	    581.342, (float)128.621, (float)-171.119, {0}, (float)423.78, (
	    float)61.9828, (float)-31.7719, {0}, (float)-80.427, (float)
	    251.589, (float)715.964, {0}, (float)-168.356, (float)101.251, (
	    float)409.326, {0}, (float)-98.5693, (float)217.724, (float)
	    623.936, {0}, (float)1.12474, (float).651821, (float).563361, {0},
	     (float).512692, (float).356655, (float).0568095, {0}, (float)
	    .565095, (float).276697, (float).0135679, {0}, (float).882116, (
	    float).232366, (float).0933443, (float).00392782, (float)1.025e-7,
	     (float).804577, (float).25533, (float).130976, (float).0136227, (
	    float)1.22022e-5, (float)3.49563, (float).881206, (float).350053, 
	    (float).0193471, (float).0715349, {0}, (float)0., (float).0133826,
	     {0}, (float)0., (float).0263324, {0}, (float)0., (float).115898, 
	    {0}, (float).117771, (float)0., {0}, (float).026417, (float)0., {
	    0}, (float).0134041, (float)0., {0}, (float)0., (float).0708109, {
	    0}, (float)0., (float).139337, {0}, (float)0., (float).613458, {0}
	    , (float).622884, (float)0., {0}, (float).139763, (float)0., {0}, 
	    (float).0709192, (float)0., {0}, (float)1., (float).984447, {0}, (
	    float)1., (float).969363, {0}, (float)1., (float).863946, {0}, (
	    float).133177, (float)0., {0}, (float).0299879, (float)0., {0}, (
	    float).0152233, (float)0., {0}, (float)0., (float)1.2298e-15, {0},
	     (float)0., (float)1.30698e-17, {0}, (float)0., (float)
	    6.88841e-18, {0}, (float).262972, (float)0., {0}, (float).0906967,
	     (float)0., {0}, (float).0502853, (float)0., {0}, (float)0., (
	    float).0271316, {0}, (float)0., (float).0187805, {0}, (float)0., (
	    float).0116385, {0}, (float)1.93321, (float)0., {0}, (float)
	    1.02732, (float)0., {0}, (float).797199, (float)0., {0}, (float)
	    1., (float).0018684, {0}, (float)1., (float).00126492, {0}, (
	    float)1., (float)7.79281e-4, {0}, (float).87751, (float)0., {0}, (
	    float).815136, (float)0., {0}, (float).752715, (float)0., {0}, (
	    float)0., (float).00771897, {0}, (float)0., (float).0200778, {0}, 
	    (float)0., (float).0257685, {0}, (float).161796, (float)0., {0}, (
	    float).0211501, (float)0., {0}, (float).00786713, (float)0., {0}, 
	    (float)0., (float).0186027, {0}, (float)0., (float).0464061, {0}, 
	    (float)0., (float).0677603, {0}, (float).347678, (float)0., {0}, (
	    float).048712, (float)0., {0}, (float).0189387, (float)0., {0}, (
	    float)0., (float)1.70004e-4, {0}, (float)0., (float)3.97168e-5, {
	    0}, (float)0., (float)1.32472e-5, {0}, (float).162566, (float)0., 
	    {0}, (float).0245786, (float)0., {0}, (float).0101498, (float)0., 
	    {0}, (float)0., (float).010595, {0}, (float)0., (float).00769337, 
	    {0}, (float)0., (float).00379276, {0}, (float).36401, (float)0., {
	    0}, (float).0826993, (float)0., {0}, (float).049237, (float)0., {
	    0}, (float)0., (float).00788518, {0}, (float)0., (float).0222081, 
	    {0}, (float)0., (float).00290003, {0}, (float).486923, (float)0., 
	    {0}, (float).0508802, (float)0., {0}, (float).0103351, (float)0., 
	    {0}, (float)0., (float).213998, {0}, (float)0., (float).58401, {0}
	    , (float)0., (float).377162, {0}, (float)3.72886, (float)0., {0}, 
	    (float).632069, (float)0., {0}, (float).151211, (float)0., {0}, (
	    float)1., (float).911175, {0}, (float)1., (float).744031, {0}, (
	    float)1., (float).403607, {0}, (float).566618, (float)0., {0}, (
	    float).232213, (float)0., {0}, (float).0760274, (float)0., {0}, (
	    float)0., (float)2.07944e-5, {0}, (float)0., (float)9.76964e-7, {
	    0}, (float)0., (float)1.88004e-7, {0}, (float).486926, (float)0., 
	    {0}, (float).050919, (float)0., {0}, (float).0103691, (float)0., {
	    0}, (float)0., (float).236362, {0}, (float)0., (float).164357, {0}
	    , (float)0., (float).0915016, {0}, (float)3.85776, (float)0., {0},
	     (float).868015, (float)0., {0}, (float).379793, (float)0., {0}, (
	    float)1., (float).416334, {0}, (float)1., (float).275956, {0}, (
	    float)1., (float).151785, {0}, (float).757316, (float)0., {0}, (
	    float).587262, (float)0., {0}, (float).433051, (float)0., {0}, (
	    float)0., (float)2.49435, (float)3.34666, {0}, (float)0., (float)
	    .11907, (float).219624, {0}, (float)0., (float).135138, (float)
	    .157002, {0}, (float).0929865, (float).12407, (float)0., {0}, (
	    float).0661234, (float).0402884, (float)0., {0}, (float).0355365, 
	    (float).0173581, (float)0., {0}, (float)0., (float)2.22302, (
	    float)2.94684, {0}, (float)0., (float).0964101, (float).167508, {
	    0}, (float)0., (float).0962927, (float).108212, {0}, (float)
	    .0655782, (float).0844926, (float)0., {0}, (float).0456643, (
	    float).0280216, (float)0., {0}, (float).0274243, (float).0135087, 
	    (float)0., {0}, (float)0., (float).0476042, (float).0837538, {0}, 
	    (float)0., (float)3.00258, (float)2.68792, {0}, (float)0., (float)
	    1.4114, (float).876317, {0}, (float).870016, (float).6974, (float)
	    0., {0}, (float).22476, (float).109065, (float)0., {0}, (float)
	    .0228322, (float).00935116, (float)0., {0}, (float)0., (float)
	    .0476042, (float).0837538, {0}, (float)0., (float).0582549, (
	    float).0943348, {0}, (float)0., (float).104636, (float).0895961, {
	    0}, (float).0886109, (float).0915334, (float)0., {0}, (float)
	    .0576913, (float).0295681, (float)0., {0}, (float).0228322, (
	    float).00935116, (float)0., {0}, (float)0., (float).0476042, (
	    float).0837538, {0}, (float)0., (float).0259416, (float).0400024, 
	    {0}, (float)0., (float).0624961, (float).0466785, {0}, (float)
	    .0695291, (float).0590188, (float)0., {0}, (float).0493283, (
	    float).0244592, (float)0., {0}, (float).0228322, (float).00935116,
	     (float)0., {0}, (float)0., (float).753662, (float).19523, {0}, (
	    float)0., (float).696362, (float).13199, {0}, (float)0., (float)
	    .650541, (float).0720655, {0}, (float).45766, (float).627631, (
	    float)0., {0}, (float).768942, (float).581809, (float)0., {0}, (
	    float)1.03122, (float).524532, (float)0., {0}, (float)15.0205, (
	    float).18509, (float)3.41287e-5, {0}, (float).218778, (float)
	    .0492726, (float)1.39916e-5, {0}, (float).13652, (float).0265448, 
	    (float)7.47042e-6, {0}, (float).114001, (float).0202154, (float)
	    5.65604e-6, {0}, (float).0871231, (float).012966, (float)
	    3.58246e-6, {0}, (float).0855016, (float).00951226, (float)
	    2.57859e-6, {0}, (float)0., (float)0., {0}, (float)0., (float)0., 
	    {0}, (float)0., (float)0., {0}, (float)0., (float)0., {0}, (float)
	    0., (float)0., (float)0., {0}, (float)0., (float)0., (float)0., {
	    0}, (float)0., (float)0., (float)0., {0}, (float)0., (float)0., (
	    float)0., {0}, (float)0., (float)0., (float)0., {0}, (float)0., (
	    float)0., (float)0., {0}, (float)9.77882e-5, (float).0145131, (
	    float)2.15393, {0}, (float).792386, (float)1.30642, (float)
	    2.15393, {0}, (float)0., (float)0., (float)0., {0}, (float)0., (
	    float)0., (float)0., {0}, (float)5.07347e-5, (float).0075297, (
	    float)1.11751, {0}, (float).272396, (float).449105, (float)
	    .740449, {0}, (float)0., (float)0., (float)0., {0}, (float)0., (
	    float)0., (float)0., {0}, (float).00320501, (float).475666, (
	    float)70.5951, {0}, (float)46.4836, (float)76.6386, (float)
	    126.356, {0}, (float)102.336, (float)62.0697, (float)37.6472, {0},
	     (float)102.336, (float).689532, (float).00464603, {0}, (float)
	    .00358884, (float).532631, (float)79.0494, {0}, (float)47.5119, (
	    float)78.3338, (float)129.151, {0}, (float)102.336, (float)
	    97.8748, (float)110.061, {0}, (float)102.336, (float)101.048, (
	    float)138.631, {0}, (float)78.0731, (float)115.783, (float)
	    133.373, {0}, (float)117.126, (float)136.113, (float)142.865, {0},
	     (float)102.336, (float)84.9992, (float)138.631, {0}, (float)
	    102.336, (float)77.3186, (float)145.441, {0}, (float)71.2616, (
	    float)78.831, (float)144.088, {0}, (float)78.0737, (float)85.6424,
	     (float)145.547, {0}, (float)101.805, (float)106.346, (float)
	    96.416, {0}, (float)101.805, (float)134.218, (float)88.7832, {0}, 
	    (float)141.637, (float)102.514, (float)189.259, {0}, (float)
	    148.416, (float)158.864, (float)189.259, {0}, (float)101.805, (
	    float)106.346, (float)96.416, {0}, (float)101.805, (float)106.238,
	     (float)74.8649, {0}, (float)127.778, (float)94.1347, (float)
	    189.259, {0}, (float)148.416, (float)158.864, (float)189.259, {0},
	     (float)101.805, (float)105.967, (float)95.3651, {0}, (float)
	    101.805, (float)128.779, (float)74.7553, {0}, (float)135.84, (
	    float)90.0044, (float)97.2743, {0}, (float)95.6588, (float)
	    88.7624, (float)97.2743, {0}, (float)101.805, (float)105.967, (
	    float)95.3651, {0}, (float)101.805, (float)100.799, (float)60.837,
	     {0}, (float)121.98, (float)81.6247, (float)97.2743, {0}, (float)
	    95.6588, (float)88.7624, (float)97.2743, {0}, (float)101.805, (
	    float)106.259, (float)96.1392, {0}, (float)101.805, (float)
	    132.945, (float)84.5884, {0}, (float)140.347, (float)99.2793, (
	    float)154.901, {0}, (float)140.534, (float)148.621, (float)
	    176.268, {0}, (float)101.805, (float)106.259, (float)96.1392, {0},
	     (float)101.805, (float)104.823, (float)69.7291, {0}, (float)
	    126.343, (float)90.2686, (float)137.958, {0}, (float)140.534, (
	    float)148.621, (float)176.268, {0}, (float).31831, (float).262711,
	     (float).210014, {0}, (float).31831, (float).136952, (float)
	    .0560376, {0}, (float).0562566, (float).0184909, (float)0., {0}, (
	    float).0194423, (float).00552188, (float)0., {0}, (float).31831, (
	    float).277499, (float).240731, {0}, (float).31831, (float).18395, 
	    (float).129291, {0}, (float).123687, (float).0835695, (float)0., {
	    0}, (float).0495581, (float).0250575, (float)0., {0}, (float)
	    .31831, (float).18902, (float).0684762, {0}, (float).31831, (
	    float).0988158, (float).0296698, {0}, (float).149335, (float)
	    .0965192, (float)0., {0}, (float).104766, (float).0654445, (float)
	    0., {0}, (float).31831, (float).153507, (float).0706614, (float)
	    .00372784, (float)2.87656e-7, (float).31831, (float).0509531, (
	    float).0209119, (float).00108815, (float)1.05921e-7, (float)
	    .0998915, (float).0367006, (float).0148545, (float)8.83317e-4, (
	    float)0., (float).0591345, (float).0231903, (float).00972307, (
	    float)5.94743e-4, (float)0., {0}, (float).31831, (float).196609, (
	    float).115478, (float).0146177, (float)3.37743e-5, (float).31831, 
	    (float).0592369, (float).0301809, (float).0038559, (float)
	    1.20858e-5, (float).0739198, (float).030023, (float).0152672, (
	    float).00238301, (float)0., (float).0132768, (float).00705566, (
	    float).00406932, (float)7.77891e-4, (float)0., {0}, (float)1.9392,
	     (float)1.66764, (float)1.48511, (float)1.34514, (float)1.48927, (
	    float)1.9392, (float)1.44453, (float)1.35009, (float)1.35131, (
	    float)1.5427, (float)1.61855, (float)1.38339, (float)1.33079, (
	    float)1.3598, (float)1.62823, (float)1.43873, (float)1.3389, (
	    float)1.32794, (float)1.37918, (float)1.62823, {0}, (float)1.9392,
	     (float)1.66764, (float)1.48511, (float)1.34514, (float)1.48927, (
	    float)1.9392, (float)1.42925, (float)1.34587, (float)1.35129, (
	    float)1.5427, (float)1.57895, (float)1.37317, (float)1.32921, (
	    float)1.35979, (float)1.62823, (float)1.43873, (float)1.3389, (
	    float)1.32794, (float)1.37918, (float)1.62823, {0}, (float)1.9392,
	     (float)1.66764, (float)1.48511, (float)1.34514, (float)1.48927, (
	    float)1.9392, (float)1.42444, (float)1.34469, (float)1.35128, (
	    float)1.5427, (float)1.56559, (float)1.37034, (float)1.32873, (
	    float)1.35979, (float)1.62823, (float)1.43873, (float)1.3389, (
	    float)1.32794, (float)1.37918, (float)1.62823 };


/* Table of constant values */

static integer c__1 = 1;
static real c_b2 = (float)0.;
static integer c__6 = 6;
static integer c__5 = 5;
static integer c__10 = 10;
static integer c__48 = 48;
static integer c__3 = 3;
static logical c_false = FALSE_;
static integer c__2 = 2;
static real c_b55 = (float).75;
static integer c__4 = 4;
static real c_b145 = (float).7;
static real c_b169 = (float).8;
static real c_b277 = (float).9;
static logical c_true = TRUE_;

/* Main program */ int MAIN__()
{
    /* Initialized data */

    static char abc[1*18+1] = "abcdefghijklmnopqr";
    static logical prnt[7] = { TRUE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,TRUE_ 
	    };
    static real accur = (float)0.;
    static char blanks[3+1] = "   ";
    static logical doprob[13] = { TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,
	    TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,TRUE_ };

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3, i__4;
    real r__1;

    /* Builtin functions */
    double asin(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(), 
	    i_indx(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double pow_ri(real *, integer *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_stop(char *
	    , ftnlen);

    /* Local variables */
    static integer j, k, lc;
    static real hl[49], pi;
    static integer iu, lu;
    static real uu[150]	/* was [10][5][3] */, u0u[50]	/* was [10][5] */;
    static integer iod;
    static real phi[3];
    static integer iss;
    static real umu[10], phi0, umu0;
    static integer icas;
    static real dfdt[5];
    static integer nphi;
    static real uavg[5], flup[5];
    static integer ntau;
    static real pmom[294]	/* was [49][6] */, utau[5];
    static integer nlyr, numu, nstr;
    static real fbeam;
    static integer ibcnd;
    static real dtauc[6], ssalb[6];
    static logical plank;
    static real rfldn[5], btemp, temis;
    static char title[100];
    static integer nprob;
    static real fisot, ttemp, cmpuu[150]	/* was [5][10][3] */, albmed[
	    10], albedo;
    static char header[127];
    static real cmpdfd[5];
    static logical lamber, deltam;
    static real cmpfdn[5], cmpfir[5], rfldir[5];
    static logical azmavg;
    static real trnmed[10], cmpfup[5], temper[7];
    extern /* Subroutine */ int getmom_(integer *, real *, integer *, real *);
    static logical usrang;
    static integer lentit;
    static logical onlyfl;
    extern /* Subroutine */ int disort_(integer *, real *, real *, real *, 
	    real *, real *, real *, logical *, integer *, real *, integer *, 
	    logical *, integer *, real *, integer *, real *, integer *, real *
	    , real *, real *, real *, logical *, real *, real *, real *, real 
	    *, real *, logical *, logical *, logical *, real *, logical *, 
	    char *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, ftnlen), errmsg_(char *, logical *, ftnlen);
    static real wvnmhi;
    extern /* Subroutine */ int prtfin_(real *, integer *, real *, integer *, 
	    real *, integer *, integer *, integer *, logical *, logical *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *, integer *, integer *);
    static real wvnmlo;
    static logical usrtau;

    /* Fortran I/O blocks */
    static icilist io___32 = { 0, header, 0, "(3A,F9.5,A,F5.2)", 127, 1 };
    static icilist io___53 = { 0, header, 0, "(3A,F5.2,A,F9.6,A,F4.2)", 127, 
	    1 };
    static icilist io___54 = { 0, header, 0, "(3A,F9.5,A,F5.2)", 127, 1 };
    static icilist io___56 = { 0, title, 0, "(3A)", 100, 1 };
    static icilist io___58 = { 0, title, 0, "(3A)", 100, 1 };
    static icilist io___59 = { 0, title, 0, "(3A)", 100, 1 };
    static icilist io___61 = { 0, title, 0, "(3A)", 100, 1 };


/*    Runs test problems for DISORT and checks answers. These */
/*    problems test almost all logical branches in DISORT. */
/*    As distributed, runs all problems.  To run just a subset, */
/*    change DATA statement for DOPROB;  e.g. to run just problem 4, */
/*    set DOPROB(4) TRUE and all others FALSE. */
/*    It is HIGHLY recommended that you use the code below as a template */
/*    for creating your own CALLs to DISORT, rather than starting from */
/*    scratch.  This will prevent a lot of mistakes and ensure that every */
/*    input argument gets a value.  Note in particular how GETMOM is */
/*    sometimes called to fill an array section of PMOM (for one layer); */
/*    several people have done this incorrectly in attempting to write it */
/*    ab initio (passing array sections for arrays that do not start at */
/*    element 1 is tricky). */
/*    Note that the ratio to the 'correct answer' may occasionally be */
/*    significantly different from unity -- even so different that */
/*    the ratio just prints as ****** rather than a number.  However, */
/*    this mostly occurs for values of flux or intensity that are very */
/*    small compared to the forcing functions (that is, small compared */
/*    to internal thermal emission and/or radiation incident at the */
/*    boundaries).  The printed number 'SERIOUSLY NON-UNIT RATIOS' */
/*    attempts to count just the cases where there is a real disagreement */
/*    and not those where quantitites are down at their noise level */
/*    (defined as 10^(-6) times their maximum value). */
/*    Further documentation can be found in the file DISOTEST.doc. */
/*  Routines called : */
/*    DISORT:   The discrete ordinates radiative transfer program */
/*    GETMOM:   Sets phase function Legendre coefficients */
/*    PRTFIN:   Prints fluxes and intensities and their ratios to */
/*              the correct values */
/*    CHEKDO:   Data block containing correct fluxes and intensities */
/*    RATIO :   Ratio of calculated to correct value with underflow */
/*              and overflow protection (kept in file DISORT.f) */
/* +---------------------------------------------------------------------+ */
/*     IMPLICIT  NONE */
/*                 ** DISORT I/O specifications ** */
/* +---------------------------------------------------------------------+ */
/*                     ** Correct answers ** */
/* +---------------------------------------------------------------------+ */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    pi = (float)2. * asin((float)1.);
    if (doprob[0]) {
/* ********************************************************************** */
/* ****  Test Problem 1:  Isotropic Scattering                       **** */
/* ****  (Compare to Ref. VH1, Table 12)                             **** */
/* ********************************************************************** */
	nstr = 16;
	nlyr = 1;
	getmom_(&c__1, &c_b2, &nstr, pmom);
	usrtau = TRUE_;
	ntau = 2;
	utau[0] = (float)0.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = (float)-1.;
	umu[1] = (float)-.5;
	umu[2] = (float)-.1;
	umu[3] = (float).1;
	umu[4] = (float).5;
	umu[5] = (float)1.;
	nphi = 1;
	phi[0] = (float)0.;
	ibcnd = 0;
	umu0 = (float).1;
	phi0 = (float)0.;
	lamber = TRUE_;
	albedo = (float)0.;
	deltam = FALSE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 6; ++icas) {
	    if (icas == 1) {
		utau[1] = (float).03125;
		ssalb[0] = (float).2;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 2) {
		utau[1] = (float).03125;
		ssalb[0] = (float)1.;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 3) {
		utau[1] = (float).03125;
		ssalb[0] = (float).99;
		fbeam = (float)0.;
		fisot = (float)1.;
	    } else if (icas == 4) {
		utau[1] = (float)32.;
		ssalb[0] = (float).2;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 5) {
		utau[1] = (float)32.;
		ssalb[0] = (float)1.;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 6) {
		utau[1] = (float)32.;
		ssalb[0] = (float).99;
		fbeam = (float)0.;
		fisot = (float)1.;
	    }
	    dtauc[0] = utau[1];
	    s_wsfi(&io___32);
	    do_fio(&c__1, "Test Case No. 1", (ftnlen)15);
	    do_fio(&c__1, abc + (icas - 1), (ftnlen)1);
	    do_fio(&c__1, ":  Isotropic Scattering, Ref. VH1, Table 12:  b =",
		     (ftnlen)49);
	    do_fio(&c__1, (char *)&utau[1], (ftnlen)sizeof(real));
	    do_fio(&c__1, ", a =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ssalb[0], (ftnlen)sizeof(real));
	    e_wsfi();
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    azmavg = TRUE_;
	    nprob = 1;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L10: */
	}
    }
    if (doprob[1]) {
/* ********************************************************************** */
/* ****  Test Problem 2:  Rayleigh Scattering, Beam Source           **** */
/* ****  (Compare To Ref. SW, Table 1)                               **** */
/* ********************************************************************** */
	nstr = 16;
	nlyr = 1;
	getmom_(&c__2, &c_b2, &nstr, pmom);
	usrtau = TRUE_;
	ntau = 2;
	utau[0] = (float)0.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = (float)-.981986;
	umu[1] = (float)-.538263;
	umu[2] = (float)-.018014;
	umu[3] = (float).018014;
	umu[4] = (float).538263;
	umu[5] = (float).981986;
	nphi = 1;
	phi[0] = (float)0.;
	ibcnd = 0;
	fbeam = pi;
	umu0 = (float).080442;
	phi0 = (float)0.;
	fisot = (float)0.;
	lamber = TRUE_;
	albedo = (float)0.;
	deltam = FALSE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	icas = 0;
	for (iod = 1; iod <= 2; ++iod) {
	    if (iod == 1) {
		utau[1] = (float).2;
	    }
	    if (iod == 2) {
		utau[1] = (float)5.;
	    }
	    dtauc[0] = utau[1];
	    for (iss = 1; iss <= 2; ++iss) {
		if (iss == 1) {
		    ssalb[0] = (float).5;
		}
		if (iss == 2) {
		    ssalb[0] = (float)1.;
		}
		++icas;
		s_wsfi(&io___53);
		do_fio(&c__1, "Test Case No. 2", (ftnlen)15);
		do_fio(&c__1, abc + (icas - 1), (ftnlen)1);
		do_fio(&c__1, ", Rayleigh Scattering, Ref. SW, Table 1:  tau\
 =", (ftnlen)47);
		do_fio(&c__1, (char *)&utau[1], (ftnlen)sizeof(real));
		do_fio(&c__1, ", mu0 =", (ftnlen)7);
		do_fio(&c__1, (char *)&umu0, (ftnlen)sizeof(real));
		do_fio(&c__1, ", ss-albedo =", (ftnlen)13);
		do_fio(&c__1, (char *)&ssalb[0], (ftnlen)sizeof(real));
		e_wsfi();
		disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
			usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &
			nphi, phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &
			lamber, &albedo, hl, &btemp, &ttemp, &temis, &deltam, 
			&plank, &onlyfl, &accur, prnt, header, &c__6, &c__5, &
			c__10, &c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg,
			 uu, u0u, albmed, trnmed, (ftnlen)127);
		azmavg = FALSE_;
		nprob = 2;
		if (nprob > 9 || icas > 8) {
		    errmsg_("Out of bounds in exact-answer arrays", &c_false, 
			    (ftnlen)36);
		}
		prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
			onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
			dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
			dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
			dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
			dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
			dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 
			+ 1) * 5 - 1405], &c__5, &c__10, &c__3);
/* L20: */
	    }
	}
    }
    if (doprob[2]) {
/* ********************************************************************** */
/* ****  Test Problem 3:  Henyey-Greenstein Scattering               **** */
/* ****  (Compare To Ref. VH2, Table 35)                             **** */
/* ********************************************************************** */
	nstr = 16;
	nlyr = 1;
	getmom_(&c__3, &c_b55, &nstr, pmom);
	usrtau = TRUE_;
	ntau = 2;
	utau[0] = (float)0.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = (float)-1.;
	umu[1] = (float)-.5;
	umu[2] = (float)-.1;
	umu[3] = (float).1;
	umu[4] = (float).5;
	umu[5] = (float)1.;
	nphi = 1;
	phi[0] = (float)0.;
	ibcnd = 0;
	umu0 = (float).1;
	phi0 = (float)0.;
	deltam = TRUE_;
	lamber = TRUE_;
	onlyfl = FALSE_;
	albedo = (float)0.;
	plank = FALSE_;
	for (icas = 1; icas <= 6; ++icas) {
	    if (icas == 1) {
		utau[1] = (float)1.;
		ssalb[0] = (float).2;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 2) {
		utau[1] = (float)1.;
		ssalb[0] = (float)1.;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 3) {
		utau[1] = (float)1.;
		ssalb[0] = (float).99;
		fbeam = (float)0.;
		fisot = (float)1.;
	    } else if (icas == 4) {
		utau[1] = (float)8.;
		ssalb[0] = (float).2;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 5) {
		utau[1] = (float)8.;
		ssalb[0] = (float)1.;
		fbeam = pi / umu0;
		fisot = (float)0.;
	    } else if (icas == 6) {
		utau[1] = (float)8.;
		ssalb[0] = (float).99;
		fbeam = (float)0.;
		fisot = (float)1.;
	    }
	    dtauc[0] = utau[1];
	    s_wsfi(&io___54);
	    do_fio(&c__1, "Test Case No. 3", (ftnlen)15);
	    do_fio(&c__1, abc + (icas - 1), (ftnlen)1);
	    do_fio(&c__1, ", Henyey-Greenstein Scattering, Ref. VH2, Table 3\
5, g = 0.75, b =", (ftnlen)65);
	    do_fio(&c__1, (char *)&utau[1], (ftnlen)sizeof(real));
	    do_fio(&c__1, ", a =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ssalb[0], (ftnlen)sizeof(real));
	    e_wsfi();
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    azmavg = TRUE_;
	    nprob = 3;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L30: */
	}
    }
    if (doprob[3]) {
/* ********************************************************************** */
/* ****  Test Problem 4:  Haze-L Scattering, Beam Source             **** */
/* ****  (Compare to Ref. GS, Tables 12-16)                          **** */
/* ********************************************************************** */
	nstr = 32;
	nlyr = 1;
	getmom_(&c__4, &c_b2, &nstr, pmom);
	dtauc[0] = (float)1.;
	usrtau = TRUE_;
	ntau = 3;
	utau[0] = (float)0.;
	utau[1] = (float).5;
	utau[2] = (float)1.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = (float)-1.;
	umu[1] = (float)-.5;
	umu[2] = (float)-.1;
	umu[3] = (float).1;
	umu[4] = (float).5;
	umu[5] = (float)1.;
	ibcnd = 0;
	fbeam = pi;
	phi0 = (float)0.;
	fisot = (float)0.;
	lamber = TRUE_;
	albedo = (float)0.;
	deltam = TRUE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    s_wsfi(&io___56);
	    do_fio(&c__1, "Test Case No. 4", (ftnlen)15);
	    do_fio(&c__1, abc + (icas - 1), (ftnlen)1);
	    do_fio(&c__1, ", Haze-L Scattering, Ref. GS, Table ", (ftnlen)36);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		ssalb[0] = (float)1.;
		nphi = 1;
		phi[0] = (float)0.;
		umu0 = (float)1.;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 3, a__1[1] = " 12";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 2) {
		ssalb[0] = (float).9;
		nphi = 1;
		phi[0] = (float)0.;
		umu0 = (float)1.;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 3, a__1[1] = " 13";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 3) {
		ssalb[0] = (float).9;
		nphi = 3;
		phi[0] = (float)0.;
		phi[1] = (float)90.;
		phi[2] = (float)180.;
		umu0 = (float).5;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 6, a__1[1] = " 14-16";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    azmavg = FALSE_;
	    nprob = 4;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L40: */
	}
    }
    if (doprob[4]) {
/* ********************************************************************** */
/* ****  Test Problem 5:  Cloud C.1 Scattering, Beam Source          **** */
/* ****  (Compare to Ref. GS, Tables 19-20)                          **** */
/* ********************************************************************** */
	nstr = 48;
	nlyr = 1;
	getmom_(&c__5, &c_b2, &nstr, pmom);
	dtauc[0] = (float)64.;
	usrtau = TRUE_;
	ntau = 3;
	usrang = TRUE_;
	numu = 6;
	umu[0] = (float)-1.;
	umu[1] = (float)-.5;
	umu[2] = (float)-.1;
	umu[3] = (float).1;
	umu[4] = (float).5;
	umu[5] = (float)1.;
	nphi = 1;
	phi[0] = (float)0.;
	ibcnd = 0;
	fbeam = pi;
	umu0 = (float)1.;
	phi0 = (float)0.;
	fisot = (float)0.;
	lamber = TRUE_;
	albedo = (float)0.;
	deltam = TRUE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    s_wsfi(&io___58);
	    do_fio(&c__1, "Test Case No. 5", (ftnlen)15);
	    do_fio(&c__1, abc + (icas - 1), (ftnlen)1);
	    do_fio(&c__1, ", Cloud C.1 Scattering, Ref. GS, Table ", (ftnlen)
		    39);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		utau[0] = (float)0.;
		utau[1] = (float)32.;
		utau[2] = (float)64.;
		ssalb[0] = (float)1.;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 3, a__1[1] = " 19";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    }
	    if (icas == 2) {
		utau[0] = (float)3.2;
		utau[1] = (float)12.8;
		utau[2] = (float)48.;
		ssalb[0] = (float).9;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 3, a__1[1] = " 20";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    azmavg = FALSE_;
	    nprob = 5;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L50: */
	}
    }
    if (doprob[5]) {
/* ********************************************************************** */
/* ****  Test Problem 6:  No Scattering, Increasingly Complex Sources**** */
/* ********************************************************************** */
	nstr = 16;
	nlyr = 1;
	ssalb[0] = (float)0.;
	wvnmlo = (float)0.;
	wvnmhi = (float)5e4;
	usrtau = TRUE_;
	usrang = TRUE_;
	numu = 4;
	umu[0] = (float)-1.;
	umu[1] = (float)-.1;
	umu[2] = (float).1;
	umu[3] = (float)1.;
	nphi = 1;
	phi[0] = (float)90.;
	ibcnd = 0;
	fbeam = (float)200.;
	umu0 = (float).5;
	phi0 = (float)0.;
	fisot = (float)0.;
	temis = (float)1.;
	onlyfl = FALSE_;
	deltam = FALSE_;
	for (icas = 1; icas <= 8; ++icas) {
	    s_wsfi(&io___59);
	    do_fio(&c__1, "Test Case No. 6", (ftnlen)15);
	    do_fio(&c__1, abc + (icas - 1), (ftnlen)1);
	    do_fio(&c__1, ": No Scattering; Source = Beam", (ftnlen)30);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		ntau = 2;
		utau[0] = (float)0.;
		utau[1] = (float)0.;
	    } else if (icas > 1) {
		ntau = 3;
		utau[0] = (float)0.;
		utau[1] = (float).5;
		utau[2] = (float)1.;
	    }
	    if (icas == 1) {
/*                                    ** Transparent medium, beam source */
		dtauc[0] = (float)0.;
		lamber = TRUE_;
		albedo = (float)0.;
		plank = FALSE_;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 19, a__1[1] = "; Bottom Albedo = 0";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 2) {
/*                                    ** Add some optical depth */
		dtauc[0] = (float)1.;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 19, a__1[1] = "; Bottom Albedo = 0";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 3) {
/*                                   ** Add some isotropic reflection */
		lamber = TRUE_;
		albedo = (float).5;
		plank = FALSE_;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 27, a__1[1] = "; Bottom Albedo=0.5 Lambert";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 4) {
/*                                   ** Use non-isotropic reflection */
		dtauc[0] = (float)1.;
		lamber = FALSE_;
		i__2 = nstr;
		for (k = 0; k <= i__2; ++k) {
		    hl[k] = pow_ri(&c_b145, &k);
/* L59: */
		}
		plank = FALSE_;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 29, a__1[1] = "; Bottom Albedo = Non-Lambert";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 5) {
/*                                   ** Add some bottom-boundary emission */
		dtauc[0] = (float)1.;
		temper[0] = (float)0.;
		temper[1] = (float)0.;
		lamber = FALSE_;
		btemp = (float)300.;
		ttemp = (float)0.;
		plank = TRUE_;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 41, a__1[1] = ", Bottom Emission; Bott Alb = Non-L\
ambert";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 6) {
/*                                   ** Add some top-boundary diffuse */
/*                                      incidence (prescribed + emitted) */
		dtauc[0] = (float)1.;
		temper[0] = (float)0.;
		temper[1] = (float)0.;
		fisot = (float)100. / pi;
		lamber = FALSE_;
		btemp = (float)300.;
		ttemp = (float)250.;
		plank = TRUE_;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 45, a__1[1] = ", Bottom+Top Emission; Bott Alb = N\
on-Lambert";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 7) {
/*                                   ** Add some internal emission */
		dtauc[0] = (float)1.;
		temper[0] = (float)250.;
		temper[1] = (float)300.;
		lamber = FALSE_;
		btemp = (float)300.;
		ttemp = (float)250.;
		plank = TRUE_;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 54, a__1[1] = ", Bottom+Top+Internal Emission; Bot\
t Alb = Non-Lambert";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 8) {
/*                                   ** Increase the optical depth */
		dtauc[0] = (float)10.;
		temper[0] = (float)250.;
		temper[1] = (float)300.;
		utau[0] = (float)0.;
		utau[1] = (float)1.;
		utau[2] = (float)10.;
		lamber = FALSE_;
		btemp = (float)300.;
		ttemp = (float)250.;
		plank = TRUE_;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 54, a__1[1] = ", Bottom+Top+Internal Emission; Bot\
t Alb = Non-Lambert";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    azmavg = FALSE_;
	    nprob = 6;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L60: */
	}
    }
    if (doprob[6]) {
/* ********************************************************************** */
/* ****  Test Problem 7:  Absorption + Scattering + All Possible     **** */
/* ****  Sources, Various Surface Reflectivities ( One Layer )       **** */
/* ********************************************************************** */
	nstr = 12;
	nlyr = 1;
	getmom_(&c__3, &c_b169, &nstr, pmom);
	dtauc[0] = (float)1.;
	ssalb[0] = (float).5;
	temper[0] = (float)300.;
	temper[1] = (float)200.;
	wvnmlo = (float)0.;
	wvnmhi = (float)5e4;
	usrtau = TRUE_;
	ntau = 3;
	utau[0] = (float)0.;
	utau[1] = (float).5;
	utau[2] = (float)1.;
	usrang = TRUE_;
	numu = 4;
	umu[0] = (float)-1.;
	umu[1] = (float)-.1;
	umu[2] = (float).1;
	umu[3] = (float)1.;
	nphi = 2;
	phi[0] = (float)0.;
	phi[1] = (float)90.;
	ibcnd = 0;
	fbeam = (float)200.;
	umu0 = (float).5;
	phi0 = (float)0.;
	fisot = (float)100.;
	btemp = (float)320.;
	ttemp = (float)100.;
	temis = (float)1.;
	deltam = TRUE_;
	plank = TRUE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    s_wsfi(&io___61);
	    do_fio(&c__1, "Test Case No. 7", (ftnlen)15);
	    do_fio(&c__1, abc + (icas - 1), (ftnlen)1);
	    do_fio(&c__1, ": Absorption + Henyey-Greenstein Scattering, All \
Sources", (ftnlen)56);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		lamber = TRUE_;
		albedo = (float)0.;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 19, a__1[1] = ", Bottom Albedo = 0";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 2) {
		lamber = TRUE_;
		albedo = (float)1.;
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 19, a__1[1] = ", Bottom Albedo = 1";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    } else if (icas == 3) {
		lamber = FALSE_;
		i__2 = nstr;
		for (k = 0; k <= i__2; ++k) {
		    hl[k] = pow_ri(&c_b145, &k);
/* L69: */
		}
/* Writing concatenation */
		i__1[0] = lentit, a__1[0] = title;
		i__1[1] = 30, a__1[1] = ", Bottom Albedo = BDR Function";
		s_cat(header, a__1, i__1, &c__2, (ftnlen)127);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    azmavg = FALSE_;
	    nprob = 7;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L70: */
	}
    }
    if (doprob[7]) {
/* ********************************************************************** */
/* ****  Test Problem 8:  Absorbing/Isotropic-Scattering Medium      **** */
/* ****  With Two Computational Layers                               **** */
/* **** (Compare Fluxes To Ref. OS, Table 1)                         **** */
/* ********************************************************************** */
	nstr = 8;
	nlyr = 2;
	getmom_(&c__1, &c_b2, &nstr, pmom);
	getmom_(&c__1, &c_b2, &nstr, &pmom[49]);
	usrtau = TRUE_;
	usrang = TRUE_;
	numu = 4;
	umu[0] = (float)-1.;
	umu[1] = (float)-.2;
	umu[2] = (float).2;
	umu[3] = (float)1.;
	nphi = 1;
	phi[0] = (float)60.;
	ibcnd = 0;
	fbeam = (float)0.;
	fisot = (float)1. / pi;
	lamber = TRUE_;
	albedo = (float)0.;
	plank = FALSE_;
	deltam = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    if (icas == 1) {
		dtauc[0] = (float).25;
		dtauc[1] = (float).25;
		ssalb[0] = (float).5;
		ssalb[1] = (float).3;
		ntau = 3;
		utau[0] = (float)0.;
		utau[1] = (float).25;
		utau[2] = (float).5;
		s_copy(header, "Test Case No. 8A:  Ref. OS, Table 1, Line 4 \
(Two Inhomogeneous Layers)", (ftnlen)127, (ftnlen)70);
	    } else if (icas == 2) {
		dtauc[0] = (float).25;
		dtauc[1] = (float).25;
		ssalb[0] = (float).8;
		ssalb[1] = (float).95;
		ntau = 3;
		utau[0] = (float)0.;
		utau[1] = (float).25;
		utau[2] = (float).5;
		s_copy(header, "Test Case No. 8b:  Ref. OS, Table 1, Line 1 \
(Two Inhomogeneous Layers)", (ftnlen)127, (ftnlen)70);
	    } else if (icas == 3) {
		dtauc[0] = (float)1.;
		dtauc[1] = (float)2.;
		ssalb[0] = (float).8;
		ssalb[1] = (float).95;
		ntau = 3;
		utau[0] = (float)0.;
		utau[1] = (float)1.;
		utau[2] = (float)3.;
		s_copy(header, "Test Case No. 8c:  Ref. OS, Table 1, Line 13\
 (Two Inhomogeneous Layers)", (ftnlen)127, (ftnlen)71);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    nprob = 8;
	    azmavg = FALSE_;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L80: */
	}
    }
    if (doprob[8]) {
/* ********************************************************************** */
/* ****  Test Problem 9:  General Emitting/Absorbing/Scattering      **** */
/* ****  Medium with Every Computational Layer Different.            **** */
/* **** (Compare 9a,b Fluxes to Ref. DGIS, Tables VI-VII, beta = 0)  **** */
/* ********************************************************************** */
	nstr = 8;
	nlyr = 6;
	i__2 = nlyr;
	for (lc = 1; lc <= i__2; ++lc) {
	    dtauc[lc - 1] = (real) lc;
	    ssalb[lc - 1] = lc * (float).05 + (float).6;
/* L86: */
	}
	usrtau = TRUE_;
	ntau = 5;
	utau[0] = (float)0.;
	utau[1] = (float)1.05;
	utau[2] = (float)2.1;
	utau[3] = (float)6.;
	utau[4] = (float)21.;
	usrang = TRUE_;
	numu = 4;
	umu[0] = (float)-1.;
	umu[1] = (float)-.2;
	umu[2] = (float).2;
	umu[3] = (float)1.;
	nphi = 1;
	phi[0] = (float)60.;
	ibcnd = 0;
	fbeam = (float)0.;
	fisot = (float)1. / pi;
	lamber = TRUE_;
	deltam = TRUE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    if (icas == 1) {
		i__2 = nlyr;
		for (lc = 1; lc <= i__2; ++lc) {
		    getmom_(&c__1, &c_b2, &nstr, &pmom[lc * 49 - 49]);
/* L87: */
		}
		albedo = (float)0.;
		plank = FALSE_;
		s_copy(header, "Test Case No. 9a:  Ref. DGIS, Tables VI-VII,\
 beta=l=0 (multiple inhomogeneous layers)", (ftnlen)127, (ftnlen)85);
	    } else if (icas == 2) {
		pmom[0] = (float)1.;
		pmom[1] = (float).66971999999999998;
		pmom[2] = (float).31267800000000001;
		pmom[3] = (float).096295714285714276;
		pmom[4] = (float).024683333333333331;
		pmom[5] = (float).0042954545454545459;
		pmom[6] = (float)5.1615384615384609e-4;
		pmom[7] = (float)4.5333333333333335e-5;
		pmom[8] = (float)2.9411764705882355e-6;
		i__2 = nlyr;
		for (lc = 2; lc <= i__2; ++lc) {
		    for (k = 0; k <= 8; ++k) {
			pmom[k + lc * 49 - 49] = pmom[k];
/* L88: */
		    }
		}
		s_copy(header, "Test Case No. 9b:  Ref. DGIS, Tables VI-VII,\
 beta=0,l=8 (multiple inhomogeneous layers)", (ftnlen)127, (ftnlen)87);
	    } else if (icas == 3) {
		temper[0] = (float)600.;
		i__2 = nlyr;
		for (lc = 1; lc <= i__2; ++lc) {
		    r__1 = (real) lc / (float)7.;
		    getmom_(&c__3, &r__1, &nstr, &pmom[lc * 49 - 49]);
		    temper[lc] = lc * (float)10. + (float)600.;
/* L89: */
		}
		nphi = 3;
		phi[0] = (float)60.;
		phi[1] = (float)120.;
		phi[2] = (float)180.;
		fbeam = pi;
		umu0 = (float).5;
		phi0 = (float)0.;
		fisot = (float)1.;
		albedo = (float).5;
		plank = TRUE_;
		wvnmlo = (float)999.;
		wvnmhi = (float)1e3;
		btemp = (float)700.;
		ttemp = (float)550.;
		temis = (float)1.;
		s_copy(header, "Test Case No. 9c:  Generalization of 9A to i\
nclude all possible complexity", (ftnlen)127, (ftnlen)74);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    nprob = 9;
	    azmavg = FALSE_;
	    if (nprob > 9 || icas > 8) {
		errmsg_("Out of bounds in exact-answer arrays", &c_false, (
			ftnlen)36);
	    }
	    prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
		    onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, &
		    dochek_1.tstfir[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfdn[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstfup[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstdfd[(icas + (nprob << 3)) * 5 - 45], &
		    dochek_1.tstuu[(((icas + (nprob << 3)) * 3 + 1) * 10 + 1) 
		    * 5 - 1405], &c__5, &c__10, &c__3);
/* L90: */
	}
    }
    if (doprob[9]) {
/* ********************************************************************** */
/* ****  Test Problem 10: Compare USRANG = True With USRANG = False  **** */
/* ****  take Problem 9c (our most general case) but only 4 Streams  **** */
/* ********************************************************************** */
	nstr = 4;
	nlyr = 6;
	temper[0] = (float)600.;
	i__2 = nlyr;
	for (lc = 1; lc <= i__2; ++lc) {
	    dtauc[lc - 1] = (real) lc;
	    ssalb[lc - 1] = lc * (float).05 + (float).6;
	    r__1 = (real) lc / (nlyr + 1);
	    getmom_(&c__3, &r__1, &nstr, &pmom[lc * 49 - 49]);
	    temper[lc] = lc * (float)10. + (float)600.;
/* L97: */
	}
	usrtau = TRUE_;
	ntau = 3;
	utau[0] = (float)0.;
	utau[1] = (float)2.1;
	utau[2] = (float)21.;
	nphi = 2;
	phi[0] = (float)60.;
	phi[1] = (float)120.;
	ibcnd = 0;
	fbeam = pi;
	umu0 = (float).5;
	phi0 = (float)0.;
	fisot = (float)1.;
	lamber = TRUE_;
	deltam = TRUE_;
	albedo = (float).5;
	plank = TRUE_;
	wvnmlo = (float)999.;
	wvnmhi = (float)1e3;
	btemp = (float)700.;
	ttemp = (float)550.;
	temis = (float)1.;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    if (icas == 1) {
		usrang = TRUE_;
		numu = 4;
		umu[0] = (float)-.788675129;
		umu[1] = (float)-.211324871;
		umu[2] = (float).211324871;
		umu[3] = (float).788675129;
		prnt[1] = TRUE_;
		prnt[4] = TRUE_;
		s_copy(header, "Test Case No. 10a:  like 9c, USRANG = True", (
			ftnlen)127, (ftnlen)42);
	    } else if (icas == 2) {
		usrang = FALSE_;
		numu = 0;
		prnt[1] = FALSE_;
		prnt[4] = FALSE_;
		s_copy(header, "Test Case No. 10b:  like 9C, USRANG = False", 
			(ftnlen)127, (ftnlen)43);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    if (icas == 1) {
/*                               ** Save results to compare to case 2 */
		i__2 = ntau;
		for (lu = 1; lu <= i__2; ++lu) {
		    cmpfir[lu - 1] = rfldir[lu - 1];
		    cmpfdn[lu - 1] = rfldn[lu - 1];
		    cmpfup[lu - 1] = flup[lu - 1];
		    cmpdfd[lu - 1] = dfdt[lu - 1];
		    i__3 = numu;
		    for (iu = 1; iu <= i__3; ++iu) {
			i__4 = nphi;
			for (j = 1; j <= i__4; ++j) {
			    cmpuu[lu + (iu + j * 10) * 5 - 56] = uu[iu + (lu 
				    + j * 5) * 10 - 61];
/* L98: */
			}
		    }
		}
	    } else if (icas == 2) {
		azmavg = FALSE_;
		prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
			onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, 
			cmpfir, cmpfdn, cmpfup, cmpdfd, cmpuu, &c__5, &c__10, 
			&c__3);
	    }
/* L100: */
	}
    }
    if (doprob[10]) {
/* ********************************************************************** */
/* ****  Test Problem 11: Single-Layer vs. Multiple Layers           **** */
/* ****  11a: Results at user levels for one computational layer     **** */
/* ****  11b: Single layer of 11a subdivided into multiple           **** */
/* ****       computational layers at the 11a user levels            **** */
/* ********************************************************************** */
	nstr = 16;
	usrang = TRUE_;
	numu = 4;
	umu[0] = (float)-1.;
	umu[1] = (float)-.1;
	umu[2] = (float).1;
	umu[3] = (float)1.;
	nphi = 2;
	phi[0] = (float)0.;
	phi[1] = (float)90.;
	ibcnd = 0;
	fbeam = (float)1.;
	umu0 = (float).5;
	phi0 = (float)0.;
	fisot = (float).5 / pi;
	lamber = TRUE_;
	albedo = (float).5;
	deltam = FALSE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    if (icas == 1) {
		nlyr = 1;
		dtauc[0] = (float)1.;
		ssalb[0] = (float).9;
		getmom_(&c__1, &c_b2, &nstr, pmom);
		usrtau = TRUE_;
		ntau = 4;
		utau[0] = (float)0.;
		utau[1] = (float).05;
		utau[2] = (float).5;
		utau[3] = (float)1.;
		prnt[1] = TRUE_;
		prnt[4] = TRUE_;
		s_copy(header, "Test Case No. 11a: One Isotropic-Scattering \
Layer", (ftnlen)127, (ftnlen)49);
	    } else if (icas == 2) {
		nlyr = ntau - 1;
		i__4 = nlyr;
		for (lc = 1; lc <= i__4; ++lc) {
		    dtauc[lc - 1] = utau[lc] - utau[lc - 1];
		    ssalb[lc - 1] = (float).9;
		    getmom_(&c__1, &c_b2, &nstr, &pmom[lc * 49 - 49]);
/* L107: */
		}
		usrtau = FALSE_;
		prnt[1] = FALSE_;
		prnt[4] = FALSE_;
		s_copy(header, "Test Case No. 11b: Same as 11a but treated a\
s multiple layers", (ftnlen)127, (ftnlen)61);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    if (icas == 1) {
/*                               ** Save results to compare to case 2 */
		i__4 = ntau;
		for (lu = 1; lu <= i__4; ++lu) {
		    cmpfir[lu - 1] = rfldir[lu - 1];
		    cmpfdn[lu - 1] = rfldn[lu - 1];
		    cmpfup[lu - 1] = flup[lu - 1];
		    cmpdfd[lu - 1] = dfdt[lu - 1];
		    i__3 = numu;
		    for (iu = 1; iu <= i__3; ++iu) {
			i__2 = nphi;
			for (j = 1; j <= i__2; ++j) {
			    cmpuu[lu + (iu + j * 10) * 5 - 56] = uu[iu + (lu 
				    + j * 5) * 10 - 61];
/* L108: */
			}
		    }
		}
	    } else if (icas == 2) {
		azmavg = FALSE_;
		prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
			onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, 
			cmpfir, cmpfdn, cmpfup, cmpdfd, cmpuu, &c__5, &c__10, 
			&c__3);
	    }
/* L110: */
	}
    }
    if (doprob[11]) {
/* ********************************************************************** */
/* ****  Test Problem 12: Test Absorption-Optical-Depth Shortcut     **** */
/* ****  compares cases where the DISORT shortcut for absorption     **** */
/* ****  optical depth .GT. 10 is not used (12a), then is used (12b) **** */
/* ****  (this shortcut is only employed when  PLANK = False.)       **** */
/* ********************************************************************** */
	nstr = 20;
	usrang = TRUE_;
	numu = 4;
	umu[0] = (float)-1.;
	umu[1] = (float)-.1;
	umu[2] = (float).1;
	umu[3] = (float)1.;
	nphi = 1;
	phi[0] = (float)0.;
	ibcnd = 0;
	fbeam = (float)1.;
	umu0 = (float)1.;
	phi0 = (float)0.;
	fisot = (float)0.;
	lamber = TRUE_;
	albedo = (float)1.;
	deltam = TRUE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    if (icas == 1) {
		nlyr = 1;
		dtauc[0] = (float)20.1;
		ssalb[0] = (float).5;
		getmom_(&c__3, &c_b277, &nstr, pmom);
		usrtau = TRUE_;
		ntau = 4;
		utau[0] = (float)0.;
		utau[1] = (float)10.;
		utau[2] = (float)19.9;
		utau[3] = (float)20.1;
		prnt[1] = TRUE_;
		prnt[4] = TRUE_;
		s_copy(header, "Test Case No. 12a:  Overhead Beam Striking A\
bsorbing/Scattering Medium", (ftnlen)127, (ftnlen)70);
	    } else if (icas == 2) {
		nlyr = ntau - 1;
		i__2 = nlyr;
		for (lc = 1; lc <= i__2; ++lc) {
		    dtauc[lc - 1] = utau[lc] - utau[lc - 1];
		    ssalb[lc - 1] = (float).5;
		    getmom_(&c__3, &c_b277, &nstr, &pmom[lc * 49 - 49]);
/* L117: */
		}
		usrtau = FALSE_;
		prnt[1] = FALSE_;
		prnt[4] = FALSE_;
		s_copy(header, "Test Case No. 12b: Same as 12a but uses shor\
tcut for absorption optical depth .GT. 10", (ftnlen)127, (ftnlen)85);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
	    if (icas == 1) {
/*                               ** Save results to compare to case 2 */
		i__2 = ntau;
		for (lu = 1; lu <= i__2; ++lu) {
		    cmpfir[lu - 1] = rfldir[lu - 1];
		    cmpfdn[lu - 1] = rfldn[lu - 1];
		    cmpfup[lu - 1] = flup[lu - 1];
		    cmpdfd[lu - 1] = dfdt[lu - 1];
		    i__3 = numu;
		    for (iu = 1; iu <= i__3; ++iu) {
			i__4 = nphi;
			for (j = 1; j <= i__4; ++j) {
			    cmpuu[lu + (iu + j * 10) * 5 - 56] = uu[iu + (lu 
				    + j * 5) * 10 - 61];
/* L118: */
			}
		    }
		}
	    } else if (icas == 2) {
		azmavg = FALSE_;
		prtfin_(utau, &ntau, umu, &numu, phi, &nphi, &c__5, &c__10, &
			onlyfl, &azmavg, rfldir, rfldn, flup, dfdt, u0u, uu, 
			cmpfir, cmpfdn, cmpfup, cmpdfd, cmpuu, &c__5, &c__10, 
			&c__3);
	    }
/* L120: */
	}
    }
    if (doprob[12]) {
/* ********************************************************************** */
/* ****  Test Problem 13: Test shortcut for flux albedo, transmission *** */
/* **** ( shortcut gives flux albedo, transmission of entire medium  **** */
/* ****   as a function of sun angle )                               **** */
/* ****  13a,c = Shortcut;  13b,d = Brute Force Method               **** */
/* ********************************************************************** */
	nstr = 16;
	nphi = 0;
	phi0 = (float)0.;
	albedo = (float).5;
	deltam = TRUE_;
	for (icas = 1; icas <= 4; ++icas) {
	    if (icas == 1) {
		ibcnd = 1;
		nlyr = 1;
		dtauc[0] = (float)1.;
		ssalb[0] = (float).99;
		getmom_(&c__3, &c_b169, &nstr, pmom);
		prnt[5] = TRUE_;
		prnt[1] = FALSE_;
		usrang = TRUE_;
		numu = 1;
		umu[0] = (float).5;
		s_copy(header, "Test Case No. 13a:  Albedo and Transmissivit\
y from Shortcut, Single Layer", (ftnlen)127, (ftnlen)73);
	    } else if (icas == 2) {
		ibcnd = 0;
		usrtau = TRUE_;
		ntau = 2;
		utau[0] = (float)0.;
		utau[1] = (float)1.;
		umu0 = (float).5;
		fbeam = (float)1. / umu0;
		fisot = (float)0.;
		lamber = TRUE_;
		plank = FALSE_;
		onlyfl = TRUE_;
		prnt[5] = FALSE_;
		prnt[1] = TRUE_;
		s_copy(header, "Test Case No. 13b:  Albedo and Transmissivit\
y by Regular Method, Single Layer", (ftnlen)127, (ftnlen)77);
	    } else if (icas == 3) {
		ibcnd = 1;
		prnt[5] = TRUE_;
		prnt[1] = FALSE_;
		nlyr = 2;
		i__4 = nlyr;
		for (lc = 1; lc <= i__4; ++lc) {
		    dtauc[lc - 1] = (float)1. / nlyr;
		    getmom_(&c__3, &c_b169, &nstr, &pmom[lc * 49 - 49]);
/* L125: */
		}
		ssalb[0] = (float).99;
		ssalb[1] = (float).5;
		s_copy(header, "Test Case No. 13c:  Albedo and Transmissivit\
y from Shortcut, Multiple Layer", (ftnlen)127, (ftnlen)75);
	    } else if (icas == 4) {
		ibcnd = 0;
		prnt[5] = FALSE_;
		prnt[1] = TRUE_;
		s_copy(header, "Test Case No. 13d:  Albedo and Transmissivit\
y by Regular Method, Multiple Layer", (ftnlen)127, (ftnlen)79);
	    }
	    disort_(&nlyr, dtauc, ssalb, pmom, temper, &wvnmlo, &wvnmhi, &
		    usrtau, &ntau, utau, &nstr, &usrang, &numu, umu, &nphi, 
		    phi, &ibcnd, &fbeam, &umu0, &phi0, &fisot, &lamber, &
		    albedo, hl, &btemp, &ttemp, &temis, &deltam, &plank, &
		    onlyfl, &accur, prnt, header, &c__6, &c__5, &c__10, &
		    c__48, &c__3, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, 
		    albmed, trnmed, (ftnlen)127);
/* L130: */
	}
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Subroutine */ int getmom_(integer *iphas, real *gg, integer *nmom, real *
	pmom)
{
    /* Initialized data */

    static real hazelm[82] = { (float)2.4126,(float)3.23047,(float)3.37296,(
	    float)3.2315,(float)2.8935,(float)2.49594,(float)2.11361,(float)
	    1.74812,(float)1.44692,(float)1.17714,(float).96643,(float).78237,
	    (float).64114,(float).51966,(float).42563,(float).34688,(float)
	    .28351,(float).23317,(float).18963,(float).15788,(float).12739,(
	    float).10762,(float).08597,(float).07381,(float).05828,(float)
	    .05089,(float).03971,(float).03524,(float).0272,(float).02451,(
	    float).01874,(float).01711,(float).01298,(float).01198,(float)
	    .00904,(float).00841,(float).00634,(float).00592,(float).00446,(
	    float).00418,(float).00316,(float).00296,(float).00225,(float)
	    .0021,(float).0016,(float).0015,(float).00115,(float).00107,(
	    float)8.2e-4,(float)7.7e-4,(float)5.9e-4,(float)5.5e-4,(float)
	    4.3e-4,(float)4e-4,(float)3.1e-4,(float)2.9e-4,(float)2.3e-4,(
	    float)2.1e-4,(float)1.7e-4,(float)1.5e-4,(float)1.2e-4,(float)
	    1.1e-4,(float)9e-5,(float)8e-5,(float)6e-5,(float)6e-5,(float)
	    5e-5,(float)4e-5,(float)4e-5,(float)3e-5,(float)3e-5,(float)2e-5,(
	    float)2e-5,(float)2e-5,(float)1e-5,(float)1e-5,(float)1e-5,(float)
	    1e-5,(float)1e-5,(float)1e-5,(float)1e-5,(float)1e-5 };
    static real cldmom[299] = { (float)2.544,(float)3.883,(float)4.568,(float)
	    5.235,(float)5.887,(float)6.457,(float)7.177,(float)7.859,(float)
	    8.494,(float)9.286,(float)9.856,(float)10.615,(float)11.229,(
	    float)11.851,(float)12.503,(float)13.058,(float)13.626,(float)
	    14.209,(float)14.66,(float)15.231,(float)15.641,(float)16.126,(
	    float)16.539,(float)16.934,(float)17.325,(float)17.673,(float)
	    17.999,(float)18.329,(float)18.588,(float)18.885,(float)19.103,(
	    float)19.345,(float)19.537,(float)19.721,(float)19.884,(float)
	    20.024,(float)20.145,(float)20.251,(float)20.33,(float)20.401,(
	    float)20.444,(float)20.477,(float)20.489,(float)20.483,(float)
	    20.467,(float)20.427,(float)20.382,(float)20.31,(float)20.236,(
	    float)20.136,(float)20.036,(float)19.909,(float)19.785,(float)
	    19.632,(float)19.486,(float)19.311,(float)19.145,(float)18.949,(
	    float)18.764,(float)18.551,(float)18.348,(float)18.119,(float)
	    17.901,(float)17.659,(float)17.428,(float)17.174,(float)16.931,(
	    float)16.668,(float)16.415,(float)16.144,(float)15.883,(float)
	    15.606,(float)15.338,(float)15.058,(float)14.784,(float)14.501,(
	    float)14.225,(float)13.941,(float)13.662,(float)13.378,(float)
	    13.098,(float)12.816,(float)12.536,(float)12.257,(float)11.978,(
	    float)11.703,(float)11.427,(float)11.156,(float)10.884,(float)
	    10.618,(float)10.35,(float)10.09,(float)9.827,(float)9.574,(float)
	    9.318,(float)9.072,(float)8.822,(float)8.584,(float)8.34,(float)
	    8.11,(float)7.874,(float)7.652,(float)7.424,(float)7.211,(float)
	    6.99,(float)6.785,(float)6.573,(float)6.377,(float)6.173,(float)
	    5.986,(float)5.79,(float)5.612,(float)5.424,(float)5.255,(float)
	    5.075,(float)4.915,(float)4.744,(float)4.592,(float)4.429,(float)
	    4.285,(float)4.13,(float)3.994,(float)3.847,(float)3.719,(float)
	    3.58,(float)3.459,(float)3.327,(float)3.214,(float)3.09,(float)
	    2.983,(float)2.866,(float)2.766,(float)2.656,(float)2.562,(float)
	    2.459,(float)2.372,(float)2.274,(float)2.193,(float)2.102,(float)
	    2.025,(float)1.94,(float)1.869,(float)1.79,(float)1.723,(float)
	    1.649,(float)1.588,(float)1.518,(float)1.461,(float)1.397,(float)
	    1.344,(float)1.284,(float)1.235,(float)1.179,(float)1.134,(float)
	    1.082,(float)1.04,(float).992,(float).954,(float).909,(float).873,
	    (float).832,(float).799,(float).762,(float).731,(float).696,(
	    float).668,(float).636,(float).61,(float).581,(float).557,(float)
	    .53,(float).508,(float).483,(float).463,(float).44,(float).422,(
	    float).401,(float).384,(float).364,(float).349,(float).331,(float)
	    .317,(float).301,(float).288,(float).273,(float).262,(float).248,(
	    float).238,(float).225,(float).215,(float).204,(float).195,(float)
	    .185,(float).177,(float).167,(float).16,(float).151,(float).145,(
	    float).137,(float).131,(float).124,(float).118,(float).112,(float)
	    .107,(float).101,(float).097,(float).091,(float).087,(float).082,(
	    float).079,(float).074,(float).071,(float).067,(float).064,(float)
	    .06,(float).057,(float).054,(float).052,(float).049,(float).047,(
	    float).044,(float).042,(float).039,(float).038,(float).035,(float)
	    .034,(float).032,(float).03,(float).029,(float).027,(float).026,(
	    float).024,(float).023,(float).022,(float).021,(float).02,(float)
	    .018,(float).018,(float).017,(float).016,(float).015,(float).014,(
	    float).013,(float).013,(float).012,(float).011,(float).011,(float)
	    .01,(float).009,(float).009,(float).008,(float).008,(float).008,(
	    float).007,(float).007,(float).006,(float).006,(float).006,(float)
	    .005,(float).005,(float).005,(float).005,(float).004,(float).004,(
	    float).004,(float).004,(float).003,(float).003,(float).003,(float)
	    .003,(float).003,(float).003,(float).002,(float).002,(float).002,(
	    float).002,(float).002,(float).002,(float).002,(float).002,(float)
	    .002,(float).001,(float).001,(float).001,(float).001,(float).001,(
	    float).001,(float).001,(float).001,(float).001,(float).001,(float)
	    .001,(float).001,(float).001,(float).001,(float).001,(float).001,(
	    float).001,(float).001 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double pow_ri(real *, integer *);

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

/*        Calculate phase function Legendre expansion coefficients */
/*        in various special cases */
/*       INPUT: IPHAS   Phase function options */
/*                      1 : Isotropic */
/*                      2 : Rayleigh */
/*                      3 : Henyey-Greenstein with asymmetry factor GG */
/*                      4 : Haze L as specified by Garcia/Siewert */
/*                      5 : Cloud C.1 as specified by Garcia/Siewert */
/*              GG      Asymmetry factor for Henyey-Greenstein case */
/*              NMOM    Index of highest Legendre coefficient needed */
/*                        ( = number of streams 'NSTR'  chosen */
/*                         for the discrete ordinate method) */
/*      OUTPUT: PMOM(K)  Legendre expansion coefficients (K=0 to NMOM) */
/*                         (be sure to dimension '0:maxval' in calling */
/*                          program) */
/*      Reference:  Garcia, R. and C. Siewert, 1985: Benchmark Results */
/*                     in Radiative Transfer, Transp. Theory and Stat. */
/*                     Physics 14, 437-484, Tables 10 And 17 */
/* ------------------------------------------------------------------ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    if (*iphas < 1 || *iphas > 5) {
	errmsg_("GETMOM--bad input variable IPHAS", &c_true, (ftnlen)32);
    }
    if (*iphas == 3 && (*gg <= (float)-1. || *gg >= (float)1.)) {
	errmsg_("GETMOM--bad input variable GG", &c_true, (ftnlen)29);
    }
    if (*nmom < 2) {
	errmsg_("GETMOM--bad input variable NMOM", &c_true, (ftnlen)31);
    }
    pmom[0] = (float)1.;
    i__1 = *nmom;
    for (k = 1; k <= i__1; ++k) {
	pmom[k] = (float)0.;
/* L10: */
    }
    if (*iphas == 2) {
/*                                       ** Rayleigh phase function */
	pmom[2] = (float).1;
    } else if (*iphas == 3) {
/*                                       ** Henyey-Greenstein phase fcn */
	i__1 = *nmom;
	for (k = 1; k <= i__1; ++k) {
	    pmom[k] = pow_ri(gg, &k);
/* L20: */
	}
    } else if (*iphas == 4) {
/*                                        ** Haze-L phase function */
	i__1 = min(82,*nmom);
	for (k = 1; k <= i__1; ++k) {
	    pmom[k] = hazelm[k - 1] / ((k << 1) + 1);
/* L30: */
	}
    } else if (*iphas == 5) {
/*                                        ** Cloud C.1 phase function */
	i__1 = min(298,*nmom);
	for (k = 1; k <= i__1; ++k) {
	    pmom[k] = cldmom[k - 1] / ((k << 1) + 1);
/* L40: */
	}
    }
    return 0;
} /* getmom_ */

/* Subroutine */ int prtfin_(real *utau, integer *ntau, real *umu, integer *
	numu, real *phi, integer *nphi, integer *maxulv, integer *maxumu, 
	logical *onlyfl, logical *azmavg, real *rfldir, real *rfldn, real *
	flup, real *dfdt, real *u0u, real *uu, real *tstfir, real *tstfdn, 
	real *tstfup, real *tstdfd, real *tstuu, integer *maxtau, integer *
	maxmu, integer *maxaz)
{
    /* Format strings */
    static char fmt_300[] = "(//,1x,45(\002=\002),/,a,i4,a,/,1x,45(\002=\002\
))";

    /* System generated locals */
    integer tstuu_dim1, tstuu_dim2, tstuu_offset, u0u_dim1, u0u_offset, 
	    uu_dim1, uu_dim2, uu_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer j, iu, lu;
    static real rat1, rat2, rat3, rat4, umax, ratv[100];
    extern doublereal ratio_(real *, real *);
    static integer numbad;
    static real fnoise, flxmax;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static real unoise;

    /* Fortran I/O blocks */
    static cilist io___78 = { 0, 6, 0, "(//,A,/,A,/,A)", 0 };
    static cilist io___79 = { 0, 6, 0, "(0P,F11.4,1P,4E15.4)", 0 };
    static cilist io___84 = { 0, 6, 0, "(11X,4( '    (',F9.4,')'))", 0 };
    static cilist io___88 = { 0, 6, 0, "(//,A,//,A,/,A,8F14.5)", 0 };
    static cilist io___89 = { 0, 6, 0, "(0P,F10.3,1P,8E14.4)", 0 };
    static cilist io___91 = { 0, 6, 0, "(10X, 8(:,'   (',F9.4,')'))", 0 };
    static cilist io___93 = { 0, 6, 0, "(//,A,//,A,/,A,/,A,8(F10.1,4X))", 0 };
    static cilist io___94 = { 0, 6, 0, "(/,0P,F10.3,F8.3,1P,8E14.4)", 0 };
    static cilist io___95 = { 0, 6, 0, "(10X,0P,F8.3, 1P,8E14.4)", 0 };
    static cilist io___96 = { 0, 6, 0, "(18X, 8(:,'   (',F9.4,')'))", 0 };
    static cilist io___97 = { 0, 6, 0, fmt_300, 0 };


/*        Print DISORT results and, directly beneath them, their */
/*        ratios to the correct answers;  print number of non-unit */
/*        ratios that occur but try to count just the cases where */
/*        there is a real disagreement and not those where flux or */
/*        intensity are down at their noise level (defined as 10^(-6) */
/*        times their maximum value).  d(flux)/d(tau) is treated the */
/*        same as fluxes in this noise estimation even though it */
/*        is a different type of quantity (although with flux units). */
/*     INPUT :   TSTFIR  correct direct flux */
/*               TSTFDN  correct diffuse down flux */
/*               TSTFUP  correct diffuse up flux */
/*               TSTDFD  correct d(flux)/d(optical depth) */
/*               TSTUU   correct intensity or azim-avg intensity */
/*               AZMAVG  TRUE, check azimuthally-averaged intensity only */
/*                       FALSE, check full intensity */
/*               (remaining input = DISORT I/O variables) */
/* -------------------------------------------------------------------- */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
    /* Parameter adjustments */
    --utau;
    --umu;
    --phi;
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2);
    uu -= uu_offset;
    u0u_dim1 = *maxumu;
    u0u_offset = 1 + u0u_dim1;
    u0u -= u0u_offset;
    --rfldir;
    --rfldn;
    --flup;
    --dfdt;
    --tstfir;
    --tstfdn;
    --tstfup;
    --tstdfd;
    tstuu_dim1 = *maxtau;
    tstuu_dim2 = *maxmu;
    tstuu_offset = 1 + tstuu_dim1 * (1 + tstuu_dim2);
    tstuu -= tstuu_offset;

    /* Function Body */
    if (*ntau > *maxtau || *numu > *maxmu || *nphi > *maxaz) {
	errmsg_("PRTFIN--out of bounds in comparator arrays", &c_true, (
		ftnlen)42);
    }
    flxmax = (float)0.;
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
/* Computing MAX */
	r__1 = flxmax, r__2 = tstfir[lu], r__1 = max(r__1,r__2), r__2 = 
		tstfdn[lu], r__1 = max(r__1,r__2), r__2 = tstfup[lu];
	flxmax = dmax(r__1,r__2);
/* L5: */
    }
    fnoise = flxmax * (float)1e-6;
    if (flxmax <= (float)0.) {
	errmsg_("PRTFIN--all fluxes zero or negative", &c_false, (ftnlen)35);
    }
    if (fnoise <= (float)0.) {
	errmsg_("PRTFIN--all fluxes near underflowing", &c_false, (ftnlen)36);
    }
    numbad = 0;
    s_wsfe(&io___78);
    do_fio(&c__1, "                  <-------------- FLUXES -------------->", 
	    (ftnlen)56);
    do_fio(&c__1, "    Optical       Downward       Downward         Upward \
   d(Net Flux)", (ftnlen)71);
    do_fio(&c__1, "      Depth         Direct        Diffuse        Diffuse \
   / d(Op Dep)", (ftnlen)71);
    e_wsfe();
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	s_wsfe(&io___79);
	do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rfldir[lu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rfldn[lu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&flup[lu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dfdt[lu], (ftnlen)sizeof(real));
	e_wsfe();
	rat1 = ratio_(&rfldir[lu], &tstfir[lu]);
	rat2 = ratio_(&rfldn[lu], &tstfdn[lu]);
	rat3 = ratio_(&flup[lu], &tstfup[lu]);
	rat4 = ratio_(&dfdt[lu], &tstdfd[lu]);
	s_wsfe(&io___84);
	do_fio(&c__1, (char *)&rat1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rat2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rat3, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rat4, (ftnlen)sizeof(real));
	e_wsfe();
	if ((rat1 < (float).99 || rat1 > (float)1.01) && (r__1 = rfldir[lu], 
		dabs(r__1)) > fnoise) {
	    ++numbad;
	}
	if ((rat2 < (float).99 || rat2 > (float)1.01) && (r__1 = rfldn[lu], 
		dabs(r__1)) > fnoise) {
	    ++numbad;
	}
	if ((rat3 < (float).99 || rat3 > (float)1.01) && (r__1 = flup[lu], 
		dabs(r__1)) > fnoise) {
	    ++numbad;
	}
	if ((rat4 < (float).99 || rat4 > (float)1.01) && (r__1 = dfdt[lu], 
		dabs(r__1)) > fnoise) {
	    ++numbad;
	}
/* L10: */
    }
    if (*onlyfl) {
	goto L100;
    }
    if (*numu > 100 || *nphi > 100) {
	errmsg_("PRTFIN--increase parameter MAXRAT", &c_true, (ftnlen)33);
    }
    if (*azmavg) {
/*                                       ** Print az-avg intensities */
	if (*numu > 8) {
	    errmsg_("PRTFIN--az-avg-intensity FORMATs inadequate", &c_false, (
		    ftnlen)43);
	}
	umax = (float)0.;
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
/* Computing MAX */
		r__1 = umax, r__2 = tstuu[lu + (iu + tstuu_dim2) * tstuu_dim1]
			;
		umax = dmax(r__1,r__2);
/* L15: */
	    }
/* L16: */
	}
	unoise = umax * (float)1e-6;
	if (umax <= (float)0.) {
	    errmsg_("PRTFIN--all az-avg intensities zero or negative", &
		    c_false, (ftnlen)47);
	}
	if (unoise <= (float)0.) {
	    errmsg_("PRTFIN--all az-avg intensities near underflowing", &
		    c_false, (ftnlen)48);
	}
	s_wsfe(&io___88);
	do_fio(&c__1, " ********* AZIMUTHALLY AVERAGED INTENSITIES  *********"
		, (ftnlen)54);
	do_fio(&c__1, "   Optical   Polar Angle Cosines", (ftnlen)32);
	do_fio(&c__1, "     Depth", (ftnlen)10);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
	}
	e_wsfe();
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    s_wsfe(&io___89);
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		do_fio(&c__1, (char *)&u0u[iu + lu * u0u_dim1], (ftnlen)
			sizeof(real));
	    }
	    e_wsfe();
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		ratv[iu - 1] = ratio_(&u0u[iu + lu * u0u_dim1], &tstuu[lu + (
			iu + tstuu_dim2) * tstuu_dim1]);
		if ((ratv[iu - 1] < (float).99 || ratv[iu - 1] > (float)1.01) 
			&& (r__1 = u0u[iu + lu * u0u_dim1], dabs(r__1)) > 
			unoise) {
		    ++numbad;
		}
/* L20: */
	    }
	    s_wsfe(&io___91);
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		do_fio(&c__1, (char *)&ratv[iu - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/* L30: */
	}
    } else {
/*                                       ** Print intensities */
	if (*nphi > 8) {
	    errmsg_("PRTFIN--intensity FORMATs inadequate", &c_false, (ftnlen)
		    36);
	}
	umax = (float)0.;
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		i__3 = *nphi;
		for (j = 1; j <= i__3; ++j) {
/* Computing MAX */
		    r__1 = umax, r__2 = tstuu[lu + (iu + j * tstuu_dim2) * 
			    tstuu_dim1];
		    umax = dmax(r__1,r__2);
/* L34: */
		}
/* L35: */
	    }
/* L36: */
	}
	unoise = umax * (float)1e-6;
	if (umax <= (float)0.) {
	    errmsg_("PRTFIN--all intensities zero or negative", &c_false, (
		    ftnlen)40);
	}
	if (unoise <= (float)0.) {
	    errmsg_("PRTFIN--all intensities near underflowing", &c_false, (
		    ftnlen)41);
	}
	s_wsfe(&io___93);
	do_fio(&c__1, " ********  I N T E N S I T I E S  *********", (ftnlen)
		43);
	do_fio(&c__1, "             Polar   Azimuthal Angles (Degrees)", (
		ftnlen)47);
	do_fio(&c__1, "   Optical   Angle", (ftnlen)18);
	do_fio(&c__1, "     Depth  Cosine", (ftnlen)18);
	i__1 = *nphi;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&phi[j], (ftnlen)sizeof(real));
	}
	e_wsfe();
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		if (iu == 1) {
		    s_wsfe(&io___94);
		    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
		    i__3 = *nphi;
		    for (j = 1; j <= i__3; ++j) {
			do_fio(&c__1, (char *)&uu[iu + (lu + j * uu_dim2) * 
				uu_dim1], (ftnlen)sizeof(real));
		    }
		    e_wsfe();
		}
		if (iu > 1) {
		    s_wsfe(&io___95);
		    do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
		    i__3 = *nphi;
		    for (j = 1; j <= i__3; ++j) {
			do_fio(&c__1, (char *)&uu[iu + (lu + j * uu_dim2) * 
				uu_dim1], (ftnlen)sizeof(real));
		    }
		    e_wsfe();
		}
		i__3 = *nphi;
		for (j = 1; j <= i__3; ++j) {
		    ratv[j - 1] = ratio_(&uu[iu + (lu + j * uu_dim2) * 
			    uu_dim1], &tstuu[lu + (iu + j * tstuu_dim2) * 
			    tstuu_dim1]);
		    if ((ratv[j - 1] < (float).99 || ratv[j - 1] > (float)
			    1.01) && (r__1 = uu[iu + (lu + j * uu_dim2) * 
			    uu_dim1], dabs(r__1)) > unoise) {
			++numbad;
		    }
/* L40: */
		}
		s_wsfe(&io___96);
		i__3 = *nphi;
		for (j = 1; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&ratv[j - 1], (ftnlen)sizeof(real));
		}
		e_wsfe();
/* L50: */
	    }
/* L60: */
	}
    }
L100:
    if (numbad > 0) {
	s_wsfe(&io___97);
	do_fio(&c__1, " ====  ", (ftnlen)7);
	do_fio(&c__1, (char *)&numbad, (ftnlen)sizeof(integer));
	do_fio(&c__1, "  SERIOUSLY NON-UNIT RATIOS    ====", (ftnlen)35);
	e_wsfe();
    }
    return 0;
} /* prtfin_ */

/* Subroutine */ int chekdo_()
{
    return 0;
} /* chekdo_ */

/*       Correct answers to test problems ( as produced by DISORT */
/*       running in 14-digit precision on a Cray computer ) */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Common blocks .. */
/*     .. */
/* ********************* Test Case 1A ********************************* */
/* ********************* Test Case 1B ********************************* */
/* ********************* Test Case 1C ********************************* */
/* ********************* Test Case 1D ********************************* */
/* ********************* Test Case 1E ********************************* */
/* ********************* Test Case 1F ********************************* */
/* ********************* Test Case 2A ********************************* */
/* ********************* Test Case 2B ********************************* */
/* ********************* Test Case 2C ********************************* */
/* ********************* Test Case 2D ********************************* */
/* ********************* Test Case 3A ********************************* */
/* ********************* Test Case 3B ********************************* */
/* ********************* Test Case 3C ********************************* */
/* ********************* Test Case 3D ********************************* */
/* ********************* Test Case 3E ********************************* */
/* ********************* Test Case 3F ********************************* */
/* ********************* Test Case 4A ********************************* */
/* ********************* Test Case 4B ********************************* */
/* ********************* Test Case 4C ********************************* */
/* ********************* Test Case 5A ********************************* */
/* ********************* Test Case 5B ********************************* */
/* ********************* Test Case 6A ********************************* */
/* ********************* Test Case 6B ********************************* */
/* ********************* Test Case 6C ********************************* */
/* ********************* Test Case 6D ********************************* */
/* ********************* Test Case 6E ********************************* */
/* ********************* Test Case 6F ********************************* */
/* ********************* Test Case 6G ********************************* */
/* ********************* Test Case 6H ********************************* */
/* ********************* Test Case 7A ********************************* */
/* ********************* Test Case 7B ********************************* */
/* ********************* Test Case 7C ********************************* */
/* ********************* Test Case 8A ********************************* */
/* ********************* Test Case 8B ********************************* */
/* ********************* Test Case 8C ********************************* */
/* ********************* Test Case 9A ********************************* */
/* ********************* Test Case 9B ********************************* */
/* ********************* Test Case 9C ********************************* */

/* Main program alias */ int testdo_ () { MAIN__ (); return 0; }
#ifdef __cplusplus
	}
#endif
