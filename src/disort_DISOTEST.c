/* DISOTEST.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Common Block Declarations */

struct dochek_1_ {
    doublereal tstfir[360]	/* was [5][8][9] */, tstfdn[360]	/* 
	    was [5][8][9] */, tstfup[360]	/* was [5][8][9] */, tstdfd[
	    360]	/* was [5][8][9] */, tstuu[10800]	/* was [5][10]
	    [3][8][9] */;
};

#define dochek_1 (*(struct dochek_1_ *) &dochek_)

/* Initialized data */

struct {
    doublereal e_1[2];
    doublereal fill_2[3];
    doublereal e_3[2];
    doublereal fill_4[3];
    doublereal e_5[2];
    doublereal fill_6[3];
    doublereal e_7[2];
    doublereal fill_8[3];
    doublereal e_9[2];
    doublereal fill_10[3];
    doublereal e_11[2];
    doublereal fill_12[13];
    doublereal e_13[2];
    doublereal fill_14[3];
    doublereal e_15[2];
    doublereal fill_16[3];
    doublereal e_17[2];
    doublereal fill_18[3];
    doublereal e_19[2];
    doublereal fill_20[23];
    doublereal e_21[2];
    doublereal fill_22[3];
    doublereal e_23[2];
    doublereal fill_24[3];
    doublereal e_25[2];
    doublereal fill_26[3];
    doublereal e_27[2];
    doublereal fill_28[3];
    doublereal e_29[2];
    doublereal fill_30[3];
    doublereal e_31[2];
    doublereal fill_32[13];
    doublereal e_33[3];
    doublereal fill_34[2];
    doublereal e_35[3];
    doublereal fill_36[2];
    doublereal e_37[3];
    doublereal fill_38[27];
    doublereal e_39[3];
    doublereal fill_40[2];
    doublereal e_41[3];
    doublereal fill_42[32];
    doublereal e_43[2];
    doublereal fill_44[3];
    doublereal e_45[3];
    doublereal fill_46[2];
    doublereal e_47[3];
    doublereal fill_48[2];
    doublereal e_49[3];
    doublereal fill_50[2];
    doublereal e_51[3];
    doublereal fill_52[2];
    doublereal e_53[3];
    doublereal fill_54[2];
    doublereal e_55[3];
    doublereal fill_56[2];
    doublereal e_57[3];
    doublereal fill_58[2];
    doublereal e_59[3];
    doublereal fill_60[2];
    doublereal e_61[3];
    doublereal fill_62[2];
    doublereal e_63[3];
    doublereal fill_64[27];
    doublereal e_65[3];
    doublereal fill_66[2];
    doublereal e_67[3];
    doublereal fill_68[2];
    doublereal e_69[3];
    doublereal fill_70[27];
    doublereal e_71[15];
    doublereal fill_72[25];
    doublereal e_73[2];
    doublereal fill_74[3];
    doublereal e_75[2];
    doublereal fill_76[3];
    doublereal e_77[2];
    doublereal fill_78[3];
    doublereal e_79[2];
    doublereal fill_80[3];
    doublereal e_81[2];
    doublereal fill_82[3];
    doublereal e_83[2];
    doublereal fill_84[13];
    doublereal e_85[2];
    doublereal fill_86[3];
    doublereal e_87[2];
    doublereal fill_88[3];
    doublereal e_89[2];
    doublereal fill_90[3];
    doublereal e_91[2];
    doublereal fill_92[23];
    doublereal e_93[2];
    doublereal fill_94[3];
    doublereal e_95[2];
    doublereal fill_96[3];
    doublereal e_97[2];
    doublereal fill_98[3];
    doublereal e_99[2];
    doublereal fill_100[3];
    doublereal e_101[2];
    doublereal fill_102[3];
    doublereal e_103[2];
    doublereal fill_104[13];
    doublereal e_105[3];
    doublereal fill_106[2];
    doublereal e_107[3];
    doublereal fill_108[2];
    doublereal e_109[3];
    doublereal fill_110[27];
    doublereal e_111[3];
    doublereal fill_112[2];
    doublereal e_113[3];
    doublereal fill_114[32];
    doublereal e_115[2];
    doublereal fill_116[3];
    doublereal e_117[3];
    doublereal fill_118[2];
    doublereal e_119[3];
    doublereal fill_120[2];
    doublereal e_121[3];
    doublereal fill_122[2];
    doublereal e_123[3];
    doublereal fill_124[2];
    doublereal e_125[3];
    doublereal fill_126[2];
    doublereal e_127[3];
    doublereal fill_128[2];
    doublereal e_129[3];
    doublereal fill_130[2];
    doublereal e_131[3];
    doublereal fill_132[2];
    doublereal e_133[3];
    doublereal fill_134[2];
    doublereal e_135[3];
    doublereal fill_136[27];
    doublereal e_137[3];
    doublereal fill_138[2];
    doublereal e_139[3];
    doublereal fill_140[2];
    doublereal e_141[3];
    doublereal fill_142[27];
    doublereal e_143[15];
    doublereal fill_144[25];
    doublereal e_145[2];
    doublereal fill_146[3];
    doublereal e_147[2];
    doublereal fill_148[3];
    doublereal e_149[2];
    doublereal fill_150[3];
    doublereal e_151[2];
    doublereal fill_152[3];
    doublereal e_153[2];
    doublereal fill_154[3];
    doublereal e_155[2];
    doublereal fill_156[13];
    doublereal e_157[2];
    doublereal fill_158[3];
    doublereal e_159[2];
    doublereal fill_160[3];
    doublereal e_161[2];
    doublereal fill_162[3];
    doublereal e_163[2];
    doublereal fill_164[23];
    doublereal e_165[2];
    doublereal fill_166[3];
    doublereal e_167[2];
    doublereal fill_168[3];
    doublereal e_169[2];
    doublereal fill_170[3];
    doublereal e_171[2];
    doublereal fill_172[3];
    doublereal e_173[2];
    doublereal fill_174[3];
    doublereal e_175[2];
    doublereal fill_176[13];
    doublereal e_177[3];
    doublereal fill_178[2];
    doublereal e_179[3];
    doublereal fill_180[2];
    doublereal e_181[3];
    doublereal fill_182[27];
    doublereal e_183[3];
    doublereal fill_184[2];
    doublereal e_185[3];
    doublereal fill_186[32];
    doublereal e_187[2];
    doublereal fill_188[3];
    doublereal e_189[3];
    doublereal fill_190[2];
    doublereal e_191[3];
    doublereal fill_192[2];
    doublereal e_193[3];
    doublereal fill_194[2];
    doublereal e_195[3];
    doublereal fill_196[2];
    doublereal e_197[3];
    doublereal fill_198[2];
    doublereal e_199[3];
    doublereal fill_200[2];
    doublereal e_201[3];
    doublereal fill_202[2];
    doublereal e_203[3];
    doublereal fill_204[2];
    doublereal e_205[3];
    doublereal fill_206[2];
    doublereal e_207[3];
    doublereal fill_208[27];
    doublereal e_209[3];
    doublereal fill_210[2];
    doublereal e_211[3];
    doublereal fill_212[2];
    doublereal e_213[3];
    doublereal fill_214[27];
    doublereal e_215[15];
    doublereal fill_216[25];
    doublereal e_217[2];
    doublereal fill_218[3];
    doublereal e_219[2];
    doublereal fill_220[3];
    doublereal e_221[2];
    doublereal fill_222[3];
    doublereal e_223[2];
    doublereal fill_224[3];
    doublereal e_225[2];
    doublereal fill_226[3];
    doublereal e_227[2];
    doublereal fill_228[13];
    doublereal e_229[2];
    doublereal fill_230[3];
    doublereal e_231[2];
    doublereal fill_232[3];
    doublereal e_233[2];
    doublereal fill_234[3];
    doublereal e_235[2];
    doublereal fill_236[23];
    doublereal e_237[2];
    doublereal fill_238[3];
    doublereal e_239[2];
    doublereal fill_240[3];
    doublereal e_241[2];
    doublereal fill_242[3];
    doublereal e_243[2];
    doublereal fill_244[3];
    doublereal e_245[2];
    doublereal fill_246[3];
    doublereal e_247[2];
    doublereal fill_248[13];
    doublereal e_249[3];
    doublereal fill_250[2];
    doublereal e_251[3];
    doublereal fill_252[2];
    doublereal e_253[3];
    doublereal fill_254[27];
    doublereal e_255[3];
    doublereal fill_256[2];
    doublereal e_257[3];
    doublereal fill_258[32];
    doublereal e_259[2];
    doublereal fill_260[3];
    doublereal e_261[3];
    doublereal fill_262[2];
    doublereal e_263[3];
    doublereal fill_264[2];
    doublereal e_265[3];
    doublereal fill_266[2];
    doublereal e_267[3];
    doublereal fill_268[2];
    doublereal e_269[3];
    doublereal fill_270[2];
    doublereal e_271[3];
    doublereal fill_272[2];
    doublereal e_273[3];
    doublereal fill_274[2];
    doublereal e_275[3];
    doublereal fill_276[2];
    doublereal e_277[3];
    doublereal fill_278[2];
    doublereal e_279[3];
    doublereal fill_280[27];
    doublereal e_281[3];
    doublereal fill_282[2];
    doublereal e_283[3];
    doublereal fill_284[2];
    doublereal e_285[3];
    doublereal fill_286[27];
    doublereal e_287[15];
    doublereal fill_288[25];
    doublereal e_289[2];
    doublereal fill_290[3];
    doublereal e_291[2];
    doublereal fill_292[3];
    doublereal e_293[2];
    doublereal fill_294[3];
    doublereal e_295[2];
    doublereal fill_296[3];
    doublereal e_297[2];
    doublereal fill_298[3];
    doublereal e_299[2];
    doublereal fill_300[123];
    doublereal e_301[2];
    doublereal fill_302[3];
    doublereal e_303[2];
    doublereal fill_304[3];
    doublereal e_305[2];
    doublereal fill_306[3];
    doublereal e_307[2];
    doublereal fill_308[3];
    doublereal e_309[2];
    doublereal fill_310[3];
    doublereal e_311[2];
    doublereal fill_312[123];
    doublereal e_313[2];
    doublereal fill_314[3];
    doublereal e_315[2];
    doublereal fill_316[3];
    doublereal e_317[2];
    doublereal fill_318[3];
    doublereal e_319[2];
    doublereal fill_320[3];
    doublereal e_321[2];
    doublereal fill_322[3];
    doublereal e_323[2];
    doublereal fill_324[123];
    doublereal e_325[2];
    doublereal fill_326[3];
    doublereal e_327[2];
    doublereal fill_328[3];
    doublereal e_329[2];
    doublereal fill_330[3];
    doublereal e_331[2];
    doublereal fill_332[3];
    doublereal e_333[2];
    doublereal fill_334[3];
    doublereal e_335[2];
    doublereal fill_336[123];
    doublereal e_337[2];
    doublereal fill_338[3];
    doublereal e_339[2];
    doublereal fill_340[3];
    doublereal e_341[2];
    doublereal fill_342[3];
    doublereal e_343[2];
    doublereal fill_344[3];
    doublereal e_345[2];
    doublereal fill_346[3];
    doublereal e_347[2];
    doublereal fill_348[123];
    doublereal e_349[2];
    doublereal fill_350[3];
    doublereal e_351[2];
    doublereal fill_352[3];
    doublereal e_353[2];
    doublereal fill_354[3];
    doublereal e_355[2];
    doublereal fill_356[3];
    doublereal e_357[2];
    doublereal fill_358[3];
    doublereal e_359[2];
    doublereal fill_360[423];
    doublereal e_361[2];
    doublereal fill_362[3];
    doublereal e_363[2];
    doublereal fill_364[3];
    doublereal e_365[2];
    doublereal fill_366[3];
    doublereal e_367[2];
    doublereal fill_368[3];
    doublereal e_369[2];
    doublereal fill_370[3];
    doublereal e_371[2];
    doublereal fill_372[123];
    doublereal e_373[2];
    doublereal fill_374[3];
    doublereal e_375[2];
    doublereal fill_376[3];
    doublereal e_377[2];
    doublereal fill_378[3];
    doublereal e_379[2];
    doublereal fill_380[3];
    doublereal e_381[2];
    doublereal fill_382[3];
    doublereal e_383[2];
    doublereal fill_384[123];
    doublereal e_385[2];
    doublereal fill_386[3];
    doublereal e_387[2];
    doublereal fill_388[3];
    doublereal e_389[2];
    doublereal fill_390[3];
    doublereal e_391[2];
    doublereal fill_392[3];
    doublereal e_393[2];
    doublereal fill_394[3];
    doublereal e_395[2];
    doublereal fill_396[123];
    doublereal e_397[2];
    doublereal fill_398[3];
    doublereal e_399[2];
    doublereal fill_400[3];
    doublereal e_401[2];
    doublereal fill_402[3];
    doublereal e_403[2];
    doublereal fill_404[3];
    doublereal e_405[2];
    doublereal fill_406[3];
    doublereal e_407[2];
    doublereal fill_408[723];
    doublereal e_409[2];
    doublereal fill_410[3];
    doublereal e_411[2];
    doublereal fill_412[3];
    doublereal e_413[2];
    doublereal fill_414[3];
    doublereal e_415[2];
    doublereal fill_416[3];
    doublereal e_417[2];
    doublereal fill_418[3];
    doublereal e_419[2];
    doublereal fill_420[123];
    doublereal e_421[2];
    doublereal fill_422[3];
    doublereal e_423[2];
    doublereal fill_424[3];
    doublereal e_425[2];
    doublereal fill_426[3];
    doublereal e_427[2];
    doublereal fill_428[3];
    doublereal e_429[2];
    doublereal fill_430[3];
    doublereal e_431[2];
    doublereal fill_432[123];
    doublereal e_433[2];
    doublereal fill_434[3];
    doublereal e_435[2];
    doublereal fill_436[3];
    doublereal e_437[2];
    doublereal fill_438[3];
    doublereal e_439[2];
    doublereal fill_440[3];
    doublereal e_441[2];
    doublereal fill_442[3];
    doublereal e_443[2];
    doublereal fill_444[123];
    doublereal e_445[2];
    doublereal fill_446[3];
    doublereal e_447[2];
    doublereal fill_448[3];
    doublereal e_449[2];
    doublereal fill_450[3];
    doublereal e_451[2];
    doublereal fill_452[3];
    doublereal e_453[2];
    doublereal fill_454[3];
    doublereal e_455[2];
    doublereal fill_456[123];
    doublereal e_457[2];
    doublereal fill_458[3];
    doublereal e_459[2];
    doublereal fill_460[3];
    doublereal e_461[2];
    doublereal fill_462[3];
    doublereal e_463[2];
    doublereal fill_464[3];
    doublereal e_465[2];
    doublereal fill_466[3];
    doublereal e_467[2];
    doublereal fill_468[123];
    doublereal e_469[2];
    doublereal fill_470[3];
    doublereal e_471[2];
    doublereal fill_472[3];
    doublereal e_473[2];
    doublereal fill_474[3];
    doublereal e_475[2];
    doublereal fill_476[3];
    doublereal e_477[2];
    doublereal fill_478[3];
    doublereal e_479[2];
    doublereal fill_480[423];
    doublereal e_481[3];
    doublereal fill_482[2];
    doublereal e_483[3];
    doublereal fill_484[2];
    doublereal e_485[3];
    doublereal fill_486[2];
    doublereal e_487[3];
    doublereal fill_488[2];
    doublereal e_489[3];
    doublereal fill_490[2];
    doublereal e_491[3];
    doublereal fill_492[122];
    doublereal e_493[3];
    doublereal fill_494[2];
    doublereal e_495[3];
    doublereal fill_496[2];
    doublereal e_497[3];
    doublereal fill_498[2];
    doublereal e_499[3];
    doublereal fill_500[2];
    doublereal e_501[3];
    doublereal fill_502[2];
    doublereal e_503[3];
    doublereal fill_504[122];
    doublereal e_505[3];
    doublereal fill_506[2];
    doublereal e_507[3];
    doublereal fill_508[2];
    doublereal e_509[3];
    doublereal fill_510[2];
    doublereal e_511[3];
    doublereal fill_512[2];
    doublereal e_513[3];
    doublereal fill_514[2];
    doublereal e_515[3];
    doublereal fill_516[22];
    doublereal e_517[3];
    doublereal fill_518[2];
    doublereal e_519[3];
    doublereal fill_520[2];
    doublereal e_521[3];
    doublereal fill_522[2];
    doublereal e_523[3];
    doublereal fill_524[2];
    doublereal e_525[3];
    doublereal fill_526[2];
    doublereal e_527[3];
    doublereal fill_528[22];
    doublereal e_529[3];
    doublereal fill_530[2];
    doublereal e_531[3];
    doublereal fill_532[2];
    doublereal e_533[3];
    doublereal fill_534[2];
    doublereal e_535[3];
    doublereal fill_536[2];
    doublereal e_537[3];
    doublereal fill_538[2];
    doublereal e_539[3];
    doublereal fill_540[772];
    doublereal e_541[3];
    doublereal fill_542[2];
    doublereal e_543[3];
    doublereal fill_544[2];
    doublereal e_545[3];
    doublereal fill_546[2];
    doublereal e_547[3];
    doublereal fill_548[2];
    doublereal e_549[3];
    doublereal fill_550[2];
    doublereal e_551[3];
    doublereal fill_552[122];
    doublereal e_553[3];
    doublereal fill_554[2];
    doublereal e_555[3];
    doublereal fill_556[2];
    doublereal e_557[3];
    doublereal fill_558[2];
    doublereal e_559[3];
    doublereal fill_560[2];
    doublereal e_561[3];
    doublereal fill_562[2];
    doublereal e_563[3];
    doublereal fill_564[1022];
    doublereal e_565[2];
    doublereal fill_566[3];
    doublereal e_567[2];
    doublereal fill_568[3];
    doublereal e_569[2];
    doublereal fill_570[3];
    doublereal e_571[2];
    doublereal fill_572[133];
    doublereal e_573[3];
    doublereal fill_574[2];
    doublereal e_575[3];
    doublereal fill_576[2];
    doublereal e_577[3];
    doublereal fill_578[2];
    doublereal e_579[3];
    doublereal fill_580[132];
    doublereal e_581[3];
    doublereal fill_582[2];
    doublereal e_583[3];
    doublereal fill_584[2];
    doublereal e_585[3];
    doublereal fill_586[2];
    doublereal e_587[3];
    doublereal fill_588[132];
    doublereal e_589[3];
    doublereal fill_590[2];
    doublereal e_591[3];
    doublereal fill_592[2];
    doublereal e_593[3];
    doublereal fill_594[2];
    doublereal e_595[3];
    doublereal fill_596[132];
    doublereal e_597[3];
    doublereal fill_598[2];
    doublereal e_599[3];
    doublereal fill_600[2];
    doublereal e_601[3];
    doublereal fill_602[2];
    doublereal e_603[3];
    doublereal fill_604[132];
    doublereal e_605[3];
    doublereal fill_606[2];
    doublereal e_607[3];
    doublereal fill_608[2];
    doublereal e_609[3];
    doublereal fill_610[2];
    doublereal e_611[3];
    doublereal fill_612[132];
    doublereal e_613[3];
    doublereal fill_614[2];
    doublereal e_615[3];
    doublereal fill_616[2];
    doublereal e_617[3];
    doublereal fill_618[2];
    doublereal e_619[3];
    doublereal fill_620[132];
    doublereal e_621[3];
    doublereal fill_622[2];
    doublereal e_623[3];
    doublereal fill_624[2];
    doublereal e_625[3];
    doublereal fill_626[2];
    doublereal e_627[3];
    doublereal fill_628[132];
    doublereal e_629[3];
    doublereal fill_630[2];
    doublereal e_631[3];
    doublereal fill_632[2];
    doublereal e_633[3];
    doublereal fill_634[2];
    doublereal e_635[3];
    doublereal fill_636[32];
    doublereal e_637[3];
    doublereal fill_638[2];
    doublereal e_639[3];
    doublereal fill_640[2];
    doublereal e_641[3];
    doublereal fill_642[2];
    doublereal e_643[3];
    doublereal fill_644[82];
    doublereal e_645[3];
    doublereal fill_646[2];
    doublereal e_647[3];
    doublereal fill_648[2];
    doublereal e_649[3];
    doublereal fill_650[2];
    doublereal e_651[3];
    doublereal fill_652[32];
    doublereal e_653[3];
    doublereal fill_654[2];
    doublereal e_655[3];
    doublereal fill_656[2];
    doublereal e_657[3];
    doublereal fill_658[2];
    doublereal e_659[3];
    doublereal fill_660[82];
    doublereal e_661[3];
    doublereal fill_662[2];
    doublereal e_663[3];
    doublereal fill_664[2];
    doublereal e_665[3];
    doublereal fill_666[2];
    doublereal e_667[3];
    doublereal fill_668[32];
    doublereal e_669[3];
    doublereal fill_670[2];
    doublereal e_671[3];
    doublereal fill_672[2];
    doublereal e_673[3];
    doublereal fill_674[2];
    doublereal e_675[3];
    doublereal fill_676[832];
    doublereal e_677[3];
    doublereal fill_678[2];
    doublereal e_679[3];
    doublereal fill_680[2];
    doublereal e_681[3];
    doublereal fill_682[2];
    doublereal e_683[3];
    doublereal fill_684[132];
    doublereal e_685[3];
    doublereal fill_686[2];
    doublereal e_687[3];
    doublereal fill_688[2];
    doublereal e_689[3];
    doublereal fill_690[2];
    doublereal e_691[3];
    doublereal fill_692[132];
    doublereal e_693[3];
    doublereal fill_694[2];
    doublereal e_695[3];
    doublereal fill_696[2];
    doublereal e_697[3];
    doublereal fill_698[2];
    doublereal e_699[3];
    doublereal fill_700[882];
    doublereal e_701[20];
    doublereal fill_702[130];
    doublereal e_703[20];
    doublereal fill_704[130];
    doublereal e_705[20];
    doublereal fill_706[30];
    doublereal e_707[20];
    doublereal fill_708[30];
    doublereal e_709[20];
    doublereal fill_710[780];
    } dochek_ = { 3.14159, 2.29844, {0}, 3.14159, 2.29844, {0}, 0., 0., {0}, 
	    3.14159, 0., {0}, 3.14159, 0., {0}, 0., 0., {0}, .252716, 
	    .0210311, {0}, .252716, .0210311, {0}, .252716, 2.56077e-28, {0}, 
	    .252716, 2.56077e-28, {0}, 3.14159, 1.42628e-4, {0}, 3.14159, 
	    1.42628e-4, {0}, 0., 0., {0}, 3.14159, 5.67011e-35, {0}, 3.14159, 
	    5.67011e-35, {0}, 0., 0., {0}, 3.14159, 1.90547, 1.15573, {0}, 
	    3.14159, 1.90547, 1.15573, {0}, 1.5708, .577864, .212584, {0}, 
	    3.14159, 3.97856e-14, 5.03852e-28, {0}, .128058, 8.67322e-6, 
	    4.47729e-21, {0}, 100., 100., {0}, 100., 36.7879, 13.5335, {0}, 
	    100., 36.7879, 13.5335, {0}, 100., 36.7879, 13.5335, {0}, 100., 
	    36.7879, 13.5335, {0}, 100., 36.7879, 13.5335, {0}, 100., 36.7879,
	     13.5335, {0}, 100., 13.5335, 2.06115e-7, {0}, 100., 36.7879, 
	    13.5335, {0}, 100., 36.7879, 13.5335, {0}, 100., 36.7879, 13.5335,
	     {0}, 0., 0., 0., {0}, 0., 0., 0., {0}, 0., 0., 0., {0}, 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 1.5708, .192354, .023555, 
	    9.65131e-6, 9.03133e-19, {0}, 0., .0794108, {0}, 0., .420233, {0},
	     3.14159, 3.04897, {0}, 0., 0., {0}, 0., .0676954, {0}, 3.14159, 
	    .00460048, {0}, 0., .0441791, {0}, 0., .106123, {0}, 0., 
	    2.51683e-4, {0}, 0., .0268008, {0}, 0., .0470303, {0}, 0., 
	    1.31455, {0}, 3.14159, 2.4966, {0}, 0., 1.99813e-5, {0}, 0., 
	    .591973, {0}, 3.14159, 1.01167, {0}, 0., 1.17401, 1.81264, {0}, 
	    0., 1.01517, 1.51554, {0}, 0., .702764, .803294, {0}, 0., 2.24768,
	     .479851, {0}, 1.74767, .233975, 6.38347e-5, {0}, 0., 0., {0}, 0.,
	     0., 0., {0}, 0., 0., 0., {0}, 0., 0., 0., {0}, 0., 0., 0., {0}, 
	    321.497, 142.493, 70.5306, {0}, 321.497, 304.775, 363.632, {0}, 
	    321.497, 255.455, 443.444, {0}, 319.83, 354.099, 301.334, {0}, 
	    319.83, 350.555, 292.063, {0}, 319.83, 353.25, 298.583, {0}, 1., 
	    .722235, .513132, {0}, 1., .795332, .650417, {0}, 1., .486157, 
	    .159984, {0}, 1., .355151, .144265, .00671445, 6.16968e-7, 1., 
	    .452357, .236473, .0276475, 7.41854e-5, 6.09217, 4.97279, 4.46616,
	     4.22731, 4.73767, {0}, .0799451, 0., {0}, .422922, 0., {0}, 
	    .0906556, 0., {0}, .259686, 0., {0}, 3.0739, 0., {0}, 2.49618, 0.,
	     {0}, .0535063, 0., {0}, .125561, 0., {0}, .062473, 0., {0}, 
	    .225915, 0., {0}, .169951, 0., {0}, 1.8269, 0., {0}, .583122, 0., 
	    {0}, .170062, 0., {0}, 2.54962, 0., {0}, 1.6875, 0., {0}, .173223,
	     .111113, 0., {0}, .123666, .0788691, 0., {0}, .225487, .123848, 
	    0., {0}, 2.66174, 1.76783, 0., {0}, .270485, .0374253, 1.02904e-5,
	     {0}, 0., 0., {0}, 0., 0., 0., {0}, 1.4845, 2.99914, 6.76676, {0},
	     .619105, 1.33245, 3.4779, {0}, 82.8836, 165.215, 359.895, {0}, 
	    85.2994, 170.343, 372.842, {0}, 336.115, 413.977, 443.05, {0}, 
	    237.35, 261.13, 456.205, {0}, 429.572, 447.018, 594.576, {0}, 
	    312.563, 268.126, 305.596, {0}, 407.04, 411.058, 530.504, {0}, 
	    .0929634, .0278952, 0., {0}, .225136, .126349, 0., {0}, .378578, 
	    .243397, 0., {0}, .227973, .0875098, .0361819, .00219291, 0., 
	    .100079, .0452015, .0241941, .00416017, 0., 4.68414, 4.24381, 
	    4.16941, 4.30667, 5.11524, {0}, 25.4067, 18.6531, {0}, 0., 0., {0}
	    , .066687, .0588936, {0}, 25.7766, 0., {0}, 0., 0., {0}, .114239, 
	    7.93633e-5, {0}, 1.6657, .189848, {0}, 0., 0., {0}, 1.67462, 
	    1.75464e-4, {0}, 0., 0., {0}, 25.9221, .0743638, {0}, 0., 0., {0},
	     .0806945, .0435243, {0}, 25.9222, 1.91941e-5, {0}, 0., 0., {0}, 
	    .100459, .0171385, {0}, 0., 0., 0., {0}, .343724, .35239, .31945, 
	    {0}, .385003, .337317, .216403, {0}, 0., 0., 0., {0}, .310129, 
	    .0452671, 1.25022e-5, {0}, 200., 200., {0}, 200., 73.5759, 
	    27.0671, {0}, 202.01, 77.9962, 40.6006, {0}, 200.899, 75.7612, 
	    36.5186, {0}, 310.552, 311.018, 679.533, {0}, 957.001, 529.253, 
	    807.923, {0}, 581.342, 128.621, -171.119, {0}, 423.78, 61.9828, 
	    -31.7719, {0}, -80.427, 251.589, 715.964, {0}, -168.356, 101.251, 
	    409.326, {0}, -98.5693, 217.724, 623.936, {0}, 1.12474, .651821, 
	    .563361, {0}, .512692, .356655, .0568095, {0}, .565095, .276697, 
	    .0135679, {0}, .882116, .232366, .0933443, .00392782, 1.025e-7, 
	    .804577, .25533, .130976, .0136227, 1.22022e-5, 3.49563, .881206, 
	    .350053, .0193471, .0715349, {0}, 0., .0133826, {0}, 0., .0263324,
	     {0}, 0., .115898, {0}, .117771, 0., {0}, .026417, 0., {0}, 
	    .0134041, 0., {0}, 0., .0708109, {0}, 0., .139337, {0}, 0., 
	    .613458, {0}, .622884, 0., {0}, .139763, 0., {0}, .0709192, 0., {
	    0}, 1., .984447, {0}, 1., .969363, {0}, 1., .863946, {0}, .133177,
	     0., {0}, .0299879, 0., {0}, .0152233, 0., {0}, 0., 1.2298e-15, {
	    0}, 0., 1.30698e-17, {0}, 0., 6.88841e-18, {0}, .262972, 0., {0}, 
	    .0906967, 0., {0}, .0502853, 0., {0}, 0., .0271316, {0}, 0., 
	    .0187805, {0}, 0., .0116385, {0}, 1.93321, 0., {0}, 1.02732, 0., {
	    0}, .797199, 0., {0}, 1., .0018684, {0}, 1., .00126492, {0}, 1., 
	    7.79281e-4, {0}, .87751, 0., {0}, .815136, 0., {0}, .752715, 0., {
	    0}, 0., .00771897, {0}, 0., .0200778, {0}, 0., .0257685, {0}, 
	    .161796, 0., {0}, .0211501, 0., {0}, .00786713, 0., {0}, 0., 
	    .0186027, {0}, 0., .0464061, {0}, 0., .0677603, {0}, .347678, 0., 
	    {0}, .048712, 0., {0}, .0189387, 0., {0}, 0., 1.70004e-4, {0}, 0.,
	     3.97168e-5, {0}, 0., 1.32472e-5, {0}, .162566, 0., {0}, .0245786,
	     0., {0}, .0101498, 0., {0}, 0., .010595, {0}, 0., .00769337, {0},
	     0., .00379276, {0}, .36401, 0., {0}, .0826993, 0., {0}, .049237, 
	    0., {0}, 0., .00788518, {0}, 0., .0222081, {0}, 0., .00290003, {0}
	    , .486923, 0., {0}, .0508802, 0., {0}, .0103351, 0., {0}, 0., 
	    .213998, {0}, 0., .58401, {0}, 0., .377162, {0}, 3.72886, 0., {0},
	     .632069, 0., {0}, .151211, 0., {0}, 1., .911175, {0}, 1., 
	    .744031, {0}, 1., .403607, {0}, .566618, 0., {0}, .232213, 0., {0}
	    , .0760274, 0., {0}, 0., 2.07944e-5, {0}, 0., 9.76964e-7, {0}, 0.,
	     1.88004e-7, {0}, .486926, 0., {0}, .050919, 0., {0}, .0103691, 
	    0., {0}, 0., .236362, {0}, 0., .164357, {0}, 0., .0915016, {0}, 
	    3.85776, 0., {0}, .868015, 0., {0}, .379793, 0., {0}, 1., .416334,
	     {0}, 1., .275956, {0}, 1., .151785, {0}, .757316, 0., {0}, 
	    .587262, 0., {0}, .433051, 0., {0}, 0., 2.49435, 3.34666, {0}, 0.,
	     .11907, .219624, {0}, 0., .135138, .157002, {0}, .0929865, 
	    .12407, 0., {0}, .0661234, .0402884, 0., {0}, .0355365, .0173581, 
	    0., {0}, 0., 2.22302, 2.94684, {0}, 0., .0964101, .167508, {0}, 
	    0., .0962927, .108212, {0}, .0655782, .0844926, 0., {0}, .0456643,
	     .0280216, 0., {0}, .0274243, .0135087, 0., {0}, 0., .0476042, 
	    .0837538, {0}, 0., 3.00258, 2.68792, {0}, 0., 1.4114, .876317, {0}
	    , .870016, .6974, 0., {0}, .22476, .109065, 0., {0}, .0228322, 
	    .00935116, 0., {0}, 0., .0476042, .0837538, {0}, 0., .0582549, 
	    .0943348, {0}, 0., .104636, .0895961, {0}, .0886109, .0915334, 0.,
	     {0}, .0576913, .0295681, 0., {0}, .0228322, .00935116, 0., {0}, 
	    0., .0476042, .0837538, {0}, 0., .0259416, .0400024, {0}, 0., 
	    .0624961, .0466785, {0}, .0695291, .0590188, 0., {0}, .0493283, 
	    .0244592, 0., {0}, .0228322, .00935116, 0., {0}, 0., .753662, 
	    .19523, {0}, 0., .696362, .13199, {0}, 0., .650541, .0720655, {0},
	     .45766, .627631, 0., {0}, .768942, .581809, 0., {0}, 1.03122, 
	    .524532, 0., {0}, 15.0205, .18509, 3.41287e-5, {0}, .218778, 
	    .0492726, 1.39916e-5, {0}, .13652, .0265448, 7.47042e-6, {0}, 
	    .114001, .0202154, 5.65604e-6, {0}, .0871231, .012966, 3.58246e-6,
	     {0}, .0855016, .00951226, 2.57859e-6, {0}, 0., 0., {0}, 0., 0., {
	    0}, 0., 0., {0}, 0., 0., {0}, 0., 0., 0., {0}, 0., 0., 0., {0}, 
	    0., 0., 0., {0}, 0., 0., 0., {0}, 0., 0., 0., {0}, 0., 0., 0., {0}
	    , 9.77882e-5, .0145131, 2.15393, {0}, .792386, 1.30642, 2.15393, {
	    0}, 0., 0., 0., {0}, 0., 0., 0., {0}, 5.07347e-5, .0075297, 
	    1.11751, {0}, .272396, .449105, .740449, {0}, 0., 0., 0., {0}, 0.,
	     0., 0., {0}, .00320501, .475666, 70.5951, {0}, 46.4836, 76.6386, 
	    126.356, {0}, 102.336, 62.0697, 37.6472, {0}, 102.336, .689532, 
	    .00464603, {0}, .00358884, .532631, 79.0494, {0}, 47.5119, 
	    78.3338, 129.151, {0}, 102.336, 97.8748, 110.061, {0}, 102.336, 
	    101.048, 138.631, {0}, 78.0731, 115.783, 133.373, {0}, 117.126, 
	    136.113, 142.865, {0}, 102.336, 84.9992, 138.631, {0}, 102.336, 
	    77.3186, 145.441, {0}, 71.2616, 78.831, 144.088, {0}, 78.0737, 
	    85.6424, 145.547, {0}, 101.805, 106.346, 96.416, {0}, 101.805, 
	    134.218, 88.7832, {0}, 141.637, 102.514, 189.259, {0}, 148.416, 
	    158.864, 189.259, {0}, 101.805, 106.346, 96.416, {0}, 101.805, 
	    106.238, 74.8649, {0}, 127.778, 94.1347, 189.259, {0}, 148.416, 
	    158.864, 189.259, {0}, 101.805, 105.967, 95.3651, {0}, 101.805, 
	    128.779, 74.7553, {0}, 135.84, 90.0044, 97.2743, {0}, 95.6588, 
	    88.7624, 97.2743, {0}, 101.805, 105.967, 95.3651, {0}, 101.805, 
	    100.799, 60.837, {0}, 121.98, 81.6247, 97.2743, {0}, 95.6588, 
	    88.7624, 97.2743, {0}, 101.805, 106.259, 96.1392, {0}, 101.805, 
	    132.945, 84.5884, {0}, 140.347, 99.2793, 154.901, {0}, 140.534, 
	    148.621, 176.268, {0}, 101.805, 106.259, 96.1392, {0}, 101.805, 
	    104.823, 69.7291, {0}, 126.343, 90.2686, 137.958, {0}, 140.534, 
	    148.621, 176.268, {0}, .31831, .262711, .210014, {0}, .31831, 
	    .136952, .0560376, {0}, .0562566, .0184909, 0., {0}, .0194423, 
	    .00552188, 0., {0}, .31831, .277499, .240731, {0}, .31831, .18395,
	     .129291, {0}, .123687, .0835695, 0., {0}, .0495581, .0250575, 0.,
	     {0}, .31831, .18902, .0684762, {0}, .31831, .0988158, .0296698, {
	    0}, .149335, .0965192, 0., {0}, .104766, .0654445, 0., {0}, 
	    .31831, .153507, .0706614, .00372784, 2.87656e-7, .31831, 
	    .0509531, .0209119, .00108815, 1.05921e-7, .0998915, .0367006, 
	    .0148545, 8.83317e-4, 0., .0591345, .0231903, .00972307, 
	    5.94743e-4, 0., {0}, .31831, .196609, .115478, .0146177, 
	    3.37743e-5, .31831, .0592369, .0301809, .0038559, 1.20858e-5, 
	    .0739198, .030023, .0152672, .00238301, 0., .0132768, .00705566, 
	    .00406932, 7.77891e-4, 0., {0}, 1.9392, 1.66764, 1.48511, 1.34514,
	     1.48927, 1.9392, 1.44453, 1.35009, 1.35131, 1.5427, 1.61855, 
	    1.38339, 1.33079, 1.3598, 1.62823, 1.43873, 1.3389, 1.32794, 
	    1.37918, 1.62823, {0}, 1.9392, 1.66764, 1.48511, 1.34514, 1.48927,
	     1.9392, 1.42925, 1.34587, 1.35129, 1.5427, 1.57895, 1.37317, 
	    1.32921, 1.35979, 1.62823, 1.43873, 1.3389, 1.32794, 1.37918, 
	    1.62823, {0}, 1.9392, 1.66764, 1.48511, 1.34514, 1.48927, 1.9392, 
	    1.42444, 1.34469, 1.35128, 1.5427, 1.56559, 1.37034, 1.32873, 
	    1.35979, 1.62823, 1.43873, 1.3389, 1.32794, 1.37918, 1.62823 };


/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b2 = 0.;
static integer c__6 = 6;
static integer c__5 = 5;
static integer c__10 = 10;
static integer c__48 = 48;
static integer c__3 = 3;
static logical c_false = FALSE_;
static integer c__2 = 2;
static doublereal c_b79 = .75;
static integer c__4 = 4;
static doublereal c_b209 = .7;
static doublereal c_b243 = .8;
static doublereal c_b457 = .9;
static logical c_true = TRUE_;

/* Main program */ MAIN__()
{
    /* Initialized data */

    static char abc[1*18+1] = "abcdefghijklmnopqr";
    static logical prnt[7] = { TRUE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,TRUE_ 
	    };
    static doublereal accur = 0.;
    static char blanks[3+1] = "   ";
    static logical doprob[13] = { TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,
	    TRUE_,TRUE_,TRUE_,TRUE_,TRUE_,TRUE_ };

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3, i__4, i__5, i__6[2];
    doublereal d__1;

    /* Builtin functions */
    double asin(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), s_rnge(char 
	    *, integer, char *, integer), e_wsfi(), i_indx(char *, char *, 
	    ftnlen, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double pow_di(doublereal *, integer *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_stop(char *
	    , ftnlen);

    /* Local variables */
    static integer icas;
    static doublereal dfdt[5];
    static integer nphi;
    static doublereal uavg[5], flup[5];
    static integer ntau;
    static doublereal pmom[294]	/* was [49][6] */, utau[5];
    static integer nlyr, numu, nstr;
    static doublereal fbeam;
    static integer j, k, ibcnd;
    static doublereal dtauc[6], ssalb[6];
    static logical plank;
    static doublereal rfldn[5], btemp, temis;
    static char title[100];
    static integer nprob;
    static doublereal fisot, ttemp, cmpuu[150]	/* was [5][10][3] */;
    static integer lc;
    static doublereal hl[49], albmed[10], albedo;
    static char header[127];
    static doublereal pi, cmpdfd[5];
    static integer iu, lu;
    static logical lamber, deltam;
    static doublereal cmpfdn[5], uu[150]	/* was [10][5][3] */, cmpfir[
	    5], rfldir[5];
    static logical azmavg;
    static doublereal trnmed[10], cmpfup[5], temper[7];
    extern /* Subroutine */ int getmom_(integer *, doublereal *, integer *, 
	    doublereal *);
    static logical usrang;
    static integer lentit;
    static logical onlyfl;
    extern /* Subroutine */ int disort_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , integer *, doublereal *, integer *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , logical *, logical *, doublereal *, logical *, char *, integer *
	    , integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, ftnlen), 
	    errmsg_(char *, logical *, ftnlen);
    static doublereal wvnmhi;
    extern /* Subroutine */ int prtfin_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    logical *, logical *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *);
    static doublereal wvnmlo;
    static logical usrtau;
    static doublereal u0u[50]	/* was [10][5] */;
    static integer iod;
    static doublereal phi[3];
    static integer iss;
    static doublereal umu[10], phi0, umu0;

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
    pi = 2. * asin(1.);
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
	utau[0] = 0.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = -1.;
	umu[1] = -.5;
	umu[2] = -.1;
	umu[3] = .1;
	umu[4] = .5;
	umu[5] = 1.;
	nphi = 1;
	phi[0] = 0.;
	ibcnd = 0;
	umu0 = .1;
	phi0 = 0.;
	lamber = TRUE_;
	albedo = 0.;
	deltam = FALSE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 6; ++icas) {
	    if (icas == 1) {
		utau[1] = .03125;
		ssalb[0] = .2;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 2) {
		utau[1] = .03125;
		ssalb[0] = 1.;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 3) {
		utau[1] = .03125;
		ssalb[0] = .99;
		fbeam = 0.;
		fisot = 1.;
	    } else if (icas == 4) {
		utau[1] = 32.;
		ssalb[0] = .2;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 5) {
		utau[1] = 32.;
		ssalb[0] = 1.;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 6) {
		utau[1] = 32.;
		ssalb[0] = .99;
		fbeam = 0.;
		fisot = 1.;
	    }
	    dtauc[0] = utau[1];
	    s_wsfi(&io___32);
	    do_fio(&c__1, "Test Case No. 1", (ftnlen)15);
	    do_fio(&c__1, abc + ((i__1 = icas - 1) < 18 && 0 <= i__1 ? i__1 : 
		    s_rnge("abc", i__1, "testdo_", (ftnlen)191)), (ftnlen)1);
	    do_fio(&c__1, ":  Isotropic Scattering, Ref. VH1, Table 12:  b =",
		     (ftnlen)49);
	    do_fio(&c__1, (char *)&utau[1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ", a =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ssalb[0], (ftnlen)sizeof(doublereal));
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)209)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)209)], &dochek_1.tstfup[(i__3 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tst\
fup", i__3, "testdo_", (ftnlen)209)], &dochek_1.tstdfd[(i__4 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__4 ? i__4 : s_rnge("tst\
dfd", i__4, "testdo_", (ftnlen)209)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    209)], &c__5, &c__10, &c__3);
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
	utau[0] = 0.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = -.981986;
	umu[1] = -.538263;
	umu[2] = -.018014;
	umu[3] = .018014;
	umu[4] = .538263;
	umu[5] = .981986;
	nphi = 1;
	phi[0] = 0.;
	ibcnd = 0;
	fbeam = pi;
	umu0 = .080442;
	phi0 = 0.;
	fisot = 0.;
	lamber = TRUE_;
	albedo = 0.;
	deltam = FALSE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	icas = 0;
	for (iod = 1; iod <= 2; ++iod) {
	    if (iod == 1) {
		utau[1] = .2;
	    }
	    if (iod == 2) {
		utau[1] = 5.;
	    }
	    dtauc[0] = utau[1];
	    for (iss = 1; iss <= 2; ++iss) {
		if (iss == 1) {
		    ssalb[0] = .5;
		}
		if (iss == 2) {
		    ssalb[0] = 1.;
		}
		++icas;
		s_wsfi(&io___53);
		do_fio(&c__1, "Test Case No. 2", (ftnlen)15);
		do_fio(&c__1, abc + ((i__1 = icas - 1) < 18 && 0 <= i__1 ? 
			i__1 : s_rnge("abc", i__1, "testdo_", (ftnlen)268)), (
			ftnlen)1);
		do_fio(&c__1, ", Rayleigh Scattering, Ref. SW, Table 1:  tau\
 =", (ftnlen)47);
		do_fio(&c__1, (char *)&utau[1], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, ", mu0 =", (ftnlen)7);
		do_fio(&c__1, (char *)&umu0, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, ", ss-albedo =", (ftnlen)13);
		do_fio(&c__1, (char *)&ssalb[0], (ftnlen)sizeof(doublereal));
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
			dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 
			45) < 360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1,
			 "testdo_", (ftnlen)287)], &dochek_1.tstfdn[(i__2 = (
			icas + (nprob << 3)) * 5 - 45) < 360 && 0 <= i__2 ? 
			i__2 : s_rnge("tstfdn", i__2, "testdo_", (ftnlen)287)]
			, &dochek_1.tstfup[(i__3 = (icas + (nprob << 3)) * 5 
			- 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tstfup", 
			i__3, "testdo_", (ftnlen)287)], &dochek_1.tstdfd[(
			i__4 = (icas + (nprob << 3)) * 5 - 45) < 360 && 0 <= 
			i__4 ? i__4 : s_rnge("tstdfd", i__4, "testdo_", (
			ftnlen)287)], &dochek_1.tstuu[(i__5 = (((icas + (
			nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 &&
			 0 <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (
			ftnlen)287)], &c__5, &c__10, &c__3);
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
	getmom_(&c__3, &c_b79, &nstr, pmom);
	usrtau = TRUE_;
	ntau = 2;
	utau[0] = 0.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = -1.;
	umu[1] = -.5;
	umu[2] = -.1;
	umu[3] = .1;
	umu[4] = .5;
	umu[5] = 1.;
	nphi = 1;
	phi[0] = 0.;
	ibcnd = 0;
	umu0 = .1;
	phi0 = 0.;
	deltam = TRUE_;
	lamber = TRUE_;
	onlyfl = FALSE_;
	albedo = 0.;
	plank = FALSE_;
	for (icas = 1; icas <= 6; ++icas) {
	    if (icas == 1) {
		utau[1] = 1.;
		ssalb[0] = .2;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 2) {
		utau[1] = 1.;
		ssalb[0] = 1.;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 3) {
		utau[1] = 1.;
		ssalb[0] = .99;
		fbeam = 0.;
		fisot = 1.;
	    } else if (icas == 4) {
		utau[1] = 8.;
		ssalb[0] = .2;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 5) {
		utau[1] = 8.;
		ssalb[0] = 1.;
		fbeam = pi / umu0;
		fisot = 0.;
	    } else if (icas == 6) {
		utau[1] = 8.;
		ssalb[0] = .99;
		fbeam = 0.;
		fisot = 1.;
	    }
	    dtauc[0] = utau[1];
	    s_wsfi(&io___54);
	    do_fio(&c__1, "Test Case No. 3", (ftnlen)15);
	    do_fio(&c__1, abc + ((i__1 = icas - 1) < 18 && 0 <= i__1 ? i__1 : 
		    s_rnge("abc", i__1, "testdo_", (ftnlen)378)), (ftnlen)1);
	    do_fio(&c__1, ", Henyey-Greenstein Scattering, Ref. VH2, Table 3\
5, g = 0.75, b =", (ftnlen)65);
	    do_fio(&c__1, (char *)&utau[1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ", a =", (ftnlen)5);
	    do_fio(&c__1, (char *)&ssalb[0], (ftnlen)sizeof(doublereal));
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)396)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)396)], &dochek_1.tstfup[(i__3 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tst\
fup", i__3, "testdo_", (ftnlen)396)], &dochek_1.tstdfd[(i__4 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__4 ? i__4 : s_rnge("tst\
dfd", i__4, "testdo_", (ftnlen)396)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    396)], &c__5, &c__10, &c__3);
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
	dtauc[0] = 1.;
	usrtau = TRUE_;
	ntau = 3;
	utau[0] = 0.;
	utau[1] = .5;
	utau[2] = 1.;
	usrang = TRUE_;
	numu = 6;
	umu[0] = -1.;
	umu[1] = -.5;
	umu[2] = -.1;
	umu[3] = .1;
	umu[4] = .5;
	umu[5] = 1.;
	ibcnd = 0;
	fbeam = pi;
	phi0 = 0.;
	fisot = 0.;
	lamber = TRUE_;
	albedo = 0.;
	deltam = TRUE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    s_wsfi(&io___56);
	    do_fio(&c__1, "Test Case No. 4", (ftnlen)15);
	    do_fio(&c__1, abc + ((i__1 = icas - 1) < 18 && 0 <= i__1 ? i__1 : 
		    s_rnge("abc", i__1, "testdo_", (ftnlen)443)), (ftnlen)1);
	    do_fio(&c__1, ", Haze-L Scattering, Ref. GS, Table ", (ftnlen)36);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		ssalb[0] = 1.;
		nphi = 1;
		phi[0] = 0.;
		umu0 = 1.;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 3, a__1[1] = " 12";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 2) {
		ssalb[0] = .9;
		nphi = 1;
		phi[0] = 0.;
		umu0 = 1.;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 3, a__1[1] = " 13";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 3) {
		ssalb[0] = .9;
		nphi = 3;
		phi[0] = 0.;
		phi[1] = 90.;
		phi[2] = 180.;
		umu0 = .5;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 6, a__1[1] = " 14-16";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)489)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)489)], &dochek_1.tstfup[(i__3 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tst\
fup", i__3, "testdo_", (ftnlen)489)], &dochek_1.tstdfd[(i__4 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__4 ? i__4 : s_rnge("tst\
dfd", i__4, "testdo_", (ftnlen)489)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    489)], &c__5, &c__10, &c__3);
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
	dtauc[0] = 64.;
	usrtau = TRUE_;
	ntau = 3;
	usrang = TRUE_;
	numu = 6;
	umu[0] = -1.;
	umu[1] = -.5;
	umu[2] = -.1;
	umu[3] = .1;
	umu[4] = .5;
	umu[5] = 1.;
	nphi = 1;
	phi[0] = 0.;
	ibcnd = 0;
	fbeam = pi;
	umu0 = 1.;
	phi0 = 0.;
	fisot = 0.;
	lamber = TRUE_;
	albedo = 0.;
	deltam = TRUE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    s_wsfi(&io___58);
	    do_fio(&c__1, "Test Case No. 5", (ftnlen)15);
	    do_fio(&c__1, abc + ((i__1 = icas - 1) < 18 && 0 <= i__1 ? i__1 : 
		    s_rnge("abc", i__1, "testdo_", (ftnlen)536)), (ftnlen)1);
	    do_fio(&c__1, ", Cloud C.1 Scattering, Ref. GS, Table ", (ftnlen)
		    39);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		utau[0] = 0.;
		utau[1] = 32.;
		utau[2] = 64.;
		ssalb[0] = 1.;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 3, a__1[1] = " 19";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    }
	    if (icas == 2) {
		utau[0] = 3.2;
		utau[1] = 12.8;
		utau[2] = 48.;
		ssalb[0] = .9;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 3, a__1[1] = " 20";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)574)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)574)], &dochek_1.tstfup[(i__3 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tst\
fup", i__3, "testdo_", (ftnlen)574)], &dochek_1.tstdfd[(i__4 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__4 ? i__4 : s_rnge("tst\
dfd", i__4, "testdo_", (ftnlen)574)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    574)], &c__5, &c__10, &c__3);
/* L50: */
	}
    }
    if (doprob[5]) {
/* ********************************************************************** */
/* ****  Test Problem 6:  No Scattering, Increasingly Complex Sources**** */
/* ********************************************************************** */
	nstr = 16;
	nlyr = 1;
	ssalb[0] = 0.;
	wvnmlo = 0.;
	wvnmhi = 5e4;
	usrtau = TRUE_;
	usrang = TRUE_;
	numu = 4;
	umu[0] = -1.;
	umu[1] = -.1;
	umu[2] = .1;
	umu[3] = 1.;
	nphi = 1;
	phi[0] = 90.;
	ibcnd = 0;
	fbeam = 200.;
	umu0 = .5;
	phi0 = 0.;
	fisot = 0.;
	temis = 1.;
	onlyfl = FALSE_;
	deltam = FALSE_;
	for (icas = 1; icas <= 8; ++icas) {
	    s_wsfi(&io___59);
	    do_fio(&c__1, "Test Case No. 6", (ftnlen)15);
	    do_fio(&c__1, abc + ((i__1 = icas - 1) < 18 && 0 <= i__1 ? i__1 : 
		    s_rnge("abc", i__1, "testdo_", (ftnlen)616)), (ftnlen)1);
	    do_fio(&c__1, ": No Scattering; Source = Beam", (ftnlen)30);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		ntau = 2;
		utau[0] = 0.;
		utau[1] = 0.;
	    } else if (icas > 1) {
		ntau = 3;
		utau[0] = 0.;
		utau[1] = .5;
		utau[2] = 1.;
	    }
	    if (icas == 1) {
/*                                    ** Transparent medium, beam source */
		dtauc[0] = 0.;
		lamber = TRUE_;
		albedo = 0.;
		plank = FALSE_;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 19, a__1[1] = "; Bottom Albedo = 0";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 2) {
/*                                    ** Add some optical depth */
		dtauc[0] = 1.;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 19, a__1[1] = "; Bottom Albedo = 0";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 3) {
/*                                   ** Add some isotropic reflection */
		lamber = TRUE_;
		albedo = .5;
		plank = FALSE_;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 27, a__1[1] = "; Bottom Albedo=0.5 Lambert";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 4) {
/*                                   ** Use non-isotropic reflection */
		dtauc[0] = 1.;
		lamber = FALSE_;
		i__1 = nstr;
		for (k = 0; k <= i__1; ++k) {
		    hl[(i__2 = k) < 49 && 0 <= i__2 ? i__2 : s_rnge("hl", 
			    i__2, "testdo_", (ftnlen)660)] = pow_di(&c_b209, &
			    k);
/* L59: */
		}
		plank = FALSE_;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 29, a__1[1] = "; Bottom Albedo = Non-Lambert";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 5) {
/*                                   ** Add some bottom-boundary emission */
		dtauc[0] = 1.;
		temper[0] = 0.;
		temper[1] = 0.;
		lamber = FALSE_;
		btemp = 300.;
		ttemp = 0.;
		plank = TRUE_;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 41, a__1[1] = ", Bottom Emission; Bott Alb = Non-L\
ambert";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 6) {
/*                                   ** Add some top-boundary diffuse */
/*                                      incidence (prescribed + emitted) */
		dtauc[0] = 1.;
		temper[0] = 0.;
		temper[1] = 0.;
		fisot = 100. / pi;
		lamber = FALSE_;
		btemp = 300.;
		ttemp = 250.;
		plank = TRUE_;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 45, a__1[1] = ", Bottom+Top Emission; Bott Alb = N\
on-Lambert";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 7) {
/*                                   ** Add some internal emission */
		dtauc[0] = 1.;
		temper[0] = 250.;
		temper[1] = 300.;
		lamber = FALSE_;
		btemp = 300.;
		ttemp = 250.;
		plank = TRUE_;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 54, a__1[1] = ", Bottom+Top+Internal Emission; Bot\
t Alb = Non-Lambert";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 8) {
/*                                   ** Increase the optical depth */
		dtauc[0] = 10.;
		temper[0] = 250.;
		temper[1] = 300.;
		utau[0] = 0.;
		utau[1] = 1.;
		utau[2] = 10.;
		lamber = FALSE_;
		btemp = 300.;
		ttemp = 250.;
		plank = TRUE_;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 54, a__1[1] = ", Bottom+Top+Internal Emission; Bot\
t Alb = Non-Lambert";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)734)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)734)], &dochek_1.tstfup[(i__3 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tst\
fup", i__3, "testdo_", (ftnlen)734)], &dochek_1.tstdfd[(i__4 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__4 ? i__4 : s_rnge("tst\
dfd", i__4, "testdo_", (ftnlen)734)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    734)], &c__5, &c__10, &c__3);
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
	getmom_(&c__3, &c_b243, &nstr, pmom);
	dtauc[0] = 1.;
	ssalb[0] = .5;
	temper[0] = 300.;
	temper[1] = 200.;
	wvnmlo = 0.;
	wvnmhi = 5e4;
	usrtau = TRUE_;
	ntau = 3;
	utau[0] = 0.;
	utau[1] = .5;
	utau[2] = 1.;
	usrang = TRUE_;
	numu = 4;
	umu[0] = -1.;
	umu[1] = -.1;
	umu[2] = .1;
	umu[3] = 1.;
	nphi = 2;
	phi[0] = 0.;
	phi[1] = 90.;
	ibcnd = 0;
	fbeam = 200.;
	umu0 = .5;
	phi0 = 0.;
	fisot = 100.;
	btemp = 320.;
	ttemp = 100.;
	temis = 1.;
	deltam = TRUE_;
	plank = TRUE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    s_wsfi(&io___61);
	    do_fio(&c__1, "Test Case No. 7", (ftnlen)15);
	    do_fio(&c__1, abc + ((i__1 = icas - 1) < 18 && 0 <= i__1 ? i__1 : 
		    s_rnge("abc", i__1, "testdo_", (ftnlen)789)), (ftnlen)1);
	    do_fio(&c__1, ": Absorption + Henyey-Greenstein Scattering, All \
Sources", (ftnlen)56);
	    e_wsfi();
	    lentit = i_indx(title, blanks, (ftnlen)100, (ftnlen)3);
	    if (icas == 1) {
		lamber = TRUE_;
		albedo = 0.;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 19, a__1[1] = ", Bottom Albedo = 0";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 2) {
		lamber = TRUE_;
		albedo = 1.;
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 19, a__1[1] = ", Bottom Albedo = 1";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
	    } else if (icas == 3) {
		lamber = FALSE_;
		i__1 = nstr;
		for (k = 0; k <= i__1; ++k) {
		    hl[(i__2 = k) < 49 && 0 <= i__2 ? i__2 : s_rnge("hl", 
			    i__2, "testdo_", (ftnlen)810)] = pow_di(&c_b209, &
			    k);
/* L69: */
		}
/* Writing concatenation */
		i__6[0] = lentit, a__1[0] = title;
		i__6[1] = 30, a__1[1] = ", Bottom Albedo = BDR Function";
		s_cat(header, a__1, i__6, &c__2, (ftnlen)127);
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)830)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)830)], &dochek_1.tstfup[(i__3 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tst\
fup", i__3, "testdo_", (ftnlen)830)], &dochek_1.tstdfd[(i__4 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__4 ? i__4 : s_rnge("tst\
dfd", i__4, "testdo_", (ftnlen)830)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    830)], &c__5, &c__10, &c__3);
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
	umu[0] = -1.;
	umu[1] = -.2;
	umu[2] = .2;
	umu[3] = 1.;
	nphi = 1;
	phi[0] = 60.;
	ibcnd = 0;
	fbeam = 0.;
	fisot = 1. / pi;
	lamber = TRUE_;
	albedo = 0.;
	plank = FALSE_;
	deltam = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    if (icas == 1) {
		dtauc[0] = .25;
		dtauc[1] = .25;
		ssalb[0] = .5;
		ssalb[1] = .3;
		ntau = 3;
		utau[0] = 0.;
		utau[1] = .25;
		utau[2] = .5;
		s_copy(header, "Test Case No. 8A:  Ref. OS, Table 1, Line 4 \
(Two Inhomogeneous Layers)", (ftnlen)127, (ftnlen)70);
	    } else if (icas == 2) {
		dtauc[0] = .25;
		dtauc[1] = .25;
		ssalb[0] = .8;
		ssalb[1] = .95;
		ntau = 3;
		utau[0] = 0.;
		utau[1] = .25;
		utau[2] = .5;
		s_copy(header, "Test Case No. 8b:  Ref. OS, Table 1, Line 1 \
(Two Inhomogeneous Layers)", (ftnlen)127, (ftnlen)70);
	    } else if (icas == 3) {
		dtauc[0] = 1.;
		dtauc[1] = 2.;
		ssalb[0] = .8;
		ssalb[1] = .95;
		ntau = 3;
		utau[0] = 0.;
		utau[1] = 1.;
		utau[2] = 3.;
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)927)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)927)], &dochek_1.tstfup[(i__3 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge("tst\
fup", i__3, "testdo_", (ftnlen)927)], &dochek_1.tstdfd[(i__4 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__4 ? i__4 : s_rnge("tst\
dfd", i__4, "testdo_", (ftnlen)927)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    927)], &c__5, &c__10, &c__3);
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
	i__1 = nlyr;
	for (lc = 1; lc <= i__1; ++lc) {
	    dtauc[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge("dtauc", 
		    i__2, "testdo_", (ftnlen)949)] = (doublereal) lc;
	    ssalb[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge("ssalb", 
		    i__2, "testdo_", (ftnlen)950)] = lc * .05 + .6;
/* L86: */
	}
	usrtau = TRUE_;
	ntau = 5;
	utau[0] = 0.;
	utau[1] = 1.05;
	utau[2] = 2.1;
	utau[3] = 6.;
	utau[4] = 21.;
	usrang = TRUE_;
	numu = 4;
	umu[0] = -1.;
	umu[1] = -.2;
	umu[2] = .2;
	umu[3] = 1.;
	nphi = 1;
	phi[0] = 60.;
	ibcnd = 0;
	fbeam = 0.;
	fisot = 1. / pi;
	lamber = TRUE_;
	deltam = TRUE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 3; ++icas) {
	    if (icas == 1) {
		i__1 = nlyr;
		for (lc = 1; lc <= i__1; ++lc) {
		    getmom_(&c__1, &c_b2, &nstr, &pmom[(i__2 = lc * 49 - 49) <
			     294 && 0 <= i__2 ? i__2 : s_rnge("pmom", i__2, 
			    "testdo_", (ftnlen)979)]);
/* L87: */
		}
		albedo = 0.;
		plank = FALSE_;
		s_copy(header, "Test Case No. 9a:  Ref. DGIS, Tables VI-VII,\
 beta=l=0 (multiple inhomogeneous layers)", (ftnlen)127, (ftnlen)85);
	    } else if (icas == 2) {
		pmom[0] = 1.;
		pmom[1] = .66971999999999998;
		pmom[2] = .31267800000000001;
		pmom[3] = .096295714285714276;
		pmom[4] = .024683333333333331;
		pmom[5] = .0042954545454545459;
		pmom[6] = 5.1615384615384609e-4;
		pmom[7] = 4.5333333333333335e-5;
		pmom[8] = 2.9411764705882355e-6;
		i__1 = nlyr;
		for (lc = 2; lc <= i__1; ++lc) {
		    for (k = 0; k <= 8; ++k) {
			pmom[(i__2 = k + lc * 49 - 49) < 294 && 0 <= i__2 ? 
				i__2 : s_rnge("pmom", i__2, "testdo_", (
				ftnlen)999)] = pmom[(i__3 = k) < 294 && 0 <= 
				i__3 ? i__3 : s_rnge("pmom", i__3, "testdo_", 
				(ftnlen)999)];
/* L88: */
		    }
		}
		s_copy(header, "Test Case No. 9b:  Ref. DGIS, Tables VI-VII,\
 beta=0,l=8 (multiple inhomogeneous layers)", (ftnlen)127, (ftnlen)87);
	    } else if (icas == 3) {
		temper[0] = 600.;
		i__1 = nlyr;
		for (lc = 1; lc <= i__1; ++lc) {
		    d__1 = (doublereal) lc / 7.;
		    getmom_(&c__3, &d__1, &nstr, &pmom[(i__2 = lc * 49 - 49) <
			     294 && 0 <= i__2 ? i__2 : s_rnge("pmom", i__2, 
			    "testdo_", (ftnlen)1008)]);
		    temper[(i__2 = lc) < 7 && 0 <= i__2 ? i__2 : s_rnge("tem\
per", i__2, "testdo_", (ftnlen)1009)] = lc * 10. + 600.;
/* L89: */
		}
		nphi = 3;
		phi[0] = 60.;
		phi[1] = 120.;
		phi[2] = 180.;
		fbeam = pi;
		umu0 = .5;
		phi0 = 0.;
		fisot = 1.;
		albedo = .5;
		plank = TRUE_;
		wvnmlo = 999.;
		wvnmhi = 1e3;
		btemp = 700.;
		ttemp = 550.;
		temis = 1.;
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
		    dochek_1.tstfir[(i__1 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__1 ? i__1 : s_rnge("tstfir", i__1, "testdo_"
		    , (ftnlen)1044)], &dochek_1.tstfdn[(i__2 = (icas + (nprob 
		    << 3)) * 5 - 45) < 360 && 0 <= i__2 ? i__2 : s_rnge("tst\
fdn", i__2, "testdo_", (ftnlen)1044)], &dochek_1.tstfup[(i__3 = (icas + (
		    nprob << 3)) * 5 - 45) < 360 && 0 <= i__3 ? i__3 : s_rnge(
		    "tstfup", i__3, "testdo_", (ftnlen)1044)], &
		    dochek_1.tstdfd[(i__4 = (icas + (nprob << 3)) * 5 - 45) < 
		    360 && 0 <= i__4 ? i__4 : s_rnge("tstdfd", i__4, "testdo_"
		    , (ftnlen)1044)], &dochek_1.tstuu[(i__5 = (((icas + (
		    nprob << 3)) * 3 + 1) * 10 + 1) * 5 - 1405) < 10800 && 0 
		    <= i__5 ? i__5 : s_rnge("tstuu", i__5, "testdo_", (ftnlen)
		    1044)], &c__5, &c__10, &c__3);
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
	temper[0] = 600.;
	i__1 = nlyr;
	for (lc = 1; lc <= i__1; ++lc) {
	    dtauc[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge("dtauc", 
		    i__2, "testdo_", (ftnlen)1066)] = (doublereal) lc;
	    ssalb[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge("ssalb", 
		    i__2, "testdo_", (ftnlen)1067)] = lc * .05 + .6;
	    d__1 = (doublereal) lc / (nlyr + 1);
	    getmom_(&c__3, &d__1, &nstr, &pmom[(i__2 = lc * 49 - 49) < 294 && 
		    0 <= i__2 ? i__2 : s_rnge("pmom", i__2, "testdo_", (
		    ftnlen)1068)]);
	    temper[(i__2 = lc) < 7 && 0 <= i__2 ? i__2 : s_rnge("temper", 
		    i__2, "testdo_", (ftnlen)1069)] = lc * 10. + 600.;
/* L97: */
	}
	usrtau = TRUE_;
	ntau = 3;
	utau[0] = 0.;
	utau[1] = 2.1;
	utau[2] = 21.;
	nphi = 2;
	phi[0] = 60.;
	phi[1] = 120.;
	ibcnd = 0;
	fbeam = pi;
	umu0 = .5;
	phi0 = 0.;
	fisot = 1.;
	lamber = TRUE_;
	deltam = TRUE_;
	albedo = .5;
	plank = TRUE_;
	wvnmlo = 999.;
	wvnmhi = 1e3;
	btemp = 700.;
	ttemp = 550.;
	temis = 1.;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    if (icas == 1) {
		usrang = TRUE_;
		numu = 4;
		umu[0] = -.788675129;
		umu[1] = -.211324871;
		umu[2] = .211324871;
		umu[3] = .788675129;
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
		i__1 = ntau;
		for (lu = 1; lu <= i__1; ++lu) {
		    cmpfir[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfir", i__2, "testdo_", (ftnlen)1131)] = 
			    rfldir[(i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : 
			    s_rnge("rfldir", i__3, "testdo_", (ftnlen)1131)];
		    cmpfdn[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfdn", i__2, "testdo_", (ftnlen)1132)] = rfldn[
			    (i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : s_rnge(
			    "rfldn", i__3, "testdo_", (ftnlen)1132)];
		    cmpfup[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfup", i__2, "testdo_", (ftnlen)1133)] = flup[(
			    i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : s_rnge(
			    "flup", i__3, "testdo_", (ftnlen)1133)];
		    cmpdfd[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpdfd", i__2, "testdo_", (ftnlen)1134)] = dfdt[(
			    i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : s_rnge(
			    "dfdt", i__3, "testdo_", (ftnlen)1134)];
		    i__2 = numu;
		    for (iu = 1; iu <= i__2; ++iu) {
			i__3 = nphi;
			for (j = 1; j <= i__3; ++j) {
			    cmpuu[(i__4 = lu + (iu + j * 10) * 5 - 56) < 150 
				    && 0 <= i__4 ? i__4 : s_rnge("cmpuu", 
				    i__4, "testdo_", (ftnlen)1137)] = uu[(
				    i__5 = iu + (lu + j * 5) * 10 - 61) < 150 
				    && 0 <= i__5 ? i__5 : s_rnge("uu", i__5, 
				    "testdo_", (ftnlen)1137)];
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
	umu[0] = -1.;
	umu[1] = -.1;
	umu[2] = .1;
	umu[3] = 1.;
	nphi = 2;
	phi[0] = 0.;
	phi[1] = 90.;
	ibcnd = 0;
	fbeam = 1.;
	umu0 = .5;
	phi0 = 0.;
	fisot = .5 / pi;
	lamber = TRUE_;
	albedo = .5;
	deltam = FALSE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    if (icas == 1) {
		nlyr = 1;
		dtauc[0] = 1.;
		ssalb[0] = .9;
		getmom_(&c__1, &c_b2, &nstr, pmom);
		usrtau = TRUE_;
		ntau = 4;
		utau[0] = 0.;
		utau[1] = .05;
		utau[2] = .5;
		utau[3] = 1.;
		prnt[1] = TRUE_;
		prnt[4] = TRUE_;
		s_copy(header, "Test Case No. 11a: One Isotropic-Scattering \
Layer", (ftnlen)127, (ftnlen)49);
	    } else if (icas == 2) {
		nlyr = ntau - 1;
		i__3 = nlyr;
		for (lc = 1; lc <= i__3; ++lc) {
		    dtauc[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge(
			    "dtauc", i__2, "testdo_", (ftnlen)1206)] = utau[(
			    i__1 = lc) < 5 && 0 <= i__1 ? i__1 : s_rnge("utau"
			    , i__1, "testdo_", (ftnlen)1206)] - utau[(i__4 = 
			    lc - 1) < 5 && 0 <= i__4 ? i__4 : s_rnge("utau", 
			    i__4, "testdo_", (ftnlen)1206)];
		    ssalb[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge(
			    "ssalb", i__2, "testdo_", (ftnlen)1207)] = .9;
		    getmom_(&c__1, &c_b2, &nstr, &pmom[(i__2 = lc * 49 - 49) <
			     294 && 0 <= i__2 ? i__2 : s_rnge("pmom", i__2, 
			    "testdo_", (ftnlen)1208)]);
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
		i__3 = ntau;
		for (lu = 1; lu <= i__3; ++lu) {
		    cmpfir[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfir", i__2, "testdo_", (ftnlen)1229)] = 
			    rfldir[(i__1 = lu - 1) < 5 && 0 <= i__1 ? i__1 : 
			    s_rnge("rfldir", i__1, "testdo_", (ftnlen)1229)];
		    cmpfdn[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfdn", i__2, "testdo_", (ftnlen)1230)] = rfldn[
			    (i__1 = lu - 1) < 5 && 0 <= i__1 ? i__1 : s_rnge(
			    "rfldn", i__1, "testdo_", (ftnlen)1230)];
		    cmpfup[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfup", i__2, "testdo_", (ftnlen)1231)] = flup[(
			    i__1 = lu - 1) < 5 && 0 <= i__1 ? i__1 : s_rnge(
			    "flup", i__1, "testdo_", (ftnlen)1231)];
		    cmpdfd[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpdfd", i__2, "testdo_", (ftnlen)1232)] = dfdt[(
			    i__1 = lu - 1) < 5 && 0 <= i__1 ? i__1 : s_rnge(
			    "dfdt", i__1, "testdo_", (ftnlen)1232)];
		    i__2 = numu;
		    for (iu = 1; iu <= i__2; ++iu) {
			i__1 = nphi;
			for (j = 1; j <= i__1; ++j) {
			    cmpuu[(i__4 = lu + (iu + j * 10) * 5 - 56) < 150 
				    && 0 <= i__4 ? i__4 : s_rnge("cmpuu", 
				    i__4, "testdo_", (ftnlen)1235)] = uu[(
				    i__5 = iu + (lu + j * 5) * 10 - 61) < 150 
				    && 0 <= i__5 ? i__5 : s_rnge("uu", i__5, 
				    "testdo_", (ftnlen)1235)];
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
	umu[0] = -1.;
	umu[1] = -.1;
	umu[2] = .1;
	umu[3] = 1.;
	nphi = 1;
	phi[0] = 0.;
	ibcnd = 0;
	fbeam = 1.;
	umu0 = 1.;
	phi0 = 0.;
	fisot = 0.;
	lamber = TRUE_;
	albedo = 1.;
	deltam = TRUE_;
	plank = FALSE_;
	onlyfl = FALSE_;
	for (icas = 1; icas <= 2; ++icas) {
	    if (icas == 1) {
		nlyr = 1;
		dtauc[0] = 20.1;
		ssalb[0] = .5;
		getmom_(&c__3, &c_b457, &nstr, pmom);
		usrtau = TRUE_;
		ntau = 4;
		utau[0] = 0.;
		utau[1] = 10.;
		utau[2] = 19.9;
		utau[3] = 20.1;
		prnt[1] = TRUE_;
		prnt[4] = TRUE_;
		s_copy(header, "Test Case No. 12a:  Overhead Beam Striking A\
bsorbing/Scattering Medium", (ftnlen)127, (ftnlen)70);
	    } else if (icas == 2) {
		nlyr = ntau - 1;
		i__1 = nlyr;
		for (lc = 1; lc <= i__1; ++lc) {
		    dtauc[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge(
			    "dtauc", i__2, "testdo_", (ftnlen)1304)] = utau[(
			    i__3 = lc) < 5 && 0 <= i__3 ? i__3 : s_rnge("utau"
			    , i__3, "testdo_", (ftnlen)1304)] - utau[(i__4 = 
			    lc - 1) < 5 && 0 <= i__4 ? i__4 : s_rnge("utau", 
			    i__4, "testdo_", (ftnlen)1304)];
		    ssalb[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge(
			    "ssalb", i__2, "testdo_", (ftnlen)1305)] = .5;
		    getmom_(&c__3, &c_b457, &nstr, &pmom[(i__2 = lc * 49 - 49)
			     < 294 && 0 <= i__2 ? i__2 : s_rnge("pmom", i__2, 
			    "testdo_", (ftnlen)1306)]);
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
		i__1 = ntau;
		for (lu = 1; lu <= i__1; ++lu) {
		    cmpfir[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfir", i__2, "testdo_", (ftnlen)1327)] = 
			    rfldir[(i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : 
			    s_rnge("rfldir", i__3, "testdo_", (ftnlen)1327)];
		    cmpfdn[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfdn", i__2, "testdo_", (ftnlen)1328)] = rfldn[
			    (i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : s_rnge(
			    "rfldn", i__3, "testdo_", (ftnlen)1328)];
		    cmpfup[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpfup", i__2, "testdo_", (ftnlen)1329)] = flup[(
			    i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : s_rnge(
			    "flup", i__3, "testdo_", (ftnlen)1329)];
		    cmpdfd[(i__2 = lu - 1) < 5 && 0 <= i__2 ? i__2 : s_rnge(
			    "cmpdfd", i__2, "testdo_", (ftnlen)1330)] = dfdt[(
			    i__3 = lu - 1) < 5 && 0 <= i__3 ? i__3 : s_rnge(
			    "dfdt", i__3, "testdo_", (ftnlen)1330)];
		    i__2 = numu;
		    for (iu = 1; iu <= i__2; ++iu) {
			i__3 = nphi;
			for (j = 1; j <= i__3; ++j) {
			    cmpuu[(i__4 = lu + (iu + j * 10) * 5 - 56) < 150 
				    && 0 <= i__4 ? i__4 : s_rnge("cmpuu", 
				    i__4, "testdo_", (ftnlen)1333)] = uu[(
				    i__5 = iu + (lu + j * 5) * 10 - 61) < 150 
				    && 0 <= i__5 ? i__5 : s_rnge("uu", i__5, 
				    "testdo_", (ftnlen)1333)];
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
	phi0 = 0.;
	albedo = .5;
	deltam = TRUE_;
	for (icas = 1; icas <= 4; ++icas) {
	    if (icas == 1) {
		ibcnd = 1;
		nlyr = 1;
		dtauc[0] = 1.;
		ssalb[0] = .99;
		getmom_(&c__3, &c_b243, &nstr, pmom);
		prnt[5] = TRUE_;
		prnt[1] = FALSE_;
		usrang = TRUE_;
		numu = 1;
		umu[0] = .5;
		s_copy(header, "Test Case No. 13a:  Albedo and Transmissivit\
y from Shortcut, Single Layer", (ftnlen)127, (ftnlen)73);
	    } else if (icas == 2) {
		ibcnd = 0;
		usrtau = TRUE_;
		ntau = 2;
		utau[0] = 0.;
		utau[1] = 1.;
		umu0 = .5;
		fbeam = 1. / umu0;
		fisot = 0.;
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
		i__3 = nlyr;
		for (lc = 1; lc <= i__3; ++lc) {
		    dtauc[(i__2 = lc - 1) < 6 && 0 <= i__2 ? i__2 : s_rnge(
			    "dtauc", i__2, "testdo_", (ftnlen)1407)] = 1. / 
			    nlyr;
		    getmom_(&c__3, &c_b243, &nstr, &pmom[(i__2 = lc * 49 - 49)
			     < 294 && 0 <= i__2 ? i__2 : s_rnge("pmom", i__2, 
			    "testdo_", (ftnlen)1408)]);
/* L125: */
		}
		ssalb[0] = .99;
		ssalb[1] = .5;
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

/* Subroutine */ int getmom_(integer *iphas, doublereal *gg, integer *nmom, 
	doublereal *pmom)
{
    /* Initialized data */

    static doublereal hazelm[82] = { 2.4126,3.23047,3.37296,3.2315,2.8935,
	    2.49594,2.11361,1.74812,1.44692,1.17714,.96643,.78237,.64114,
	    .51966,.42563,.34688,.28351,.23317,.18963,.15788,.12739,.10762,
	    .08597,.07381,.05828,.05089,.03971,.03524,.0272,.02451,.01874,
	    .01711,.01298,.01198,.00904,.00841,.00634,.00592,.00446,.00418,
	    .00316,.00296,.00225,.0021,.0016,.0015,.00115,.00107,8.2e-4,
	    7.7e-4,5.9e-4,5.5e-4,4.3e-4,4e-4,3.1e-4,2.9e-4,2.3e-4,2.1e-4,
	    1.7e-4,1.5e-4,1.2e-4,1.1e-4,9e-5,8e-5,6e-5,6e-5,5e-5,4e-5,4e-5,
	    3e-5,3e-5,2e-5,2e-5,2e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5 }
	    ;
    static doublereal cldmom[299] = { 2.544,3.883,4.568,5.235,5.887,6.457,
	    7.177,7.859,8.494,9.286,9.856,10.615,11.229,11.851,12.503,13.058,
	    13.626,14.209,14.66,15.231,15.641,16.126,16.539,16.934,17.325,
	    17.673,17.999,18.329,18.588,18.885,19.103,19.345,19.537,19.721,
	    19.884,20.024,20.145,20.251,20.33,20.401,20.444,20.477,20.489,
	    20.483,20.467,20.427,20.382,20.31,20.236,20.136,20.036,19.909,
	    19.785,19.632,19.486,19.311,19.145,18.949,18.764,18.551,18.348,
	    18.119,17.901,17.659,17.428,17.174,16.931,16.668,16.415,16.144,
	    15.883,15.606,15.338,15.058,14.784,14.501,14.225,13.941,13.662,
	    13.378,13.098,12.816,12.536,12.257,11.978,11.703,11.427,11.156,
	    10.884,10.618,10.35,10.09,9.827,9.574,9.318,9.072,8.822,8.584,
	    8.34,8.11,7.874,7.652,7.424,7.211,6.99,6.785,6.573,6.377,6.173,
	    5.986,5.79,5.612,5.424,5.255,5.075,4.915,4.744,4.592,4.429,4.285,
	    4.13,3.994,3.847,3.719,3.58,3.459,3.327,3.214,3.09,2.983,2.866,
	    2.766,2.656,2.562,2.459,2.372,2.274,2.193,2.102,2.025,1.94,1.869,
	    1.79,1.723,1.649,1.588,1.518,1.461,1.397,1.344,1.284,1.235,1.179,
	    1.134,1.082,1.04,.992,.954,.909,.873,.832,.799,.762,.731,.696,
	    .668,.636,.61,.581,.557,.53,.508,.483,.463,.44,.422,.401,.384,
	    .364,.349,.331,.317,.301,.288,.273,.262,.248,.238,.225,.215,.204,
	    .195,.185,.177,.167,.16,.151,.145,.137,.131,.124,.118,.112,.107,
	    .101,.097,.091,.087,.082,.079,.074,.071,.067,.064,.06,.057,.054,
	    .052,.049,.047,.044,.042,.039,.038,.035,.034,.032,.03,.029,.027,
	    .026,.024,.023,.022,.021,.02,.018,.018,.017,.016,.015,.014,.013,
	    .013,.012,.011,.011,.01,.009,.009,.008,.008,.008,.007,.007,.006,
	    .006,.006,.005,.005,.005,.005,.004,.004,.004,.004,.003,.003,.003,
	    .003,.003,.003,.002,.002,.002,.002,.002,.002,.002,.002,.002,.001,
	    .001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,.001,
	    .001,.001,.001,.001 };

    /* System generated locals */
    integer pmom_dim1, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    double pow_di(doublereal *, integer *);

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
    /* Parameter adjustments */
    pmom_dim1 = *nmom - 0 + 1;

    /* Function Body */
    if (*iphas < 1 || *iphas > 5) {
	errmsg_("GETMOM--bad input variable IPHAS", &c_true, (ftnlen)32);
    }
    if (*iphas == 3 && (*gg <= -1. || *gg >= 1.)) {
	errmsg_("GETMOM--bad input variable GG", &c_true, (ftnlen)29);
    }
    if (*nmom < 2) {
	errmsg_("GETMOM--bad input variable NMOM", &c_true, (ftnlen)31);
    }
    pmom[(i__1 = 0) < 1 * pmom_dim1 ? i__1 : s_rnge("pmom", i__1, "getmom_", (
	    ftnlen)1554)] = 1.;
    i__1 = *nmom;
    for (k = 1; k <= i__1; ++k) {
	pmom[(i__2 = k) < 1 * pmom_dim1 && 0 <= i__2 ? i__2 : s_rnge("pmom", 
		i__2, "getmom_", (ftnlen)1556)] = 0.;
/* L10: */
    }
    if (*iphas == 2) {
/*                                       ** Rayleigh phase function */
	pmom[(i__1 = 2) < 1 * pmom_dim1 ? i__1 : s_rnge("pmom", i__1, "getmo\
m_", (ftnlen)1562)] = .1;
    } else if (*iphas == 3) {
/*                                       ** Henyey-Greenstein phase fcn */
	i__1 = *nmom;
	for (k = 1; k <= i__1; ++k) {
	    pmom[(i__2 = k) < 1 * pmom_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "pmom", i__2, "getmom_", (ftnlen)1567)] = pow_di(gg, &k);
/* L20: */
	}
    } else if (*iphas == 4) {
/*                                        ** Haze-L phase function */
	i__1 = min(82,*nmom);
	for (k = 1; k <= i__1; ++k) {
	    pmom[(i__2 = k) < 1 * pmom_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "pmom", i__2, "getmom_", (ftnlen)1573)] = hazelm[(i__3 = 
		    k - 1) < 82 && 0 <= i__3 ? i__3 : s_rnge("hazelm", i__3, 
		    "getmom_", (ftnlen)1573)] / ((k << 1) + 1);
/* L30: */
	}
    } else if (*iphas == 5) {
/*                                        ** Cloud C.1 phase function */
	i__1 = min(298,*nmom);
	for (k = 1; k <= i__1; ++k) {
	    pmom[(i__2 = k) < 1 * pmom_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "pmom", i__2, "getmom_", (ftnlen)1579)] = cldmom[(i__3 = 
		    k - 1) < 299 && 0 <= i__3 ? i__3 : s_rnge("cldmom", i__3, 
		    "getmom_", (ftnlen)1579)] / ((k << 1) + 1);
/* L40: */
	}
    }
    return 0;
} /* getmom_ */

/* Subroutine */ int prtfin_(doublereal *utau, integer *ntau, doublereal *umu,
	 integer *numu, doublereal *phi, integer *nphi, integer *maxulv, 
	integer *maxumu, logical *onlyfl, logical *azmavg, doublereal *rfldir,
	 doublereal *rfldn, doublereal *flup, doublereal *dfdt, doublereal *
	u0u, doublereal *uu, doublereal *tstfir, doublereal *tstfdn, 
	doublereal *tstfup, doublereal *tstdfd, doublereal *tstuu, integer *
	maxtau, integer *maxmu, integer *maxaz)
{
    /* Format strings */
    static char fmt_300[] = "(//,1x,45(\002=\002),/,a,i4,a,/,1x,45(\002=\002\
))";

    /* System generated locals */
    integer tstuu_dim1, tstuu_dim2, tstuu_dim3, tstuu_offset, u0u_dim1, 
	    u0u_offset, uu_dim1, uu_dim2, uu_offset, i__1, i__2, i__3, i__4, 
	    i__5;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal umax, ratv[100];
    static integer j;
    extern doublereal ratio_(doublereal *, doublereal *);
    static integer iu, lu, numbad;
    static doublereal fnoise, flxmax;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static doublereal unoise, rat1, rat2, rat3, rat4;

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
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2 * 1);
    u0u_dim1 = *maxumu;
    u0u_offset = 1 + u0u_dim1 * 1;
    tstuu_dim1 = *maxtau;
    tstuu_dim2 = *maxmu;
    tstuu_dim3 = *maxaz;
    tstuu_offset = 1 + tstuu_dim1 * (1 + tstuu_dim2 * 1);

    /* Function Body */
    if (*ntau > *maxtau || *numu > *maxmu || *nphi > *maxaz) {
	errmsg_("PRTFIN--out of bounds in comparator arrays", &c_true, (
		ftnlen)42);
    }
    flxmax = 0.;
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
/* Computing MAX */
	d__1 = flxmax, d__2 = tstfir[lu - 1], d__1 = max(d__1,d__2), d__2 = 
		tstfdn[lu - 1], d__1 = max(d__1,d__2), d__2 = tstfup[lu - 1];
	flxmax = max(d__1,d__2);
/* L5: */
    }
    fnoise = flxmax * 1e-6;
    if (flxmax <= 0.) {
	errmsg_("PRTFIN--all fluxes zero or negative", &c_false, (ftnlen)35);
    }
    if (fnoise <= 0.) {
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
	do_fio(&c__1, (char *)&utau[lu - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rfldir[lu - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rfldn[lu - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&flup[lu - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dfdt[lu - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
	rat1 = ratio_(&rfldir[lu - 1], &tstfir[lu - 1]);
	rat2 = ratio_(&rfldn[lu - 1], &tstfdn[lu - 1]);
	rat3 = ratio_(&flup[lu - 1], &tstfup[lu - 1]);
	rat4 = ratio_(&dfdt[lu - 1], &tstdfd[lu - 1]);
	s_wsfe(&io___84);
	do_fio(&c__1, (char *)&rat1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rat2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rat3, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rat4, (ftnlen)sizeof(doublereal));
	e_wsfe();
	if ((rat1 < .99 || rat1 > 1.01) && (d__1 = rfldir[lu - 1], abs(d__1)) 
		> fnoise) {
	    ++numbad;
	}
	if ((rat2 < .99 || rat2 > 1.01) && (d__1 = rfldn[lu - 1], abs(d__1)) 
		> fnoise) {
	    ++numbad;
	}
	if ((rat3 < .99 || rat3 > 1.01) && (d__1 = flup[lu - 1], abs(d__1)) > 
		fnoise) {
	    ++numbad;
	}
	if ((rat4 < .99 || rat4 > 1.01) && (d__1 = dfdt[lu - 1], abs(d__1)) > 
		fnoise) {
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
	umax = 0.;
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
/* Computing MAX */
		d__1 = umax, d__2 = tstuu[(i__3 = lu + (iu + tstuu_dim2) * 
			tstuu_dim1 - tstuu_offset) < 1 * tstuu_dim1 * 
			tstuu_dim2 * tstuu_dim3 && 0 <= i__3 ? i__3 : s_rnge(
			"tstuu", i__3, "prtfin_", (ftnlen)1716)];
		umax = max(d__1,d__2);
/* L15: */
	    }
/* L16: */
	}
	unoise = umax * 1e-6;
	if (umax <= 0.) {
	    errmsg_("PRTFIN--all az-avg intensities zero or negative", &
		    c_false, (ftnlen)47);
	}
	if (unoise <= 0.) {
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
	    do_fio(&c__1, (char *)&umu[iu - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    s_wsfe(&io___89);
	    do_fio(&c__1, (char *)&utau[lu - 1], (ftnlen)sizeof(doublereal));
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		do_fio(&c__1, (char *)&u0u[iu + lu * u0u_dim1 - u0u_offset], (
			ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		ratv[(i__3 = iu - 1) < 100 && 0 <= i__3 ? i__3 : s_rnge("ratv"
			, i__3, "prtfin_", (ftnlen)1736)] = ratio_(&u0u[iu + 
			lu * u0u_dim1 - u0u_offset], &tstuu[(i__4 = lu + (iu 
			+ tstuu_dim2) * tstuu_dim1 - tstuu_offset) < 1 * 
			tstuu_dim1 * tstuu_dim2 * tstuu_dim3 && 0 <= i__4 ? 
			i__4 : s_rnge("tstuu", i__4, "prtfin_", (ftnlen)1736)]
			);
		if ((ratv[(i__3 = iu - 1) < 100 && 0 <= i__3 ? i__3 : s_rnge(
			"ratv", i__3, "prtfin_", (ftnlen)1737)] < .99 || ratv[
			(i__3 = iu - 1) < 100 && 0 <= i__3 ? i__3 : s_rnge(
			"ratv", i__3, "prtfin_", (ftnlen)1737)] > 1.01) && (
			d__1 = u0u[iu + lu * u0u_dim1 - u0u_offset], abs(d__1)
			) > unoise) {
		    ++numbad;
		}
/* L20: */
	    }
	    s_wsfe(&io___91);
	    i__2 = *numu;
	    for (iu = 1; iu <= i__2; ++iu) {
		do_fio(&c__1, (char *)&ratv[(i__3 = iu - 1) < 100 && 0 <= 
			i__3 ? i__3 : s_rnge("ratv", i__3, "prtfin_", (ftnlen)
			1741)], (ftnlen)sizeof(doublereal));
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
	umax = 0.;
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    i__3 = *numu;
	    for (iu = 1; iu <= i__3; ++iu) {
		i__2 = *nphi;
		for (j = 1; j <= i__2; ++j) {
/* Computing MAX */
		    d__1 = umax, d__2 = tstuu[(i__4 = lu + (iu + j * 
			    tstuu_dim2) * tstuu_dim1 - tstuu_offset) < 1 * 
			    tstuu_dim1 * tstuu_dim2 * tstuu_dim3 && 0 <= i__4 
			    ? i__4 : s_rnge("tstuu", i__4, "prtfin_", (ftnlen)
			    1757)];
		    umax = max(d__1,d__2);
/* L34: */
		}
/* L35: */
	    }
/* L36: */
	}
	unoise = umax * 1e-6;
	if (umax <= 0.) {
	    errmsg_("PRTFIN--all intensities zero or negative", &c_false, (
		    ftnlen)40);
	}
	if (unoise <= 0.) {
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
	    do_fio(&c__1, (char *)&phi[j - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    i__3 = *numu;
	    for (iu = 1; iu <= i__3; ++iu) {
		if (iu == 1) {
		    s_wsfe(&io___94);
		    do_fio(&c__1, (char *)&utau[lu - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&umu[iu - 1], (ftnlen)sizeof(
			    doublereal));
		    i__2 = *nphi;
		    for (j = 1; j <= i__2; ++j) {
			do_fio(&c__1, (char *)&uu[iu + (lu + j * uu_dim2) * 
				uu_dim1 - uu_offset], (ftnlen)sizeof(
				doublereal));
		    }
		    e_wsfe();
		}
		if (iu > 1) {
		    s_wsfe(&io___95);
		    do_fio(&c__1, (char *)&umu[iu - 1], (ftnlen)sizeof(
			    doublereal));
		    i__2 = *nphi;
		    for (j = 1; j <= i__2; ++j) {
			do_fio(&c__1, (char *)&uu[iu + (lu + j * uu_dim2) * 
				uu_dim1 - uu_offset], (ftnlen)sizeof(
				doublereal));
		    }
		    e_wsfe();
		}
		i__2 = *nphi;
		for (j = 1; j <= i__2; ++j) {
		    ratv[(i__4 = j - 1) < 100 && 0 <= i__4 ? i__4 : s_rnge(
			    "ratv", i__4, "prtfin_", (ftnlen)1784)] = ratio_(&
			    uu[iu + (lu + j * uu_dim2) * uu_dim1 - uu_offset],
			     &tstuu[(i__5 = lu + (iu + j * tstuu_dim2) * 
			    tstuu_dim1 - tstuu_offset) < 1 * tstuu_dim1 * 
			    tstuu_dim2 * tstuu_dim3 && 0 <= i__5 ? i__5 : 
			    s_rnge("tstuu", i__5, "prtfin_", (ftnlen)1784)]);
		    if ((ratv[(i__4 = j - 1) < 100 && 0 <= i__4 ? i__4 : 
			    s_rnge("ratv", i__4, "prtfin_", (ftnlen)1785)] < 
			    .99 || ratv[(i__4 = j - 1) < 100 && 0 <= i__4 ? 
			    i__4 : s_rnge("ratv", i__4, "prtfin_", (ftnlen)
			    1785)] > 1.01) && (d__1 = uu[iu + (lu + j * 
			    uu_dim2) * uu_dim1 - uu_offset], abs(d__1)) > 
			    unoise) {
			++numbad;
		    }
/* L40: */
		}
		s_wsfe(&io___96);
		i__2 = *nphi;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, (char *)&ratv[(i__4 = j - 1) < 100 && 0 <= 
			    i__4 ? i__4 : s_rnge("ratv", i__4, "prtfin_", (
			    ftnlen)1789)], (ftnlen)sizeof(doublereal));
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
