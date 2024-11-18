import pyarts
import numpy as np

f = pyarts.arts.convert.kaycm2freq(np.linspace(1, 21000, 101))

atm = pyarts.arts.AtmPoint()
atm.pressure = 1e4
atm.temperature = 250
atm["H2O"] = 1e-2
atm["O2"] = 0.21
atm["N2"] = 0.79
atm["CO2"] = 400e-6

data = pyarts.arts.PredefinedModelData.fromcatalog("predef/", ["H2O-ForeignContCKDMT400", "H2O-SelfContCKDMT400"])

self_abs400 = pyarts.arts.predef.get_self_h2o_ckdmt400(f, atm, data)

self_ref400 = np.array([9.19895412e-08, 1.56724440e-04, 2.82856717e-05, 6.95598923e-06,
       3.26021741e-06, 1.30146968e-06, 1.26394130e-06, 1.39463521e-05,
       3.35576519e-05, 3.15502097e-06, 3.31108597e-07, 1.18282771e-07,
       7.97020190e-08, 7.04029075e-08, 6.86534100e-08, 2.64100534e-07,
       3.57178971e-07, 4.70763060e-06, 3.33191920e-05, 3.09874591e-06,
       3.40522091e-07, 1.04499034e-07, 5.66934848e-08, 3.97640748e-08,
       1.25834549e-07, 2.11353623e-06, 1.54490489e-06, 1.13040350e-07,
       1.38123499e-08, 4.68425892e-09, 1.77269284e-09, 2.89023528e-09,
       2.93522724e-08, 2.04101182e-07, 1.11615892e-06, 1.39595298e-06,
       3.75762635e-08, 9.39283668e-09, 3.13359038e-09, 1.89102931e-09,
       6.33411889e-09, 1.54190618e-08, 1.43451992e-07, 1.79495815e-08,
       1.83742214e-09, 3.03435257e-10, 1.25666488e-10, 1.79143920e-10,
       4.69540665e-10, 9.64452638e-09, 3.99187711e-08, 5.57703268e-08,
       9.77638588e-09, 2.68146778e-09, 2.45941656e-10, 8.07501003e-11,
       1.71587142e-10, 8.45303813e-10, 5.96415982e-09, 7.73688897e-10,
       4.56682514e-10, 2.27232918e-11, 8.54905969e-12, 1.24676074e-11,
       9.41541712e-11, 1.25828334e-09, 5.26022416e-09, 5.78700198e-10,
       9.11749205e-10, 7.98467201e-11, 1.16361446e-11, 5.73396010e-12,
       4.19352927e-11, 4.01748629e-10, 7.64266594e-11, 5.24648428e-11,
       1.22751530e-11, 1.83485032e-12, 2.63307495e-12, 2.50027287e-11,
       4.32678532e-10, 2.84325413e-10, 2.86877282e-11, 7.72808631e-11,
       2.91940239e-11, 2.84419669e-12, 3.47592137e-12, 4.45288667e-11,
       5.03700993e-11, 3.79234832e-12, 5.05179066e-12, 4.30220255e-12,
       1.27623466e-12, 1.77847173e-11, 2.35445259e-10, 6.98084008e-11,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00])

assert np.allclose(self_abs400, self_ref400)

foreign_abs400 = pyarts.arts.predef.get_foreign_h2o_ckdmt400(f, atm, data)

foreign_ref400 = np.array([7.91178291e-08, 6.80489234e-04, 1.95375691e-05, 1.45158347e-06,
       1.02882514e-07, 1.04559006e-08, 2.68733334e-07, 9.90478517e-05,
       1.92478048e-04, 6.54931113e-06, 1.67555770e-07, 7.84695663e-09,
       6.10529657e-10, 1.75279155e-08, 4.43039844e-08, 9.02314267e-07,
       5.56607315e-07, 1.85288042e-05, 1.49567131e-04, 3.67707520e-06,
       2.46794603e-07, 1.19519114e-08, 5.18898089e-09, 5.09342108e-09,
       9.83767065e-08, 1.26554748e-05, 9.25099562e-06, 1.24665026e-07,
       8.91821896e-09, 5.33090610e-10, 1.74704276e-10, 1.62460746e-09,
       1.56069467e-07, 1.00743829e-06, 8.32086844e-06, 1.09905920e-05,
       1.22257995e-07, 7.66766175e-09, 3.13000466e-10, 4.89363751e-09,
       2.27558118e-08, 4.42963910e-08, 7.56031119e-07, 3.90825216e-08,
       2.81932805e-09, 1.14925235e-10, 3.19843627e-11, 6.09308177e-10,
       8.24749116e-10, 3.68946553e-08, 2.39701702e-07, 3.80742419e-07,
       4.18647353e-08, 3.99639786e-08, 3.89150014e-10, 3.00630833e-11,
       5.67983757e-10, 3.31871500e-09, 3.08575917e-08, 1.63418074e-09,
       2.35560720e-09, 4.09097686e-11, 2.81590732e-12, 9.96138288e-12,
       1.76302443e-10, 6.35694883e-09, 2.64495452e-08, 1.14134679e-09,
       4.61546160e-09, 3.99178717e-10, 2.53580859e-11, 3.68156251e-12,
       1.43691049e-10, 1.96479163e-09, 2.00176976e-10, 2.78103007e-10,
       1.63006623e-10, 1.68714297e-12, 1.10670804e-12, 4.98782142e-11,
       2.68005970e-09, 1.34450762e-09, 4.88713078e-11, 4.08641405e-10,
       7.82169609e-11, 1.64257321e-12, 3.89621735e-12, 1.63291803e-10,
       2.16510087e-10, 2.08946096e-12, 1.65459032e-11, 9.32021707e-12,
       3.50607065e-13, 1.97331201e-11, 4.40758127e-10, 5.57756762e-11,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00])

assert np.allclose(foreign_abs400, foreign_ref400)
