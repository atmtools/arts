ABS use_mytran jobname=master_b
()
()
(=================================================)

oh    h2o   hf    hcn   hnc   co    h2co  no    o2    ho2   h2o2  hcl
110   0     110   0     110   0     0     0     0     110     0     0  
n2o   no2   o3_16 o3_17 o3_18 ch3cl clo   hocl  ocs   hno3  so2   cof2
0     0     0     0     0     0     0     0     0     0     0     0  
oclo  hbr   bro   clno3 h2so4 clooclh2o_l  dry-cont H2O-model
110   110   110   110   110   110   110     1       1
abs_offs
110
h2o         ../fascod/midlatitude-summer.H2O.dat
hcn         ../fascod/midlatitude-summer.HCN.dat
co          ../fascod/midlatitude-summer.CO.dat
h2co        ../fascod/midlatitude-summer.H2CO.dat
no          ../fascod/midlatitude-summer.NO.dat
o2          ../fascod/midlatitude-summer.O2.dat
h2o2        ../fascod/midlatitude-summer.H2O2.dat
hcl         ../fascod/midlatitude-summer.HCL.dat
n2o         ../fascod/midlatitude-summer.N2O.dat
no2         ../fascod/midlatitude-summer.NO2.dat
o3_16       ../fascod/midlatitude-summer.O3.dat
o3_17       ../fascod/midlatitude-summer.O3.dat
o3_18       ../fascod/midlatitude-summer.O3.dat
ch3cl       ../fascod/midlatitude-summer.CH3CL.dat
clo         ../fascod/midlatitude-summer.CLO.dat
hocl        ../fascod/midlatitude-summer.HOCL.dat
ocs         ../fascod/midlatitude-summer.OCS.dat
hno3        ../fascod/midlatitude-summer.HNO3.dat
so2         ../fascod/midlatitude-summer.SO2.dat
cof2        ../fascod/midlatitude-summer.COF2.dat

  RH      STVMR       HB      rel. humid.    strat. VMR    breakpoint [km]
  0       0          0       

   HL       HU      DH         lower and upper alt. limit, step width [km]
  0.000   100.000   1.000   

MOD_ATM [ MONTH  DAY [LAT.]]  [ FILENAME [p1]]  [ T       hpr       p    ]
13      ../fascod/midlatitude-summer.tp.dat

   FRMIN     FRMAX        [ CATNAME ]                      SENSITIVITY [K]
   0.0       0.0          master_b.my2                     0.0

FREQUENCIES FILENAME                     f_min    f_max     df   [GHz]
master_b.fre

END

SCAN

TBR z(obs) z/-p (tan) th-zen e-hor e-ver dTbrmax mode   intpol back refr source weight
820.000    25.000      0.0   0.00 0.00  1.0E+003 limb   combi  yes  no   Planck
absorption file  -->
antenna file     -->pencil_beam
T-bright file    -->./master_b.tbr
END

angle incr. (pencil) : 0.020000
tangent altitude     : 0.000000 50.000000 1.000000

END
