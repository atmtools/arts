ABS jobname=o2_spin_rot
(use a jpl catalogue, the my2 calculates oxygen lines twice, once the)
(one in the my2 file, once the from the comp_o2_abs routine)
(=================================================)

oh    h2o   hf    hcn   hnc   co    h2co  no    o2    ho2   h2o2  hcl
110   200   110   110   110   110   110   110   0     110   110   110  
n2o   no2   o3_16 o3_17 o3_18 ch3cl clo   hocl  ocs   hno3  so2   cof2
110   110   110   110   110   110   110   110   110   110   110   110  
oclo  hbr   bro   clno3 h2so4 clooclh2o_l  dry-cont H2O-model
110   110   110   110   110   110   110    110       1
abs_offs
110
h2o         ../fascod/midlatitude-summer.H2O.dat
o2          ../fascod/midlatitude-summer.O2.dat

  RH      STVMR       HB      rel. humid.    strat. VMR    breakpoint [km]
  0       0          0       

   HL       HU      DH         lower and upper alt. limit, step width [km]
  0.000   100.000   1.000   

MOD_ATM [ MONTH  DAY [LAT.]]  [ FILENAME [p1]]  [ T       hpr       p    ]
13      ../fascod/midlatitude-summer.tp.dat

   FRMIN     FRMAX        [ CATNAME ]                      SENSITIVITY [K]
   0.0       0.0        o2_spin_rot.cat                     0.0

FREQUENCIES FILENAME                     f_min    f_max     df   [GHz]
o2_spin_rot.fre

END

SCAN

TBR z(obs) z/-p (tan) th-zen e-hor e-ver dTbrmax mode   intpol back refr source weight
820.000    25.000      0.0   0.00 0.00  1.0E+003 limb   combi  yes  no   Planck
absorption file  -->
antenna file     -->pencil_beam
T-bright file    -->./o2_spin_rot.tbr
END

angle incr. (pencil) : 0.020000
tangent altitude     : 0.000000 90.000000 1.000000

END
