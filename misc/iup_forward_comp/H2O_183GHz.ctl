ABS  use_mytran  jobname=H2O_183GHz
(VMR-Model Nrs.: 0: specify filename to read, 1..9: div., 10: CalTechMod.)
(No calculation, no loading        :   100 + Model Nr.                   )
(No calculation, but load ModelNr. :   200 + Model Nr.                   )

oh    h2o   hf    hcn   hnc   co    h2co  no    o2    ho2   h2o2  hcl
110   0     111   113   110   110   110   110   110   110   110   110   
n2o   no2   o3_16 o3_17 o3_18 ch3cl clo   hocl  hno3  so2   oclo  hbr
110   110   110   110   110   110   111   110   110   113   101   113   
bro   clno3 h2so4 clooclN2-cont H2O-model
111   110   110   111   101     1     
h2o           H2O.midlatitude-summer.arts.vmr

  RH      STVMR       HB      rel. humid.    strat. VMR    breakpoint [km]
 0.00   0.00E+000    0.000

   HL       HU      DH         lower and upper alt. limit, step width [km]
  0.000  95.000   1.000

MOD_ATM [ MONTH  DAY [LAT.]]  [ FILENAME [p1]]  [ T       hpr       p    ]
  13     midlatitude-summer.iup_for.atm

   FRMIN     FRMAX        [ CATNAME ]                      SENSITIVITY [K]
    0.000 41200.000        H2O_183GHz.my2                        1.00E-001

FREQUENCIES FILENAME                     f(low)   f(high)   delta-f  [GHz]
H2O_183GHz.fre                           0.000      0.000   0.000000

END

SCAN

TBR z(obs) z/-p (tan) th-zen e-hor e-ver dTbrmax mode   intpol back refr source weight
 300.000    20.000    106.957 0.00 0.00 1.0E+003 limb   linear  yes  no   Planck
absorption file  -->
antenna file     -->pencil_beam
T-bright file    -->./H2O_183GHz.tbr
END

angle incr. (pencil) : 0.03
tangent altitude     : 0.00 90.00 1.00

END
