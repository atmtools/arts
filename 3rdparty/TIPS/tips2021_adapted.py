from pyarts.classes.PartitionFunctionsData import PartitionFunctionsData
import numpy as np
import pickle
import sys

# MIT LICENSE
#
# COPYRIGHT (C) [2021] [ROBERT GAMACHE]
#
# PERMISSION IS HEREBY GRANTED, FREE OF CHARGE, TO ANY PERSON OBTAINING A COPY
# OF THIS SOFTWARE AND ASSOCIATED DOCUMENTATION FILES (THE "SOFTWARE"), TO DEAL
# IN THE SOFTWARE WITHOUT RESTRICTION, INCLUDING WITHOUT LIMITATION THE RIGHTS
# TO USE, COPY, MODIFY, MERGE, PUBLISH, DISTRIBUTE, SUBLICENSE, AND/OR SELL
# COPIES OF THE SOFTWARE, AND TO PERMIT PERSONS TO WHOM THE SOFTWARE IS
# FURNISHED TO DO SO, SUBJECT TO THE FOLLOWING CONDITIONS:
#
# THE ABOVE COPYRIGHT NOTICE AND THIS PERMISSION NOTICE SHALL BE INCLUDED IN ALL
# COPIES OR SUBSTANTIAL PORTIONS OF THE SOFTWARE.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#  --  UPDATES  --
#    C    FOR UPDATES SEE SEE PUBLICATION 
#    GAMACHE ET AL., Total Internal Partition Sums for the HITRAN2020 database, JQSRT, ??, ??, 2021. 
#
#    THIS PROGRAM CALCULATES THE TOTAL INTERNAL
#    PARTITION SUM (TIPS) FOR A GIVEN MOLECULE,ISOTOPOLOGUE, AND
#    TEMPERATURE.  CURRENT LIMITATIONS ARE THE MOLECULAR SPECIES ON THE
#    HITRAN MOLECULAR DATABASE PLUS A FEW ADDITIONAL MOLECULES AND THE TEMPERATURE RANGE IS GENERALLY 1 - 5000 K.

# ADAPTED FOR ARTS
#
# This file has been adapted from the original to output ARTS style data
# The data must still be imported from the referenced article.
#
# To use:
# python THIS_FILE INPUT_DATA_DIRECTORY OUTPUT_DATA_DIRECTORY
#
# The INPUT_DATA_DIRECTORY must be a directory that contains all of the *.QTpy
# files associated with the TIPS2021 data.
#
# The OUTPUT_DATA_DIRECTORY must exist and will have many XML-files written to it
#
# (Richard Larsson, ric.larsson@gmail.com, 2021-06-08)

mol_id = ['1 = H2O','2 = CO2','3 = O3','4 = N2O','5 = CO','6 = CH4','7 = O2',
'8 = NO','9 = SO2','10 = NO2','11 = NH3','12 = HNO3','13 = OH','14 = HF','15 = HCl',
'16 = HBr','17 = HI','18 = ClO','19 = OCS','20 = H2CO','21 = HOCl','22 = N2','23 = HCN',
'24 = CH3Cl','25 = H2O2','26 = C2H2','27 = C2H6','28 = PH3','29 = COF2','30 = SF6','31 = H2S',
'32 = HCOOH','33 = HO2','no 34     ','35 = ClONO2','36 = NO+','37 = HOBr','38 = C2H4','39 = CH3OH',
'40 = CH3Br','41 = CH3CN','42 = CF4','43 = C4H2','44 = HC3N','45 = H2','46 = CS','47 = SO3','48=C2N2',
'49=COCl2','50=SO','51=CH3F','52=GeH4','53=CS2','54=CH3I','55=NF3','56=C3H4','57=CH3',]

molecules = [' ','H2O   ','CO2   ','O3    ','N2O   ','CO    ','CH4   ','O2    ',
'NO    ','SO2   ','NO2   ','NH3   ','HNO3  ','OH    ','HF    ','HCl   ',
'HBr   ','HI    ','ClO   ','OCS   ','H2CO  ','HOCl  ','N2    ','HCN   ',
'CH3Cl ','H2O2  ','C2H2  ','C2H6  ','PH3   ','COF2  ','SF6   ','H2S   ',
'HCOOH ','HO2   ','O     ','ClONO2','NO+   ','HOBr  ','C2H4  ','CH3OH ',
'CH3Br ','CH3CN ','CF4   ','C4H2  ','HC3N  ','H2    ','CS    ','SO3   ','C2N2  ', 
'COCl2 ','SO    ','CH3F  ','GeH4  ','CS2   ','CH3I   ','NF3    ','C3H4  ','CH3   ']

niso = [0,9,13,18,5,9,4,6,3,4,2,2,2,3,2,4,4,2,2,6,3,2,3,3,2,1,3,3,1,
2,1,3,1,1,1,2,1,2,3,1,2,4,1,1,6,2,4,1,2,2,3,1,5,4,2,1,1,1]


# H2O
Tmax = [0, 5000.,5000.,5000.,5000.,5000.,5000.,6000.,6000.,6000.,
# CO2
5000.,5000.,3500.,3500.,3500.,3500.,5000.,3500.,5000.,5000.,3500.,5000.,5000.,
# O3
1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,
# N2O,                           CO
5000.,5000.,5000.,5000.,5000.,   9000.,9000.,9000.,9000.,9000.,9000.,9000.,9000.,9000.,
# CH4,                     O2,                                    NO,                  SO2
2500.,2500.,2500.,2500.,   7500.,7500.,7500.,7500.,7500.,7500.,   5000.,5000.,5000.,   5000.,5000.,5000.,5000.,
# NO2,         NH3,           HNO3,          OH,                  HF,            HCl
1000.,1000.,   6000.,6000.,   3500.,3500.,   9000.,5000.,5000.,   6000.,6000.,   6000.,6000.,6000.,6000.,
# HBr,                     HI,            ClO,           OCS,                                   H2CO
6000.,6000.,6000.,6000.,   6000.,6000.,   5000.,5000.,   5000.,5000.,5000.,5000.,5000.,5000.,   3500.,5000.,5000.,
# HOCl,        N2,                  HCN,                 CH3Cl,         H2O2,    C2H2, 
5000.,5000.,   9000.,9000.,9000.,   3500.,3500.,3500.,   5000.,5000.,   6000.,   5000.,5000.,5000.,   
# C2H6,              PH3,     COF2,          SF6,     H2S,                 HCOOH,   HO2,   O atom, ClONO2,
5000.,5000.,5000.,   4500.,   3500.,3500.,   5000.,   4000.,5000.,5000.,   5000.,   5000.,   0.,   5000.,5000.,   
#  NO+,  HOBr,          C2H4,                CH3OH,   CH3Br,         CH3CN,                     CF4
5000.,   5000.,5000.,   5000.,5000.,5000.,   3500.,   5000.,5000.,   5000.,5000.,5000.,5000.,   3010.,
# C4H2,  HC3N,                                  H2,            CS,                        SO3, 
5000.,   5000.,5000.,5000.,5000.,5000.,5000.,   6000.,6000.,   5000.,5000.,5000.,5000.,   3500.,
# C2N2,        COCl2,        SO,                 CH3F    GeH4                          
5000.,5000.,   5000.,5000.,  5000.,5000.,5000.,  5000.,  5000.,5000.,5000.,5000.,5000.,
# CS2                       CH3I,         NF3     C3H4,    CH3,  
  5000.,5000.,5000.,5000.,  5000.,5000.,  5000.,  5000.,   5000.]


isotopologue = list()
isotopologue = ['     ',\
#1  H2O
'1=161, 2=181, 3=171, 4=162, 5=182, 6=172, 7=262,\
 8=282, 9=272',
#2  CO2
' 1=626, 2=636, 3=628, 4=627,  5=638,  6=637,  7=828,  8=827,\
 9=727, 10=838, 11=837, 12=737, 13=646',
#3  O3
' 1=666, 2=668, 3=686, 4=667, 5=676, 6=886, 7=868, 8=678, 9=768,\
 10 = 786, 11=776, 12=767, 13=888, 14=887, 15=878, 16=778,\
  17=787, 18=777',
#4  N2O
' 1=446, 2=456, 3=546, 4=448, 5=447',
#5   CO
' 1=26, 2=36, 3=28, 4=27, 5=38, 6=37, 7=46, 8=48, 9=47',
#6   CH4
' 1=211, 2=311, 3=212, 4=312',
#7   O2
' 1=66, 2=68, 3=67, 4=88, 5=87, 6=77',
#8   NO
' 1=46, 2=56, 3=48',
#9   SO2
' 1=626, 2=646, 3=636, 4=628',
#10   NO2
' 1=646, 2=656',
#11   NH3
' 1=4111, 2=5111',
#12   HNO3
' 1=146, 2=156',
#13   OH
' 1=61, 2=81, 3=62',
#14   HF
' 1=19, 2=29',
#15   HCl
' 1=15, 2=17, 3=25, 4=27',
#16   HBr
' 1=19, 2=11, 3=29, 4=21',
#17   HI
' 1=17, 2=27',
#18   ClO
' 1=56, 2=76',
#19   OCS
' 1=622, 2=624, 3=632, 4=623, 5=822, 6=634',
#20   H2CO
' 1=126, 2=136, 3=128',
#21   HOCl
' 1=165, 2=167',
#22   N2
' 1=44, 2=45, 3=55',
#23   HCN
' 1=124, 2=134, 3=125',
#24   CH3Cl
' 1=215, 2=217',
#25   H2O2
' 1=1661',
#26   C2H2
' 1=1221, 2=1231, 3=1222',
#27   C2H6
' 1=1221, 2=1231, 3=1222',
#28   PH3
 '1=1111',
#29   COF2
' 1=269, 2=369',
#30   SF6
' 1=29',
#31   H2S
' 1=121, 2=141, 3=131',
#32   HCOOH
' 1=126',
#33   HO2
' 1=166',
#34   O
' No values',
#35   ClONO2
'1=5646, 2=7646',
#36   NO+
' 1=46',
#37   HOBr
' 1=169, 2=161',
#38   C2H4
' 1=112211, 2=112311, 3=112212',
#39   CH3OH
' 1=2161',
#40   CH3Br
' 1=219, 2=211',
#41   CH3CN
' 1=2124, 2=3124, 3=2134, 4=3134',
#42   CF4
' 1=29',
#43   C4H2
' 1=2211',
#44   HC3N
' 1=12224, 2=12225, 3=12234, 4=12324, 5=13224, 6=22224',
#45   H2
' 1=11, 2=12',
#46   CS
' 1=22, 2=24, 3=32, 4=23',
#47   SO3
' 1=26',
#48   C2N2
' 1=4224, 2=5225',
#49   COCl2
'1=2655, 2=2657',
#50   SO
'1=26, 2=46, 3=28',
#51   CH3F
'1=219',
#52 GeH4
'1=7411, 2=7211, 3=7011, 4=7311, 5=7611',
#53   CS2
'1=222, 2=244, 3=223, 4=232',
#54   CH3I
'1=271, 2=317',
#55   NF3
'1=49',
#56   C3H4
'1=1221',
#57   CH3
'1=2111']

assert len(sys.argv) == 3, f"Call should be: python {sys.argv[0]} INDIR OUTDIR"
idir = sys.argv[1]
odir = sys.argv[2]

def extract_isonames(mol):
    out = []
    for isot in isotopologue[mol].split(","):
        isot = isot.strip()
        l = isot.split('=')
        if len(l) > 1:
            out.append([int(l[0]), l[1]])
            for i in range(2, len(l)):
                out.append([int(l[0]), l[i]])
    return out


for mol in range(1, len(molecules)):
    molname = molecules[mol].strip()
    isonames = extract_isonames(mol)
    for iso, isoname in isonames:
        fullname = f"{molname.strip()}-{isoname.strip()}"
        file = f'{idir}/{mol}_{iso}.QTpy'

        QTdict = {}
        with open(file, 'rb') as handle:
            QTdict = pickle.loads(handle.read())
    
        global_ID = 0
        for I in range(1,int(mol)):
            global_ID = global_ID + niso[I]
        global_ID = global_ID + int(iso)
        
        Ts = np.linspace(1., Tmax[global_ID], int(Tmax[global_ID]))
        Qs = np.array([float(QTdict[str(int(T))]) for T in Ts])
        data = np.append(Ts, Qs).reshape(2, len(Ts)).T
        PartitionFunctionsData("Interp", data).savexml(f"{odir}/{fullname}.xml")
