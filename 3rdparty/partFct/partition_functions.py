# this file contains following routines:
#
#

import numpy as N

########## main functions
def do_TIPS_all(infile='TIPS_70-3000K.DAT',lowT=150.,highT=300.):
  '''
  Complete chain for deriving partition function polynomial coefficients for
   full set of TIPS species.
  
  Output is printed to screen in a format to allow direct pasting into ARTS'
   partition_function_data.cc.

  Parameters
  ----------
  infile: str
      name & location of file containing the partition function table
  lowT : float
      lower limit of temperature range to fit over
  highT : float
      upper limit of temperature range to fit over
  '''
  qtt,tt,header = readQT_TIPS(infile)
  
  fit_all(qtt,tt,header=header,lowT=lowT,highT=highT)
  


def do_JPL(tag='44013',lowT=9.,highT=300.,do_vib=False):
#def do_JPL(tag='44013',lowT=150.,highT=300.,do_vib=True):
  '''
  Complete chain for deriving partition function polynomial coefficients for
   one JPL species.

  Output is printed to screen in a format to allow direct pasting into ARTS'
   partition_function_data.cc.
  For TIPS we usually use lowT=150K and highT=300K, respectively. However, 
   for a sufficient fit with JPL 9-300K looks more reasonable and is used here
   (that has the disadvantage, though, that at low T Q(T) can fall below zero).
  Usually JPL partition functions need to be corrected for vibrational
   contribution. However, for some species this seems to already be included,
   namely: C3H8 (d044013.cat), CH3OH (d032003.cat).

  Parameters
  ----------
  tag: str
      (JPL) species tag to identify species to process.
  lowT : float
      lower limit of temperature range to fit over
  highT : float
      upper limit of temperature range to fit over
  do_vib: boolean
      flag whether do perform vibrational correction of JPL partition functions
  '''
  T = get_T_JPL()
  QT,nn = get_data_JPL_1spec(tag,do_vib=do_vib)
  tt = set_tt(lowT,highT)
  qtt = interpol_qtt(QT,T,tt,nn)
  
  fit_all(qtt.reshape(-1,1),tt,header=None,lowT=lowT,highT=highT)
  


########## helper functions

#### general stuff
def fit_all(qtt,tt,header=None,lowT=150.,highT=300.):
  '''
  Fits and outputs polynomial coefficients for a bunch of species.
  
  E.g., for getting polynomial coefficients for all TIPS species at once. Output
   format is such that data can be directly pasted into ARTS'
   partition_function_data.cc.
  
  Parameters
  ----------
  qtt: array of float (dim: (Nt,Nspec))
      tabulated partition functions
  tt : array of float (dim: Nt)
      temperatures corresponding to the qtt
  header: list of str (dim: Nspec)
      species tags (mol_iso) from TIPS. following HITRAN conventions.
  lowT: float
      lower limit of temperature range to fit over
  highT: float
      upper limit of temperature range to fit over

  '''
  for i in N.arange(qtt.shape[1]):
    x = do_poly3fit(qtt[(tt>=lowT) & (tt<=highT),i], tt[(tt>=lowT) & (tt<=highT)])
    if (header is not None) and (len(header)==qtt.shape[1]):
      print('%s\nQcoeff(%.6e,  %.6e,  %.6e,  %.6e ));' %(header[i],x[0],x[1],x[2],x[3]))
    else:
      print('%4i - Qcoeff(%.6e,  %.6e,  %.6e,  %.6e ));' %(i,x[0],x[1],x[2],x[3]))  
  


def do_poly3fit(y,k):
  '''
  Least-square fit of 3rd order polynomial coefficients.
  
  Used for deriving partition functions coefficients for ARTS from tabulated
   TIPS or JPL-based values.

  Parameters
  ----------
  y : array (dim: Ny)
      measurement vector, i.e., the data the polynomial is fitted to.
      for partition function fits this are the tabulated partition functions.
  k : array (dim: Ny)
      the base of the polynomial. used for deriving the kernel matrix with the
       derivations of the polynomial regarding its coefficients.
      for partition function fits this is the temperature grid corresponding to
       the tabulated partition functions.

  Returns
  -------
  x : array (dim: 4)
      the polynomial coefficients (x[0]+x[1]*a+x[2]*a^2+x[3]*a^3)
  '''
  K=N.ones((y.size,4))
  K[:,1]=k
  K[:,2]=k**2
  K[:,3]=k**3
  KtK=N.dot(K.T,K)
  KtKi=N.linalg.inv(KtK)
  KtKiKt=N.dot(KtKi,K.T)
  x=N.dot(KtKiKt,y)
  return x
  


#### TIPS stuff
def readQT_TIPS(infile):
  '''
  Reads tabulated partition functions from modified-TIPS output.
  
  Partition fucntion tabulated over temperature for a bunch of species.
  
  Parameters
  ----------
  infile: str
      name & location of file containing the partition function table

  Returns
  -------
  qtt: array of float (dim: (Nt,Nspec))
      tabulated partition functions
  tt : array of float (dim: Nt)
      temperatures corresponding to the qtt
  header: list of str (dim: Nspec)
      species tags (mol_iso) from TIPS. following HITRAN conventions.
  '''
  
  #first we read and disentangle the header which has the species info
  file = open(infile,'r')
  header = file.readline()
  file.close()
  #TIPS outputted header needs some adjustment...
  while '_ ' in header: header = header.replace('_ ','_')
  #now separate one-liner into array of species tags and remove temp tag
  header = header.split()[1:]

  qtt = N.loadtxt(infile,skiprows=1) #first row contains already derived header info
  tt = qtt[:,0]
  qtt = qtt[:,1:]
  
  return qtt,tt,header
  


#### JPL stuff
def get_T_JPL():
  '''
  Inititializes JPL's partition function temperature grid.
  
  Returns
  -------
  tt : array of float (dim: 7)
      temperature grid of JPL's tabulated partition functions
  '''
  tt=N.array([9.375,18.75,37.50,75.00,150.0,225.0,300.0])
  return tt


#def readQT_JPL_1spec(tag,catdir=''):
#  '''
#  '''


def get_data_JPL_1spec(tag,do_vib=True):
  '''
  Returns JPL's tabulated partition functions, optionally corrected for
   vibrational contribution, and temperature exponent of Q(T) dependency for the
   selected species.
  
  This provides the partition functions on JPL's 7-element temperature set
   (order: increasing T).
  Currently only C3H8 implemented. If further are needed, put them here (the
   QT from the d0XXXXX.cat file and degF from the c0XXXXX.cat, where XXXXX
   represents the JPL species tag), or use readQT_JPL_1spec (well, once it is
   implemented).

  Parameters
  ----------
  tag: str
      JPL species tag

  Returns
  -------
  qtt : array of float (dim: 7)
      JPL partition functions, tabulated on JPL fixed temperature set
  '''
  T = get_T_JPL()
  nT = T.size
  if tag=='44013': #C3H8
    degF = 3
    QT = N.array([3656,10276,28986,83718,290290,721927,1575246],dtype=N.float)
    #that's the catdir version, which are the log10 of the Q(T), and sorted for decreasing T.
    #QT=N.array([3.0190, 2.8322, 2.5696, 2.1242, 1.6863, 1.2685, 0.9367])
    #QT=10**N.fliplr(QT.reshape(1,-1)).reshape(-1)
  elif tag=='32003': #CH3OH (this is in TIPS, but with Q(T)=const=0 which is bullshit)
    degF = 3
    QT = N.array([19.5433,68.7464,230.2391,731.0698,2437.7654,5267.8635,9473.1198],dtype=N.float)
  elif tag=='16001': #O (this is in TIPS, but with Q(T)=const=0 which is bullshit)
    degF = 0
    QT = N.array([5.000,5.000,5.007,5.156,5.770,6.324,6.741],dtype=N.float)
  #elif tag=='':
  else:
    assert 1==0, \
      'No partition function data available for tag %s (tag is unknown).' %tag
  assert QT.size==nT, \
    'Something went wrong. QT must have %i elemnts, but has %i' %(nT,QT.size)
  if do_vib:
    vib = get_vib_levels(tag)
    QT = add_vib(QT,T,vib)
  nn = get_Texponent(degF)
  return QT,nn
  


def get_vib_levels(tag):
  '''
  Return the vibrational energy levels of a species.
  
  Those are required for correcting partition functions, e.g., from JPL for the
   vibrational part.
  Currently only implemented for C3H8. If other species are needed, add those
   (for the vib levels, have a look in partition_function.pro. there are plenty
   of them implemented already.)

  Parameters
  ----------
  tag: str
      species tag. following JPL convention, as this is basically needed for
       JPL data.

  Returns
  -------
  vib : array of float
      vibrational energy levels of species TAG. If not indicated otherwise,
      vib levels are taken from webbook.nist.gov/chemistry/form-ser.html
  '''
  if tag=='44013': #C3H8
    #table had some levels listed more than once with different mode type,
    # hence we assume they are two different transistions and have to be
    # handled additively.
    vib=N.array([  216,  268,  369,  748,  869,  922,  940, 1054, 1158, 1192,
                  1278, 1338, 1378, 1392, 1451, 1462, 1464, 1472, 1476, 2887,
                  2887, 2962, 2967, 2968, 2968, 2973, 2977],dtype=N.float)
  elif tag=='32003': #CH3OH
    vib=N.array([  200,  295, 1033, 1060, 1165, 1345, 1455, 1477, 1477, 2844,
                  2960, 3000, 3681])
  elif tag=='16001': #O; has no vib levels
    vib=N.array([])
  #elif tag=='':
  else:
    print('No vib data available for tag %s.' %tag)
    vib=N.array([],dtype=N.float)
  return vib
  


def get_Texponent(degF):
  '''
  Sets the exponent for temperature dependence of partition functions from the
   degree of freedom of a species.

  Parameters
  ----------
  degF: integer
      degree of freedom of the species (lin:2, non-lin:3)

  Returns
  -------
  nn : float
      temperature exponent of Q(T)-relation   
  '''
  if degF==2:
    nn = 1.
  elif degF==3:
    nn = 1.5
  elif degF==0:
    #that's for atomic species. no temp.exp. found in literature. linear in T
    # seems safest.
    nn = 1.
  else:
    assert 1==0, \
      'Sorry, temperature exponent for degree of freedom = %i is unknown.' %degF
  return nn
  


def add_vib(qtt,tt,vib):
  '''
  Add vibrational corrections to partition functions.
  
  This is, e.g., necessary for [most of] JPL delivered Q(T).
  Code taken from partition_function.pro and converted to python.

  Parameters
  ----------
  qtt: array of float(dim: Nt)
      uncorrected partition functions
  tt : array of float(dim: Nt)
      temperatures corresponding to the qtt
  vib : array of float (dim: Nvib)
      vibrational energy levels (in cm^-1)

  Returns
  -------
  qttv : array of float (dim: Nt)
      vibrationally corrected partition functions   
  '''
  
  qttv=N.copy(qtt)
  if (vib.size>0):
    cm2t= 1.43875 #conversion factor cm^-1 -> K
    for i in N.arange(qtt.size):
      for j in N.arange(vib.size):
        a=-vib[j]*cm2t/tt[i]
        if a>=-80: qttv[i]=qttv[i]/(1-N.exp(a))
  return qttv
  


def set_tt(lowT,highT,step=1.):
  '''
  Sets a regular grid between lowT and highT with stepwidth step
  
  Used for setting the (fine) temperature grid for partition functions to
   calculate on from JPL data

  Parameters
  ----------
  lowT : float
      lower limit of temperature grid
  highT : float
      upper limit of temperature grid
  step : float
      step width of regular grid

  Returns
  -------
  tt : array of float (dim: Nt)
      regular spaced (temperature) grid   
  '''
  tt = N.arange(lowT,highT+step,step)
  return tt
  


def interpol_qtt(QT,T,tt,nn):
  '''
  Interpolates partition function over temperature
  
  Following JPL rules for temperature dependence, i.e., partition functions
   proportional to T^1.0 for linear molecules and T^1.5 for non-linear molecules.


  Parameters
  ----------
  QT: array of float(dim: NT)
      partition functions on given temperature grid
  T : array of float(dim: NT)
      temperature grid of given partition functions
  tt : array of float(dim: Nt)
      temperatures to interpolate partition functions to
  nn : float
      temperature exponent

  Returns
  -------
  qtt : array of float (dim: Nt)
      interpolated partition functions   

  '''
  qtt = N.zeros_like(tt)
  Tn = T**nn
  ttn = tt**nn
  qtt[tt<=T[1]] = QT[0] + (QT[1]-QT[0])/(Tn[1]-Tn[0]) * (ttn[tt<=T[1]]-Tn[0])
  for i in N.arange(2,T.size):
    qtt[(tt>T[i-1]) & (tt<=T[i])] = QT[i-1] + \
      (QT[i]-QT[i-1])/(Tn[i]-Tn[i-1]) * (ttn[(tt>T[i-1]) & (tt<=T[i])]-Tn[i-1])
  qtt[tt>T[-1]] = QT[-2] + (QT[-1]-QT[-2])/(Tn[-1]-Tn[-2]) * (ttn[tt>T[-1]]-Tn[-2])
  return qtt
