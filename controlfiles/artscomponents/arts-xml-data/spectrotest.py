import numpy as N
import pylab as p

def get_data(fa3h=None,fa4h=None,fa4p=None):
  '''
  Reads data (cross sections, frequency grid, pressure grid, absorption species)
   from ARTS lookup table files.

  Parameters
  ----------
  fa3h: string
      name of ARTS3 catalogue file from HITRAN data
  fa4h: string
      name of ARTS4 catalogue file from HITRAN data
  fa4p:
      name of ARTS4 catalogue file from toolbox spectroscopic line data
  NOTE 1: not all these files have to be given. 'None' will create ampty
   cross-section output. If none of them is given, also the other output will be
   empty.
  NOTE 2: all of them might be used for other ARTS LUT data, too. Basically we
   can process 3 LUT files, whatever they contain (but all of them should be for
   identical frequency, pressure, and species grids!!!).

  Returns
  -------
  freq: array of float (dim: nfreq)
      Frequency grid
  pres: array of float (dim: npres)
      Pressure grid
  absspec: list of string (dim: nspec)
      Species name tags
  axs_a3h: array of float (dim: (ntemp,nspec,nfreq,npres))
      absorption cross section tensor from ARTS3 HITRAN data
  axs_a4h: array of float (dim: (ntemp,nspec,nfreq,npres))
      absorption cross section tensor from ARTS4 HITRAN data
  axs_a4p: array of float (dim: (ntemp,nspec,nfreq,npres))
      absorption cross section tensor from ARTS4 Toolbox data
  NOTE: ntemp is assumed to be 1.
  '''
  if fa3h is not None:
    f=open(fa3h)
  elif fa4h is not None:
    f=open(fa4h)
  elif fa4p is not None:
    f=open(fa4p)
  else:
    absspec = []
    freq = N.array([])
    pres = N.copy(freq)
    axs_a3h = N.copy(freq)
    axs_a4h = N.copy(freq)
    axs_a4p = N.copy(freq)
    return absspec, freq, pres, axs_a3h, axs_a4h, axs_a4p

  lines=f.readlines()
  f.close()

  absspec=[x for x in lines if '<SpeciesTag>' in x]
  for i,l in enumerate(absspec): absspec[i]=l.partition('"')[2].partition('-')[0]

  ind = first_substring(lines,"FrequencyGrid")
  nfreq = N.int(lines[ind].partition('nelem="')[-1].partition('"')[0])
  freq = N.zeros(nfreq)
  for i in N.arange(nfreq):
    freq[i] = N.float(lines[ind+i+1])
  
  ind = first_substring(lines,"PressureGrid")
  npres = N.int(lines[ind].partition('nelem="')[-1].partition('"')[0])
  pres = N.zeros(npres)
  for i in N.arange(npres):
    pres[i]=N.float(lines[ind+i+1])

  ind = first_substring(lines,"AbsorptionCrossSections")
  if fa3h is None:
    axs_a3h = N.array([])
  else:
    axs_a3h=N.loadtxt(fa3h,skiprows=ind+1,comments='<').reshape((-1,len(absspec),nfreq,npres))
  if fa4h is None:
    axs_a4h = N.array([])
  else:
    axs_a4h=N.loadtxt(fa4h,skiprows=ind+1,comments='<').reshape((-1,len(absspec),nfreq,npres))
  if fa4p is None:
    axs_a4p = N.array([])
  else:
    axs_a4p=N.loadtxt(fa4p,skiprows=ind+1,comments='<').reshape((-1,len(absspec),nfreq,npres))

  return absspec, freq, pres, axs_a3h, axs_a4h, axs_a4p


def fast_stats(absspec,axs1,axs2,axs3=None):
  '''
  Calculates and prints out some deviation statistics (max and mean relative
   deviation) between two (or three) sets of absorption cross sections.

  Parameters
  ----------
  absspec: string
      Species name tags corresponding to species in axs1/2 (2nd dimension)
  axs1/axs2/axs3: array of float (dim: (ntemp,nspec,nfreq,npres))
      The absorption cross section tensors. Sets axs1 and axs2 are mandatory. If
       axs3 is present, it is compared to the axs1.    

  Returns
  -------
  no returns
  '''
  if axs3 is not None:
    for i in N.arange(axs1.shape[1]):
      print('%10.1e %10.1e %10.1e %10.1e'
            %(#absspec[i],
              N.nanmax(2. * abs(axs1[0,i,:,:]-axs2[0,i,:,:]) /
                            (axs1[0,i,:,:]+axs2[0,i,:,:]) ),
              N.nansum(2. * abs(axs1[0,i,:,:]-axs2[0,i,:,:])/
                            (axs1[0,i,:,:]+axs2[0,i,:,:])) /
                N.isfinite(1./(axs1[0,i,:,:]+axs2[0,i,:,:])).sum(),
              N.nanmax(2. * abs(axs1[0,i,:,:]-axs3[0,i,:,:]) /
                            (axs1[0,i,:,:]+axs3[0,i,:,:]) ),
              N.nansum(2. * abs(axs1[0,i,:,:]-axs3[0,i,:,:])/
                            (axs1[0,i,:,:]+axs3[0,i,:,:])) /
                N.isfinite(1./(axs1[0,i,:,:]+axs3[0,i,:,:])).sum()  ) )
  else:
    for i in N.arange(axs1.shape[1]):
      print('%10s %10.1e %10.1e'
            %(absspec[i],
              N.nanmax(2. * abs(axs1[0,i,:,:]-axs2[0,i,:,:]) /
                            (axs1[0,i,:,:]+axs2[0,i,:,:]) ),
              N.nansum(2. * abs(axs1[0,i,:,:]-axs2[0,i,:,:])/
                            (axs1[0,i,:,:]+axs2[0,i,:,:])) /
                N.isfinite(1./(axs1[0,i,:,:]+axs2[0,i,:,:])).sum()  ) )



def plotdiff(axs1,axs2,freq,absspec,
             start=0,ende=None,
             newfig=False,save=False,
             title='HITRAN vs. Toolbox for ',
             casename='H-vs-TB'):
  '''
  Plots deviations of two sets of absorption cross sections.
  
  Plotting is done over frequency for all pressure levels present into one plot.
   Separate plots for separate species.
  Simultaneously does four different type of plots:
  (A) relative deviations
    (A1) log plot of abs(axs1-axs2) / 0.5*(axs1+axs2)
    (A2) linear plot of (axs1-axs2) / 0.5*(axs1+axs2)
  (B) absolute deviations
    (B1) linear plot of (axs1-axs2)
    (B2) log plot of abs(axs1-axs2)

  Parameters
  ----------
  axs1/axs2: array of float (dim: (ntemp,nspec,nfreq,npres))
      The absorption cross section tensors.
  freq: array of float (dim: nfreq)
      Frequency grid.
  absspec: list of string (dim: nspec)
      Species name tags corresponding to axs1 and axs2 2nd dimension.
  start: integer
      Index of absspec species to start with (earlier entries are ignored).
  ende: integer
      Index of last absspec species to be plotted (later entries are ignored).
  newfig: boolean
      Flag whether to open a new figure for each plot (=True; recommended when
       plotting on screen). If False, figure content will be erased before each
       new plot (recommended when saving the figure, avoiding to open plenty of
       on-screen figure windows).
  save: boolean
      Flag whether to save figures as graphic (png) files.
  title: string
      Basic string of each figure title. Will be complemented with the absspec
       tag corresponging to the respective figure.
  casename: string
      Basic string for graphics file names. Graphics are saved in figures/ with
       full name constructed from a plot type indentifier, this casename and the
       abspsec tag of the corresponding species.

  Returns
  -------
  no returns
  '''
  if ende is None: ende=axs1.shape[1]-1
  for i in N.arange(start,ende+1):                     
    if ((axs1[0,i,:,:]+axs2[0,i,:,:]).max()!=0.):
        if newfig: p.figure()
        else: p.clf()
        p.semilogy(freq*1e-9,0.5*abs(axs1[0,i,:,:]-axs2[0,i,:,:])/
                             (axs1[0,i,:,:]+axs2[0,i,:,:])+1e-16)
        p.title('%s %s' %(title,absspec[i]))
        p.xlabel('Frequency [GHz]')
        p.ylabel('relative difference [-]')
        p.tight_layout()
        if save:
          p.savefig('figures/drel_logabs_%s_%s.png' %(casename,absspec[i]))
        if newfig: p.figure()
        else: p.clf()
        p.plot(freq*1e-9,0.5*(axs1[0,i,:,:]-axs2[0,i,:,:])/
                         (axs1[0,i,:,:]+axs2[0,i,:,:]))
        p.title('%s %s' %(title,absspec[i]))
        p.xlabel('Frequency [GHz]')
        p.ylabel('relative difference [-]')
        p.tight_layout()
        if save:
          p.savefig('figures/drel_lin_%s_%s.png' %(casename,absspec[i]))
        if newfig: p.figure()
        else: p.clf()
        p.plot(freq*1e-9,axs1[0,i,:,:]-axs2[0,i,:,:])
        p.xlabel('Frequency [GHz]')
        p.title('%s %s' %(title,absspec[i]))
        p.ylabel('absolute difference [m2]')
        p.tight_layout()
        if save:
          p.savefig('figures/dabs_lin_%s_%s.png' %(casename,absspec[i]))
        if newfig: p.figure()
        else: p.clf()
        p.semilogy(freq*1e-9,abs(axs1[0,i,:,:]-axs2[0,i,:,:])+1e-50)
        p.title('%s %s' %(title,absspec[i]))
        p.xlabel('Frequency [GHz]')
        p.ylabel('absolute difference [m2]')
        p.tight_layout()
        if save:
          p.savefig('figures/dabs_logabs_%s_%s.png' %(casename,absspec[i]))
    else:
        print('No abs xs !=0 (no lines? vmr=0?) for %s' %absspec[i])
    p.show()


def first_substring(strings, substring):
  '''
  Derives the first occurence of substring within a list of strings.
  
  Parameters
  ----------
  strings: list of strings
      list with strings that shall be searched
  substring: string
      the string whose first occurence in strings shall be determined (can be a
       substring of the entries in strings)
      
  Returns
  -------
  index: integer
      the index of the first entry in strings that contains substring
  '''
  return (i for i, string in enumerate(strings) if substring in string).next()

if __name__ == '__main__':
  fa3h='ARTS3-from-HITRAN08_LUT_WithFascode.xml'
  fa4h='ARTS4-from-HITRAN08_LUT_WithFascode.xml'
  fa4p='Perrin_LUT_WithFascode.xml'
  fa3h2='ARTS3-from-HITRAN08_LUT_NoFascode.xml'
  fa4h2='ARTS4-from-HITRAN08_LUT_NoFascode.xml'
  fa4p2='Perrin_LUT_NoFascode.xml'

  absspec,freq,pres,axs_a3h,axs_a4h,axs_a4p = get_data(fa3h,fa4h,fa4p)
  absspec2,freq2,pres2,axs_a3h2,axs_a4h2,axs_a4p2 = get_data(fa3h2,fa4h2,fa4p2)

  fast_stats(absspec,axs_a3h, axs_a4h, axs_a4p)
  fast_stats(absspec2,axs_a3h2, axs_a4h2, axs_a4p2)


