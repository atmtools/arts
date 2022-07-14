# -*- coding: utf-8 -*-

import os
import pyarts
import numpy as np

def all_cia_files(path):
    out = []
    for fil in os.listdir(path):
        fil = os.path.abspath(os.path.join(path, fil))
        if fil.endswith('.cia'):
            out.append(fil)
    return out


class header:
    def __init__(self, line):
        assert len(line) >= 20+10+10+7+7+10+6+27+3, \
            f"This is not a header:\n'{line}'"

        self.chem_sym = line[:20]
        self.low_lim = line[20:20+10]
        self.upp_lim = line[20+10:20+10+10]
        self.num = line[20+10+10:20+10+10+7]
        self.T = line[20+10+10+7:20+10+10+7+7]
        self.max_cia = line[20+10+10+7+7:20+10+10+7+7+10]
        self.res = line[20+10+10+7+7+10:20+10+10+7+7+10+6]
        self.comment = line[20+10+10+7+7+10+6:20+10+10+7+7+10+6+27]
        self.ref = line[20+10+10+7+7+10+6+27:20+10+10+7+7+10+6+27+3]

        self.chem_sym = self.chem_sym.strip()
        self.comment = self.comment.strip()

        self.low_lim = float(self.low_lim)
        self.upp_lim = float(self.upp_lim)
        self.T = float(self.T)
        self.max_cia = float(self.max_cia)
        self.res = float(self.res)

        self.num = int(self.num)
        self.ref = int(self.ref)

    def __repr__(self):
        return f"{self.chem_sym} {self.low_lim}-{self.upp_lim} cm-1 at {self.T} K"


class data:
    def __init__(self, num):
        self.i = 0
        self.n = num
        self.v = np.zeros(shape=num) * np.nan
        self.k = np.zeros(shape=num) * np.nan

    def __call__(self, line):
        s = line.split()
        self.v[self.i] = float(s[0])
        self.k[self.i] = float(s[1])
        self.i += 1
        return self.i < self.n

    def __repr__(self):
        return f"{self.n} datapoints"


def parse_cia_file(path):
    a = open(path, 'r')

    out = []

    line = a.readline()
    while line:
        h = header(line)
        d = data(h.num)
        while d(a.readline()):
            pass
        out.append((h, d))
        line = a.readline()

    a.close()
    return out


def parse_cia_files(path):
    out = []
    for fil in all_cia_files(path):
        out.append(parse_cia_file(fil))
    return out


class cia_record:
    def __init__(self, h, d):
        self.t = [h.T]
        self.k = [d.k]
        self.v = d.v

        self.h = [h]

    def __call__(self, h, d):
        if len(d.v) != len(self.v) or np.any(d.v != self.v):
            return False
        self.t.append(h.T)
        self.k.append(d.k)
        self.h.append(h)
        return True

    def any(self): return len(self.t)

    def fix(self):
        x = np.argsort(self.t)
        self.t = np.array(self.t)[x]
        self.k = np.array(self.k)[x]
        self.h = np.array(self.h)[x]
        return self

    def try_merge(self, other):
        count = 0
        for n in self.v:
            count += n in other.v

        # Must be complete overlap
        if count != len(self.v) and count != len(other.v):
            return other, False

        print(f"Merge by extending zeros {self.h} and {other.h}")

        k = np.zeros((len(self.t)+len(other.t),
                     max(len(self.v), len(other.v))))
        v = self.v if count == len(self.v) else other.v

        for i in range(len(self.t)):
            for j in range(len(self.v)):
                if self.v[j] in v:
                    k[i, j] = self.k[i, j]
        for i in range(len(other.t)):
            for j in range(len(other.v)):
                if other.v[j] in v:
                    k[i, j] = other.k[i, j]

        self.v = self.v if len(self.v) > len(other.v) else other.v
        self.k = k
        self.t = np.append(self.t, other.t)
        self.h = np.append(self.h, other.h)

        return self.fix(), True

    def try_interp_extrap_merge(self, other):
        sm = np.median(self.v)
        om = np.median(other.v)

        overlap = (sm > other.v[0] and sm < other.v[-1]
                   ) or (om > self.v[0] and om < self.v[-1])

        if not overlap:
            return other, False

        print(
            f"Merge by interpolation and zero extension {self.h} and {other.h}")

        x = np.unique(np.append(self.v, other.v))

        k = []
        for i in range(len(self.t)):
            k.append(np.interp(x, self.v, self.k[i], left=0, right=0))
        for i in range(len(other.t)):
            k.append(np.interp(x, other.v, other.k[i], left=0, right=0))

        self.v = x
        self.k = k
        self.t = np.append(self.t, other.t)
        self.h = np.append(self.h, other.h)

        return self.fix(), True

    def __repr__(self):
        return f"{self.h[0].chem_sym} ({len(self.t)}x{len(self.v)}) {self.v[0]}-{self.v[-1]} cm-1; {self.t[0]}-{self.t[-1]} K"


class false_cia_record:
    def __init__(self): pass
    def __call__(*args): return False
    def any(self): return False


def sort_cia_data(datalist):
    out = []
    data = false_cia_record()
    for x in datalist:
        if not data(x[0], x[1]):
            if data.any():
                out.append(data.fix())
            data = cia_record(x[0], x[1])
    if data.any():
        out.append(data.fix())
    return out


def count(data):
    out = {}
    for spec in data:
        for f in spec:
            out[f.spec] = out.get(f.spec, 0) + 1
    return out


def fix_missing(data, first=True):
    n = len(data)

    out = [data.pop(0)]
    while len(data):
        next = data.pop(0)
        x, v = out[-1].try_merge(next)

        if not v:
            x, v = out[-1].try_interp_extrap_merge(next)

        if v:
            out[-1] = x
        else:
            out.append(x)

    if first or n != len(out):
        out = fix_missing(out[::-1], False)[::-1]
    return out


def remove_negatives(data):
    had = 0
    for x in data:
        count = (x.k < 0).sum()
        if count:
            f"Remove {count} negative values {x.h}"
        x.k[x.k < 0] = 0
    return data

def hitran2arts(path):
    """ Reads all Hitran .cia files in a folder and tries to convert them to
    Arts format
    
    When the Hitran file data contains two datasets and one data set has a
    pure subset of the frequency grid of the other, these datasets are merged
    with the missing data set to zero.  A print() expression is called
    informing that this has happened, starting with the words
    "Merge by extending zeros" and finishing with a copy of the two datasets'
    names/headers
    
    When the Hitran file data contains two data sets that clearly overlap but
    where one data set is not a complete subset of the other, a new frequency
    grid as the combination of these two datasets' grids is formed.  The 
    two sets of data are then interpolated onto the new frequency grid by
    linear means.  Extrapolation is not allowed so any data outside of the
    range of the old frequency grid is set to zero.  A print() expression is
    called informing that this has happened, starting with the words
    "Merge by interpolation and zero extension" and finishing with a copy of
    the two datasets' names/headers
    
    If the Hitran data contains negative values, these are set to zero.
    print() informs how many negative values are removed.
    
    If generating the CIARecord fails, a warning about this is printed.

    Parameters
    ----------
    path : Path line
        A path to a folder with Hitran cia data

    Returns
    -------
    aoc : dict
        A list of all CIARecord that could be generated

    """
    x = [sort_cia_data(x) for x in parse_cia_files(path)]
    
    for i in range(len(x)):
        x[i] = remove_negatives(fix_missing(x[i]))

    aoc = {}
    for specs in x:
        out = []

        spec = specs[0].h[0].chem_sym.replace("Air", "Bath")
        for data in specs:
            out.append(pyarts.arts.GriddedField2([pyarts.arts.convert.kaycm2freq(data.v), data.t],
                                                 100**-5 * data.k.T, ["Frequency", "Temperature"]))

        try:
            aoc[spec.replace('-', '-CIA-')] = pyarts.arts.CIARecord(
                out, spec.split('-')[0], spec.split('-')[1])
        except:
            print("Fail for", spec)
    return aoc