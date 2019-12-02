# -*- coding: utf-8 -*-
"""
Implementation of classes to handle various ARTS internal structures.

"""


import arts.constants as constants
import arts.spectroscopy as spectroscopy
import arts
    
import numpy as np
import scipy.interpolate as _ip
from scipy.special import wofz as _Faddeeva_
from fractions import Fraction as _R
from numpy.polynomial import Polynomial as _P

import os

__all__ = ['ARTSCAT5',
           'Rational',
           'PressureBroadening',
           'LineMixing',
           'PartitionFunctions',
           'LineFunctionsData',
           'read_hitran_online'
           ]

def read_hitran_online(hitran_file, fmin=0, fmax=1e9999,
                       hit_tmp='HITRAN2012.par', remove_hit_tmp=True,
                       reset_qn=True):
    """ Reads catalog from HITRAN online
    
    This function is meant for specific input, so failure is an option.  This
    code might also break in the future should HITRAN update their format again
    
    The format has to be [.par line, qns', qns''] in that order
    
    There is a temporary file generated.  This must be named so that the ARTS
    function abs_linesReadFromHitran can read the file.  After the temporary
    file has been read by ARTS, the catalog is accessed for each line to set
    their quantum numbers of qns' and qns''.  Note that these are added as 
    additions to existing numbers.  Caveats follow.
    
    The following qunatum numbers are ignored at present:
        kronigParity
        
        ElecStateLabel
        
    The following are reinterpreted for ARTS:
        v to v1
        
        nuclearSpinRef to F
        
    The following are given values despite having none:
        parity=- is set as -1
        
        parity=+ is set as +1
        
    All other variables are written to verbatim.  This means that if ARTS does
    not have the quantum number keys prescribed, the catalog that is produced
    will not be possible to read with ARTS.  Take care!
    
    Input:
        hitran_file:
            a file generated from hitran.org
        
        fmin:
            fmin in abs_linesReadFromHitran
        
        fmax:
            fmax in abs_linesReadFromHitran
        
        hit_tmp:
            temporary filename for abs_linesReadFromHitran
        
        remove_hit_tmp:
            Flag to remove or not the file hit_tmp using os.remove
    
    Output:
        ARTSCAT5 catalog of lines
    """
    assert not hitran_file == hit_tmp, "Must have separate files"

    # setup a quiet workspace
    from arts.workspace import Workspace
    arts = Workspace(verbosity=0)

    # read a file from hitran of format [par, qn upper, qn lower]
    f = open(hitran_file, 'r')
    par_data = ''
    up_data = []
    lo_data = []
    line = f.readline().split('\t')
    while len(line) > 2:
        par_data += line[0]
        up_data.append(line[1].split(';'))
        lo_data.append(line[2].split(';'))
        lo_data[-1][-1].replace('\n', '')
        line = f.readline().split('\t')
        if len(line) > 2:
            par_data += '\n'
    f.close()

    # Create a temporary file and read it into arts
    f = open(hit_tmp, 'w')
    f.write(par_data)
    f.close()
    arts.abs_linesReadFromHitran(filename=hit_tmp, fmin=float(fmin),
                                 fmax=float(fmax))

    # delete previous file
    if remove_hit_tmp:
        os.remove(hit_tmp)

    # replace quantum numbers by
    arts.abs_linesNewestArtscatFromLegacyCatalog()
    cat = arts.abs_lines.value.as_ARTSCAT5()

    for i in range(len(cat)):
        if reset_qn:
            cat._dictionaries[i]['QN'] = arts.catalogues.QuantumNumberRecord(arts.catalogues.QuantumNumbers(''), arts.catalogues.QuantumNumbers(''))
        for qn in up_data[i]:
            key, data = qn.split('=')
            if 'ElecStateLabel' in key:
                pass
            elif 'nuclearSpinRef' in key:
                cat.quantumnumbers(i)['UP']['F'] = Rational(data)
            elif 'kronigParity' in key:
                pass
            elif 'parity' in key:
                if data == '-':
                    cat.quantumnumbers(i)['UP'][key] = Rational(-1)
                elif data == '+':
                    cat.quantumnumbers(i)['UP'][key] = Rational(+1)
            elif 'v' == key:
                cat.quantumnumbers(i)['UP']['v1'] = Rational(data)
            else:
                cat.quantumnumbers(i)['UP'][key] = Rational(data)
        for qn in lo_data[i]:
            key, data = qn.split('=')
            if 'ElecStateLabel' in key:
                pass
            elif 'nuclearSpinRef' in qn:
                cat.quantumnumbers(i)['LO']['F'] = Rational(data)
            elif 'kronigParity' in key:
                pass
            elif 'parity' in key:
                if data == '-':
                    cat.quantumnumbers(i)['LO'][key] = Rational(-1)
                elif data == '+':
                    cat.quantumnumbers(i)['LO'][key] = Rational(+1)
            elif 'v' == key:
                cat.quantumnumbers(i)['LO']['v1'] = Rational(data)
            else:
                cat.quantumnumbers(i)['LO'][key] = Rational(data)
    return cat


class LineFunctionsData:
    def __init__(self):
        self.LS = None
        self.LM = None
        self.species = None
        self.data = None
        
    @staticmethod
    def len_of_key(key):
        if key in ["LM_AER"]:
            return 12
        elif key in ["#"]:
            return 0
        elif key in ["T0"]:
            return 1
        elif key in ["T1", "T3", "T5"]:
            return 2
        elif key in ["T2", "T4"]:
            return 3
    
    def fill_data(self, array, ndata, start=0):
        pos = 1*start
        i = 0
        while i < self.species:
            self.data[i]['spec'] = array[pos]
            pos += 1
            j = 0
            while j < ndata:
                self.data[i]['data'][j]['key'] = array[pos]
                pos += 1
                for k in range(self.len_of_key(self.data[i]['data'][j]['key'])):
                    self.data[i]['data'][j]['val'].append(array[pos])
                    pos += 1
                j += 1
            i += 1
        return pos

    def read_as_part_of_artscat5(self, array, i):
        self.LS = array[i+0]
        self.LM = array[i+1]
        self.species = int(array[i+2])
        j = i + 3
        return self.fill_data(array, self.set_data_shape(), j)
    
    def len_of_data(self):
        if self.LS in ['DP']:
            shape_len = 0
        elif self.LS in ['LP', 'VP']:
            shape_len = 2
        elif self.LS in ['SDVP']:
            shape_len = 4
        elif self.LS in ['HTP']:
            shape_len = 6
        else:
            raise RuntimeError(f'Unknown LS value: {self.LS}')
        if self.LM in ['#']:
            mixing_len = 0
        elif self.LM in ['LM1', 'INT', 'ConstG']:
            mixing_len = 1
        elif self.LM in ['LM2']:
            mixing_len = 3
        else:
            raise RuntimeError(f'Unknown LM value: {self.LM}')
        return mixing_len + shape_len
        
    def set_data_shape(self):
        ndata = self.len_of_data()
        
        self.data = {}
        i = 0
        while i < self.species:
            self.data[i] = {'spec': None, 'data': {}}
            j = 0
            while j < ndata:
                self.data[i]['data'][j] = {'key': None, 'val': []}
                j += 1
            i += 1
        return ndata
    
    def __repr__(self):
        ndata = self.len_of_data()
        st = ''
        st += self.LS + ' ' + self.LM + ' ' + str(self.species) + ' '
        for x in range(self.species):
            st += self.data[x]['spec'] + ' '
            for y in range(ndata):
                st += self.data[x]['data'][y]['key'] + ' '
                for z in self.data[x]['data'][y]['val']:
                    st += z + ' '
        return st
    __str__=__repr__


class ARTSCAT5:
    """Class to contain ARTSCAT entries that can be accessed and  manipulated

        Access this data as
        (N, I, F0, S0, T0, E0, A, GL, GU, PB, QN, LM)  = ARTSCAT5[line_nr],
        where N is the name of the species, I is the AFGL isotopological code,
        F0 is the central frequency, S0 is the line strength at temperature T0,
        E0 is the lower energy state, A is the einstein coefficient, GL is the
        lower population constant, GU is the upper population constant, PB
        is a dictionary of the pressurebroadening scheme, QN is a
        QuantumNumberRecord, and LM is a line-mixing dictionary.  The
        dictionaries have keys corresponding to the ARTSCAT tags.  line_nr is
        an index absolutely less than len(self)

        Note:  Must be ARTSCAT5 line type or this will leave the class data
        in disarray.

        Future tech debt 1: The reading of tagged data will fail if major tags
        ever match minor tags (ex. PB is major, N2 is minor for PB, if LM ever
        gets the PB minor tag, then the method below will fail).

        Future tech debt 2: To add version number and make the class general,
        ARTSCAT3 only have minor tag N2 for line mixing and ARTSCAT4 has AP.
        These tags are however implicit and not written.  A
    """

    _spec_ind = 0
    _iso_ind = 1
    _freq_ind = 2
    _str_ind = 3
    _t0_ind = 4
    _elow_ind = 5
    _ein_ind = 6
    _glow_ind = 7
    _gupp_ind = 8
    _pb_ind = 9
    _qn_ind = 10
    _lm_ind = 11
    _ze_ind = 12
    _lf_ind = 13
    _lsm_ind = 14

    def __init__(self, init_data=None):
        self._dictionaries = np.array([], dtype=dict)
        self._n = 0
        self.LineRecordData = {
                'freq': np.array([]),
                'afgl': np.array([], dtype='int'),
                'str': np.array([]),
                'glow': np.array([]),
                'gupp': np.array([]),
                'elow': np.array([]),
                'spec': np.array([], dtype='str'),
                'ein': np.array([]),
                't0': np.array([])}

        if init_data is None:
            return

        self.append(init_data, sort=False)

    def _append_linestr_(self, linerecord_str):
        """Takes an arts-xml catalog string and appends info to the class data
        """
        lr = linerecord_str.split()
        len_lr = len(lr)
        if len_lr == 0:
            return
        assert len_lr > 9, "Cannot recognize line data"

        self._dictionaries = np.append(self._dictionaries,
                                       {"QN": QuantumNumberRecord(),
                                        "PB": {"Type": None,
                                               "Data": np.array([])},
                                        "LM": {"Type": None,
                                               "Data": np.array([])},
                                        "LF": LineFunctionsData(),
                                        "ZE": None,
                                        "LSM": {}})

        spec = lr[1].split('-')
        self.LineRecordData['spec'] = np.append(self.LineRecordData['spec'],
                                                spec[self._spec_ind])
        self.LineRecordData['afgl'] = np.append(self.LineRecordData['afgl'],
                                                int(spec[self._iso_ind]))
        self.LineRecordData['freq'] = np.append(self.LineRecordData['freq'],
                                                float(lr[self._freq_ind]))
        self.LineRecordData['str'] = np.append(self.LineRecordData['str'],
                                               float(lr[self._str_ind]))
        self.LineRecordData['t0'] = np.append(self.LineRecordData['t0'],
                                              float(lr[self._t0_ind]))
        self.LineRecordData['elow'] = np.append(self.LineRecordData['elow'],
                                                float(lr[self._elow_ind]))
        self.LineRecordData['ein'] = np.append(self.LineRecordData['ein'],
                                               float(lr[self._ein_ind]))
        self.LineRecordData['glow'] = np.append(self.LineRecordData['glow'],
                                                float(lr[self._glow_ind]))
        self.LineRecordData['gupp'] = np.append(self.LineRecordData['gupp'],
                                                float(lr[self._gupp_ind]))
        self._n += 1

        key = lr[9]
        i = 10
        qnr = ''
        ze = {"POL": None}
        while i < len_lr:
            this = lr[i]
            if this in ['QN', 'PB', 'LM', 'ZE', 'LSM', 'LF']:
                key = this
            elif key == 'QN':
                qnr += ' ' + this
            elif key == 'ZE':
                ze = {"POL": lr[i], "GU": float(lr[i+1]),
                      "GL": float(lr[i+2])}
                i += 2
            elif key == 'LSM':
                x = int(lr[i])
                i += 1
                for nothing in range(x):
                    self._dictionaries[-1]['LSM'][lr[i]] = lr[i+1]
                    i += 2
            elif key == 'LF':
                i=self._dictionaries[-1]['LF'].read_as_part_of_artscat5(lr, i)-1
            else:
                try:
                    self._dictionaries[-1][key]["Data"] = \
                        np.append(self._dictionaries[-1][key]["Data"],
                                  float(this))
                except:
                    self._dictionaries[-1][key]["Type"] = this
            i += 1
        self._dictionaries[-1]['QN'] = QuantumNumberRecord.from_str(qnr)
        self._dictionaries[-1]['LM'] = LineMixing(self._dictionaries[-1]['LM'])
        self._dictionaries[-1]['PB'] = \
            PressureBroadening(self._dictionaries[-1]['PB'])
        self._dictionaries[-1]['ZE'] = ze

    def _append_line_(self, line):
        """Appends a line from data
        """
        self.LineRecordData['spec'] = np.append(self.LineRecordData['spec'],
                                                str(line[self._spec_ind]))
        self.LineRecordData['afgl'] = np.append(self.LineRecordData['afgl'],
                                                int(line[self._iso_ind]))
        self.LineRecordData['freq'] = np.append(self.LineRecordData['freq'],
                                                line[self._freq_ind])
        self.LineRecordData['str'] = np.append(self.LineRecordData['str'],
                                               line[self._str_ind])
        self.LineRecordData['t0'] = np.append(self.LineRecordData['t0'],
                                              line[self._t0_ind])
        self.LineRecordData['elow'] = np.append(self.LineRecordData['elow'],
                                                line[self._elow_ind])
        self.LineRecordData['ein'] = np.append(self.LineRecordData['ein'],
                                               line[self._ein_ind])
        self.LineRecordData['glow'] = np.append(self.LineRecordData['glow'],
                                                line[self._glow_ind])
        self.LineRecordData['gupp'] = np.append(self.LineRecordData['gupp'],
                                                line[self._gupp_ind])
        self._dictionaries = np.append(self._dictionaries,
                                       {'PB': line[self._pb_ind],
                                        'QN': line[self._qn_ind],
                                        'LM': line[self._lm_ind],
                                        'ZE': line[self._ze_ind],
                                        'LF': line[self._lf_ind],
                                        'LSM': line[self._lsm_ind]})
        self._n += 1

    @property
    def F0(self):
        return self.LineRecordData['freq']

    @F0.setter
    def F0(self, nums):
        self.LineRecordData['freq'] = nums

    @property
    def S0(self):
        return self.LineRecordData['str']

    @S0.setter
    def S0(self, nums):
        self.LineRecordData['str'] = nums

    @property
    def Species(self):
        return self.LineRecordData['spec']

    @property
    def Iso(self):
        return self.LineRecordData['afgl']

    @property
    def T0(self):
        return self.LineRecordData['t0']

    @property
    def A(self):
        return self.LineRecordData['ein']

    @property
    def g00(self):
        return self.LineRecordData['gupp']

    @property
    def g0(self):
        return self.LineRecordData['glow']

    @property
    def E0(self):
        return self.LineRecordData['elow']

    def _append_ArrayOfLineRecord_(self, array_of_linerecord):
        """Appends lines in ArrayOfLineRecord to ARTSCAT5
        """
        assert array_of_linerecord.version == 'ARTSCAT-5', "Only for ARTSCAT-5"
        for l in array_of_linerecord:
            self._append_linestr_(l)

    def _append_ARTSCAT5_(self, artscat5):
        """Appends all the lines of another artscat5 to this
        """
        for line in artscat5:
            self._append_line_(line)

    def set_testline(self, i_know_what_i_am_doing=False):
        assert(i_know_what_i_am_doing)
        self._n = 1
        self.LineRecordData = {
                        'freq': np.array([100e9], dtype='float'),
                        'afgl': np.array([626], dtype='int'),
                        'str': np.array([1], dtype='float'),
                        'glow': np.array([0], dtype='float'),
                        'gupp': np.array([3], dtype='float'),
                        'elow': np.array([0], dtype='float'),
                        'spec': np.array(['CO2'], dtype='str'),
                        'ein': np.array([1], dtype='float'),
                        't0': np.array([300], dtype='float')}
        self._dictionaries = np.array([
                {"QN": QuantumNumberRecord(as_quantumnumbers("J 1"),
                                           as_quantumnumbers("J 0")),
                 "PB": PressureBroadening([10e3, 0.8, 20e3, 0.8, 1e3,
                                           -1, -1, -1, -1, -1]),
                 "LM": LineMixing([300, 1e-10, 0.8])}])

    def append(self, other, sort=True):
        """Appends data to ARTSCAT5.  Used at initialization

        Parameters:
            other (str, ARTSCAT5, ArrayOfLineRecord, tuple): Data to append,
            Must fit with internal structures.  Easiest to guarantee if other
            is another ARTSCAT5 or an ArrayOfLineRecord containing ARTSCAT-5
            data.

            sort: Sorts the lines by frequency if True
        """
        if type(other) is str:
            self._append_linestr_(other)
        elif type(other) is ARTSCAT5:
            self._append_ARTSCAT5_(other)
        elif type(other) is tuple:  # For lines --- this easily fails
            self._append_line_(other)
        elif type(other) is ArrayOfLineRecord:
            self._append_ArrayOfLineRecord_(other)
        elif type(other) in [list, np.ndarray]:
            for x in other:
                self.append(x)
        else:
            assert False, "Unknown type"
        self._assert_sanity_()

        if sort:
            self.sort()

    def sort(self, kind='freq', ascending=True):
        """Sorts the ARTSCAT5 data by kind.  Set ascending to False for
        descending sorting

        Parameters:
            kind (str): The key to LineRecordData

            ascending (bool): True sorts ascending, False sorts descending

        Examples:
            Sort by descending frequnecy

            >>> cat = arts.xml.load('C2H2.xml').as_ARTSCAT5()
            >>> cat.LineRecordData['freq']
            array([  5.94503434e+10,   1.18899907e+11,   1.78347792e+11, ...,
                     2.25166734e+12,   2.31051492e+12,   2.36933091e+12])
            >>> cat.sort(ascending=False)
            >>> cat.LineRecordData['freq']
            array([  2.36933091e+12,   2.31051492e+12,   2.25166734e+12, ...,
                     1.78347792e+11,   1.18899907e+11,   5.94503434e+10])

            Sort by line strength

            >>> cat = arts.xml.load('C2H2.xml').as_ARTSCAT5()
            >>> cat.LineRecordData['str']
            array([  9.02281290e-21,   7.11410308e-20,   2.34380510e-19, ...,
                     4.77325112e-19,   3.56443438e-19,   2.63222798e-19])
            >>> cat.sort(kind='str')
            >>> cat.LineRecordData['str']
            array([  9.02281290e-21,   7.11410308e-20,   2.34380510e-19, ...,
                     1.09266008e-17,   1.10644138e-17,   1.10939452e-17])

        """
        assert kind in self.LineRecordData, "kind must be in LineRecordData"

        i = np.argsort(self.LineRecordData[kind])
        if not ascending:
            i = i[::-1]

        for key in self.LineRecordData:
            self.LineRecordData[key] = self.LineRecordData[key][i]
        self._dictionaries = self._dictionaries[i]

    def remove(self, spec=None, afgl=None,
               upper_limit=None, lower_limit=None, kind='freq'):
        """Removes lines not within limits of kind

        This loops over all lines in self and only keeps those fulfilling

        .. math::
            l \\leq x \\leq u,

        where l is a lower limit, u is an upper limit, and x is a parameter
        in self.LineRecordData

        If spec and/or afgl are given, then the species and the afgl code of
        the line must match as well as the limit critera in order for the line
        to be removed

        Parameters:
            spec (str): Species string

            afgl (int): AFGL isotopologue code

            upper_limit (float): value to use for upper limit [-]

            lower_limit (float): value to use for lower limit [-]

            kind (str): keyword for determining x.  Must be key in
            self.LineRecordData

        Returns:
            None: Only changes self

        Examples:
            Remove lines below 1 THz and above 1.5 THz

            >>> cat = arts.xml.load('C2H2.xml').as_ARTSCAT5()
            >>> cat
            ARTSCAT-5 with 40 lines. Species: ['C2H2']
            >>> cat.remove(lower_limit=1000e9, upper_limit=1500e9)
            >>> cat
            ARTSCAT-5 with 9 lines. Species: ['C2H2']

            Remove weak lines

            >>> cat = arts.xml.load('C2H2.xml').as_ARTSCAT5()
            >>> cat
            ARTSCAT-5 with 40 lines. Species: ['C2H2']
            >>> cat.remove(lower_limit=1e-18, kind='str')
            >>> cat
            ARTSCAT-5 with 31 lines. Species: ['C2H2']
        """
        assert upper_limit is not None or lower_limit is not None, \
            "Cannot remove lines when the limits are undeclared"
        assert kind in self.LineRecordData, "Needs kind in LineRecordData"

        if spec is not None:
            if spec not in self.LineRecordData['spec']:
                return  # Nothing to remove

        if afgl is not None:
            if afgl not in self.LineRecordData['afgl']:
                return  # Nothing to remove

        remove_these = []
        for i in range(self._n):
            if spec is not None:
                if spec != self.LineRecordData['spec'][i]:
                    continue
            if afgl is not None:
                if afgl != self.LineRecordData['afgl'][i]:
                    continue
            if lower_limit is not None:
                if self.LineRecordData[kind][i] < lower_limit:
                    remove_these.append(i)
                    continue
            if upper_limit is not None:
                if self.LineRecordData[kind][i] > upper_limit:
                    remove_these.append(i)
                    continue

        for i in remove_these[::-1]:
            self.remove_line(i)

    def __repr__(self):
        return "ARTSCAT-5 with " + str(self._n) + " lines. Species: " + \
            str(np.unique(self.LineRecordData['spec']))

    def __len__(self):
        return self._n

    def _assert_sanity_(self):
        """Helper to assert that the data is good
        """
        assert self._n == len(self.LineRecordData['freq']) and \
            self._n == len(self.LineRecordData['spec']) and    \
            self._n == len(self.LineRecordData['afgl']) and    \
            self._n == len(self.LineRecordData['str']) and     \
            self._n == len(self.LineRecordData['t0']) and      \
            self._n == len(self.LineRecordData['elow']) and    \
            self._n == len(self.LineRecordData['ein']) and     \
            self._n == len(self.LineRecordData['glow']) and    \
            self._n == len(self.LineRecordData['gupp']) and    \
            self._n == len(self._dictionaries),                \
            self._error_in_length_message_()

    def __getitem__(self, index):
        """Returns a single line as tuple --- TODO: create LineRecord class?
        """
        return (self.LineRecordData['spec'][index],
                self.LineRecordData['afgl'][index],
                self.LineRecordData['freq'][index],
                self.LineRecordData['str'][index],
                self.LineRecordData['t0'][index],
                self.LineRecordData['elow'][index],
                self.LineRecordData['ein'][index],
                self.LineRecordData['glow'][index],
                self.LineRecordData['gupp'][index],
                self.pressurebroadening(index),
                self.quantumnumbers(index),
                self.linemixing(index),
                self.zeemandata(index),
                self.linefunctionsdata(index),
                self.lineshapemodifiers(index))

    def get_arts_str(self, index):
        """Returns the arts-xml catalog string for line at index
        """
        l = self[index]
        s = '@ ' + l[self._spec_ind] + '-' + str(l[self._iso_ind])
        s += ' ' + str(l[self._freq_ind])
        s += ' ' + str(l[self._str_ind])
        s += ' ' + str(l[self._t0_ind])
        s += ' ' + str(l[self._elow_ind])
        s += ' ' + str(l[self._ein_ind])
        s += ' ' + str(l[self._glow_ind])
        s += ' ' + str(l[self._gupp_ind])
        text = str(self.pressurebroadening(index))
        if len(text) > 0:
            s += ' PB ' + text
        text = str(self.quantumnumbers(index))
        if len(text) > 0:
            s += ' QN ' + text
        text = str(self.linemixing(index))
        if len(text) > 0:
            s += ' LM ' + text
        if self.zeemandata(index)['POL']:
            s += ' ZE ' + str(self.zeemandata(index)['POL']) + ' '
            s += str(self.zeemandata(index)['GU']) + ' '
            s += str(self.zeemandata(index)['GL'])
        text = str(self.linefunctionsdata(index))
        if len(text) > 0:
            s += ' LF ' + text

        if len(self.lineshapemodifiers(index)):
            s += ' LSM ' + str(len(self.lineshapemodifiers(index)))
            for i in self.lineshapemodifiers(index):
                s += ' ' + i + ' ' + str(self.lineshapemodifiers(index)[i])
        return s

    def pressurebroadening(self, index):
        """Return pressure broadening entries for line at index
        """
        return self._dictionaries[index]['PB']

    def quantumnumbers(self, index):
        """Return quantum number entries for line at index
        """
        return self._dictionaries[index]['QN']

    def linemixing(self, index):
        """Return line mixing entries for line at index
        """
        return self._dictionaries[index]['LM']

    def zeemandata(self, index):
        """Return Zeeman entries for line at index
        """
        return self._dictionaries[index]['ZE']

    def linefunctionsdata(self, index):
        """Return line function entries for line at index
        """
        return self._dictionaries[index]['LF']

    def lineshapemodifiers(self, index):
        """Return line mixing entries for line at index
        """
        return self._dictionaries[index]['LSM']

    def _error_in_length_message(self):
        return "Mis-matching length of vectors/lists storing line information"

    def as_ArrayOfLineRecord(self, index=None):
        """Turns ARTSCAT5 into ArrayOfLineRecord that can be stored to file
        """
        out = []
        if index is None:
            for i in range(self._n):
                out.append(self.get_arts_str(i))
        else:
            out.append(self.get_arts_str(index))
        if out == []:
            return ArrayOfLineRecord(data=[''], version='ARTSCAT-5')
        else:
            return ArrayOfLineRecord(data=out, version='ARTSCAT-5')

    def changeForQN(self, kind='change', qid=None,
                    spec=None, afgl=None, qns=None, information=None):
        """Change information of a line according to identifiers

        Input:
            kind (str): kind of operation to be applied, either 'change' for
            overwriting, 'add' for addition (+=), 'sub' for subtraction (-=),
            'remove' to remove matching lines, or 'keep' to only keep matching
            lines

            qid (QuantumIdentifier): Identifier to the transition or energy
            level as defined in ARTS

            spec (str or NoneType):  Name of species for which the operation
            applies.  None means for all species.  Must be None if qid is given

            afgl (int or NoneType):  AFGL isotopologue integer for which the
            operation applies.  None means for all isotopologue.  Must be None
            if qid is given

            qns (dict, None, QuantumNumberRecord, QuantumNumbers):  The quantum
            numbers for which the operation applies.  None means all quantum
            numbers.  Can be level or transition.  Must be None if qid is given

            information (dict or NoneType):  None for kind 'remove'. dict
            otherwise.  Keys in ARTSCAT5.LineRecordData for non-dictionaries.
            Use 'QN' for quantum numbers, 'LM' for line mixing, and 'PB' for
            pressure-broadening.  If level QN-key, the data is applied for both
            levels if they match (only for 'QN'-data)

        Output:
            None, only changes the class instance itself

        Examples:
            Add S = 1 to both levels quantum numbers by adding information to
            all lines

            >>> cat = arts.xml.load('O2.xml').as_ARTSCAT5()
            >>> cat.quantumnumbers(0)
            UP v1 0 J 32 F 61/2 N 32 LO v1 0 J 32 F 59/2 N 32
            >>> cat.changeForQN(information={'QN': {'S': 1}}, kind='add')
            >>> cat.quantumnumbers(0)
            UP S 1 v1 0 J 32 F 61/2 N 32 LO S 1 v1 0 J 32 F 59/2 N 32

            Remove all lines not belonging to a specific isotopologue and band
            by giving the band quantum numbers

            >>> cat = arts.xml.load('O2.xml').as_ARTSCAT5()
            >>> cat
            ARTSCAT-5 with 6079 lines. Species: ['O2']
            >>> cat.changeForQN(kind='keep', afgl=66, qns={'LO': {'v1': 0},
                                                           'UP': {'v1': 0}})
            >>> cat
            ARTSCAT-5 with 187 lines. Species: ['O2']

            Change the frequency of the 119 GHz line to 3000 THz by giving a
            full and unique quantum number match

            >>> cat = arts.xml.load('O2.xml').as_ARTSCAT5()
            >>> cat.sort()
            >>> cat.LineRecordData['freq']
            array([  9.00e+03,   2.35e+04,   4.01e+04, ...,   2.99e+12,
                     2.99e+12,   2.99e+12])
            >>> cat.changeForQN(afgl=66, qns={'LO': {'v1': 0, 'J': 0, 'N': 1},
                                              'UP': {'v1': 0, 'J': 1, 'N': 1}},
                                information={'freq': 3000e9})
            >>> cat.sort()
            >>> cat.LineRecordData['freq']
            array([  9.00e+03,   2.35e+04,   4.01e+04, ...,   2.99e+12,
                     2.99e+12,   3.00e+12])
        """
        if qid is not None:
            assert spec is None and afgl is None and qns is None, \
                "Only one of qid or spec, afgl, and qns combinations allowed"
            spec = qid.species
            afgl = qid.afgl
            qns = as_quantumnumbers(qns)
        else:
            qns = as_quantumnumbers(qns)

        if kind == 'remove':
            assert information is None, "information not None for 'remove'"
            remove_these = []
            remove = True
            change = False
            add = False
            sub = False
            keep = False
        elif kind == 'change':
            assert type(information) is dict, "information is not dictionary"
            remove = False
            change = True
            add = False
            sub = False
            keep = False
        elif kind == 'add':
            assert type(information) is dict, "information is not dictionary"
            remove = False
            change = False
            add = True
            sub = False
            keep = False
        elif kind == 'sub':
            assert type(information) is dict, "information is not dictionary"
            remove = False
            change = False
            add = False
            sub = True
            keep = False
        elif kind == 'keep':
            assert information is None, "information not None for 'keep'"
            remove_these = []
            remove = False
            change = False
            add = False
            sub = False
            keep = True
        else:
            raise RuntimeError(f'Invalid value for kind: {kind}')
        assert remove or change or add or keep or sub, "Invalid kind"

        # Check if the quantum number information is for level or transition
        if type(qns) is QuantumNumberRecord:
            for_transitions = True
        else:
            for_transitions = False

        if information is not None:
            for key in information:
                assert key in self.LineRecordData or \
                       key in ['QN', 'LM', 'PB'], \
                       "Unrecognized key"

        # Looping over all the line data
        for i in range(self._n):

            # If spec is None, then all species, otherwise this should match
            if spec is not None:
                if not spec == self.LineRecordData['spec'][i]:
                    if keep:
                        remove_these.append(i)
                    continue

            # If afgl is None, then all isotopes, otherwise this should match
            if afgl is not None:
                if not afgl == self.LineRecordData['afgl'][i]:
                    if keep:
                        remove_these.append(i)
                    continue

            # Test which levels match and which do not --- partial matching
            test = self.quantumnumbers(i) >= qns
            if for_transitions:
                test = [test]  # To let all and any work

            # Append lines to remove later (so indexing is not messed up)
            if remove and all(test):
                remove_these.append(i)
                continue
            elif keep and not all(test):
                remove_these.append(i)
                continue
            elif keep or remove:
                continue

            # Useless to continue past this point if nothing matches
            if not all(test) and for_transitions:
                continue
            elif not any(test):
                continue

            # There should only be matches remaining but only QN info is level-
            # based so all other infromation must be perfect match
            for info_key in information:
                info = information[info_key]
                if info_key == 'QN':
                    if for_transitions:
                        if change:
                            self._dictionaries[i]['QN'] = info
                        elif add:
                            self._dictionaries[i]['QN'] += info
                        elif sub:
                            self._dictionaries[i]['QN'] -= info
                        else:
                            assert False, "Programmer error?"
                    else:
                        if test[0]:
                            if change:
                                self._dictionaries[i]['QN']['UP'] = info
                            elif add:
                                self._dictionaries[i]['QN']['UP'] += info
                            elif sub:
                                self._dictionaries[i]['QN']['UP'] -= info
                            else:
                                assert False, "Programmer error?"
                        if test[1]:
                            if change:
                                self._dictionaries[i]['QN']['LO'] = info
                            elif add:
                                self._dictionaries[i]['QN']['LO'] += info
                            elif sub:
                                self._dictionaries[i]['QN']['LO'] -= info
                            else:
                                assert False, "Programmer error?"
                elif info_key in ['PB', 'LM']:
                    if not all(test):
                        continue
                    if change:
                        self._dictionaries[i][info_key] = info
                    elif add:
                        assert info.kind == \
                            self._dictionaries[i][info_key].kind, \
                            "Can only add to matching type"
                        self._dictionaries[i][info_key].data += info
                    elif sub:
                        assert info.kind == \
                            self._dictionaries[i][info_key].kind, \
                            "Can only sub from matching type"
                        self._dictionaries[i][info_key].data -= info
                    else:
                        assert False, "Programmer error?"
                else:
                    if not all(test):
                        continue
                    if change:
                        self.LineRecordData[info_key][i] = info
                    elif add:
                        self.LineRecordData[info_key][i] += info
                    elif sub:
                        self.LineRecordData[info_key][i] -= info
                    else:
                        assert False, "Programmer error?"

        # Again, to not get into index problems, this loop is reversed
        if remove or keep:
            for i in remove_these[::-1]:
                self.remove_line(i)

    def remove_line(self, index):
        """Remove line at index from line record
        """
        for key in self.LineRecordData:
            t1 = self.LineRecordData[key][:index]
            t2 = self.LineRecordData[key][(index+1):]
            self.LineRecordData[key] = np.append(t1, t2)

        _t = self._dictionaries
        t1 = _t[:index]
        t2 = _t[(index+1):]
        self._dictionaries = np.append(t1, t2)

        self._n -= 1
        self._assert_sanity_()

    def cross_section(self, temperature=None, pressure=None,
                      vmrs=None, mass=None, isotopologue_ratios=None,
                      partition_functions=None, f=None):
        """Provides an estimation of the cross-section in the provided
        frequency range

        Computes the following estimate (summing over all lines):

        .. math::
            \\sigma(f) = \\sum_{k=0}^{k=n-1}
            r_k S_{0, k}(T_0) K_1 K_2 \\frac{Q(T_0)}{Q(T)}
            \\frac{1 + g_k \\; p^2 + iy_k \\; p }{\\gamma_{D,k}\\sqrt{\\pi}}
            \\; F\\left(\\frac{f - f_{0,k} -  \\Delta f_k \\; p^2 -
            \\delta f_kp + i\\gamma_{p,k}p} {\\gamma_{D,k}}\\right),

        where there are n lines,
        r is the isotopologue ratio,
        S_0 is the line strength,
        K_1 is the boltzman level statistics,
        K_2 is the stimulated emission,
        Q is the partition sum, G is the second
        order line mixing coefficient,
        p is pressure,
        Y is the first order line mixing coefficient,
        f_0 is the line frequency,
        Delta-f is the second order line mixing frequency shift,
        delta-f is the first order pressure shift,
        gamma_p is the pressure broadening half width,
        gamma_D is the Doppler half width, and
        F is assumed to be the Faddeeva function

        Note 1: this is only meant for quick-and-dirty estimates.  If data is
        lacking, very simplistic assumptions are made to complete the
        calculations.
        Lacking VMR for a species assumes 1.0 vmr of the line species itself,
        lacking mass assumes dry air mass,
        lacking isotopologue ratios means assuming a ratio of unity,
        lacking partition functions means the calculations are performed at
        line temperatures,
        lacking frequency means computing 1000 frequencies from lowest
        frequency line minus its pressure broadening to the highest frequency
        line plus its pressure broadening,
        lacking pressure means computing at 1 ATM, and
        lacking temperature assumes atmospheric temperatures the same as the
        first line temperature.  If input f is None then the return
        is (f, sigma), else the return is (sigma)

        Warning: Use only as an estimation, this function is neither optimized
        and is only tested for a single species in arts-xml-data to be within
        1% of the ARTS computed value

        Parameters:
            temperature (float): Temperature [Kelvin]

            pressure (float): Pressure [Pascal]

            vmrs (dict-like): Volume mixing ratios.  See PressureBroadening for
            use [-]

            mass (dict-like): Mass of isotopologue [kg]

            isotopologue_ratios (dict-like):  Isotopologue ratios of the
            different species [-]

            partition_functions (dict-like):  Partition function estimator,
            should compute partition function by taking temperature as the only
            argument [-]

            f (ndarray): Frequency [Hz]

        Returns:
            (f, xsec) or xsec depending on f

        Examples:
            Plot cross-section making no assumptions on the atmosphere or
            species, i.e., isotopologue ratios is 1 for all isotopologue
            (will not agree with ARTS)

            >>> import matplotlib.pyplot as plt
            >>> cat = arts.xml.load('O2.xml').as_ARTSCAT5()
            >>> (f, x) = cat.cross_section()
            >>> plt.plot(f, x)

            Plot cross-sections by specifying limited information on the
            species (will agree reasonably with ARTS)

            >>> import matplotlib.pyplot as plt
            >>> cat = arts.xml.load('O2.xml').as_ARTSCAT5()
            >>> cat.changeForQN(afgl=66, kind='keep')
            >>> f, x = cat.cross_section(mass={"O2-66": 31.9898*constants.amu},
                                         isotopologue_ratios={"O2-66": 0.9953})
            >>> plt.plot(f, x)

        """
        if self._n == 0:
            if f is None:
                return 0, 0
            else:
                return np.zeros_like(f)

        if temperature is None:
            temperature = self.LineRecordData['t0'][0]

        if pressure is None:
            pressure = constants.atm

        if vmrs is None:
            vmrs = {}

        if mass is None:
            mass = {}

        if f is None:
            return_f = True
            f0 = self.pressurebroadening(0).compute_pressurebroadening_params(
                    temperature, self.LineRecordData['t0'][0],
                    pressure, vmrs)[0]
            f0 = self.LineRecordData['freq'][0] - f0
            f1 = self.pressurebroadening(-1).compute_pressurebroadening_params(
                    temperature, self.LineRecordData['t0'][-1],
                    pressure, vmrs)[0]
            f1 = self.LineRecordData['freq'][-1] - f1
            f = np.linspace(f0, f1, num=1000)
        else:
            return_f = False

        if isotopologue_ratios is None:
            isotopologue_ratios = {}

        if partition_functions is None:
            partition_functions = {}

        # Cross-section
        sigma = np.zeros_like(f)

        for i in range(self._n):
            spec_key = self.LineRecordData['spec'][i] + '-' + \
                          str(self.LineRecordData['afgl'][i])

            if spec_key in mass:
                m = mass[spec_key]
            else:
                m = constants.molar_mass_dry_air / constants.avogadro
            gamma_D = \
                spectroscopy.doppler_broadening(temperature,
                                                self.LineRecordData['freq'][i],
                                                m)
            (G, Df,
             Y) = self.linemixing(i).compute_linemixing_params(temperature)

            (gamma_p,
             delta_f) = \
                self.pressurebroadening(i).compute_pressurebroadening_params(
                temperature, self.LineRecordData['t0'][i], pressure, vmrs)

            K1 = spectroscopy.boltzmann_level(self.LineRecordData['elow'][i],
                                              temperature,
                                              self.LineRecordData['t0'][i])
            K2 = spectroscopy.stimulated_emission(
                    self.LineRecordData['freq'][i],
                    temperature,
                    self.LineRecordData['t0'][i])

            if spec_key in partition_functions:
                Q = partition_functions[spec_key]
            else:
                Q = np.ones_like

            if spec_key in isotopologue_ratios:
                r = isotopologue_ratios[spec_key]
            else:
                r = 1.0

            S = r * self.LineRecordData['str'][i] * K1 * K2 * \
                Q(self.LineRecordData['t0'][i]) / Q(temperature)

            lm = 1 + G * pressure**2 + 1j * Y * pressure
            z = (f - self.LineRecordData['freq'][i] -
                 delta_f - Df * pressure**2 + 1j * gamma_p) / gamma_D
            sigma += (S * (lm * _Faddeeva_(z) / np.sqrt(np.pi) / gamma_D)).real
        if return_f:
            return f, sigma
        else:
            return sigma

    def write_xml(self, xmlwriter, attr=None):
        """Write an ARTSCAT5 object to an ARTS XML file.
        """
        tmp = self.as_ArrayOfLineRecord()
        tmp.write_xml(xmlwriter, attr=attr)


class Rational(_R):
    """Rational number

    This is a copy of fractions.Fraction with only the __repr__ function over-
    written to match ARTS style.  That is 3/2 is represented as such rather
    than as "Fraction(3, 2)".  See original class for more information,
    limitations, and options
    """
    def __init__(self, *args):
        super(Rational, self).__init__()

    def __repr__(self):
        return str(self.numerator) + '/' + str(self.denominator)
    _R.__repr__ = __repr__


class LineMixing:
    """LineMixing data as in ARTS

    Not fully feature-complete

    Used by ARTSCAT5 to estimated ARTSCAT-5 style line mixing in cross_section
    """

    _none = None
    _first_order = "L1"
    _second_order = "L2"
    _lblrtm = "LL"
    _lblrtm_nonresonant = "NR"
    _for_band = "BB"
    _possible_kinds = [_none, _first_order, _second_order, _lblrtm,
                       _lblrtm_nonresonant, _for_band]

    def __init__(self, data=None, kind=None):
        self.data = data
        if kind is not None:
            self.kind = kind
        self._manually_changed_data = False

        self._assert_sanity_()
        self._make_data_as_in_arts_()

    def _assert_sanity_(self):
        if self._type is self._none:
            assert len(self._data) == 0, "Data available for none-type"
        elif self._type is self._first_order:
            assert len(self._data) == 3, "Data mismatching first order"
        elif self._type is self._second_order:
            assert len(self._data) == 10, "Data mismatching second order"
        elif self._type is self._lblrtm:
            assert len(self._data) == 12, "Data mismatching LBLRTM data"
        elif self._type is self._lblrtm_nonresonant:
            assert len(self._data) == 1, "Data mismatching LBLRTM data"
        elif self._type is self._for_band:
            assert len(self._data) == 1, "Data mismatching band data"
        else:
            assert False, "Cannot recognize data type at all"

    def _make_data_as_in_arts_(self):
        if self._type in [self._none, self._for_band]:
            return
        elif self._type is self._first_order:
            self._t0 = self._data[0]
            self._y0 = self._data[1]
            self._ey = self._data[2]
        elif self._type is self._second_order:
            self._t0 = self._data[6]

            self._y0 = self._data[0]
            self._y1 = self._data[1]
            self._ey = self._data[7]

            self._g0 = self._data[2]
            self._g1 = self._data[3]
            self._eg = self._data[8]

            self._f0 = self._data[4]
            self._f1 = self._data[5]
            self._ef = self._data[9]
        elif self._type is self._lblrtm:
            self._y = _ip.interp1d(self._data[:4], self._data[4:8])
            self._g = _ip.interp1d(self._data[:4], self._data[8:])
        else:
            assert False, "Unknown data type"

    def _revert_from_arts_to_data_(self):
        if self._manually_changed_data:
            return

        if self._type in [self._none, self._for_band]:
            return
        elif self._type is self._first_order:
            self._data[0] = self._t0
            self._data[1] = self._y0
            self._data[2] = self._ey
        elif self._type is self._second_order:
            self._data[6] = self._t0

            self._data[0] = self._y0
            self._data[1] = self._y1
            self._data[7] = self._ey

            self._data[2] = self._g0
            self._data[3] = self._g1
            self._data[8] = self._eg

            self._data[4] = self._f0
            self._data[5] = self._f1
            self._data[9] = self._ef
        elif self._type is self._lblrtm:
            assert all(self._y.x == self._g.x), "Mismatch between y and g"
            self._data[:4] = self._y.x
            self._data[4:8] = self._y.y
            self._data[8:] = self._g.y
        else:
            assert False, "Unknown data type"

    def __repr__(self):
        self._revert_from_arts_to_data_()
        out = ''
        if self._type is self._none:
            return "No Line-Mixing"
        elif self._type in self._possible_kinds:
            out += self._type
        else:
            assert False, "Cannot recognize kind"
        for i in self.data:
            out += ' ' + str(i)
        return out

    def __str__(self):
        self._revert_from_arts_to_data_()
        out = ''
        if self._type is self._none:
            return out
        elif self._type in self._possible_kinds:
            out += self._type
        else:
            assert False, "Cannot recognize kind"
        for i in self.data:
            out += ' ' + str(i)
        return out

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, val):
        self.data[index] = val
        self._make_data_as_in_arts_()

    @property
    def data(self):
        return self._data

    @property
    def kind(self):
        return self._type

    @kind.setter
    def kind(self, val):
        found = False
        for i in self._possible_kinds:
            if i == val:
                self._type = i
                found = True
                break
        assert found, "Cannot recognize kind"

    @data.setter
    def data(self, val):
        self._data = val
        if self._data is None:
            self._data = np.array([], dtype=float)
            self._type = self._none
        elif type(self._data) is dict:
            if self._data['Type'] is None:
                self._type = self._none
            else:
                self.kind = self._data['Type']
            self._data = self._data['Data']
        else:
            if len(self._data) == 10:
                self._type = self._second_order
            elif len(self._data) == 3:
                self._type = self._first_order
            elif len(self._data) == 12:
                self._type = self._lblrtm
            elif len(self._data) == 0:
                self._type = self._none
            else:
                assert False, "Cannot recognize data type automatically"
        self._manually_changed_data = True

    def compute_linemixing_params(self, temperature):
        """Returns the line mixing parameters for given temperature(s)

        Cross-section is found from summing all lines

        .. math::
            \\sigma(f) \\propto \\sum_{k=0}^{k=n-1}
            \\left[1 + G_k \\; p^2 + iY_k \\; p\\right] \\;
            F\\left(\\frac{f - f_{0,k} -  \\Delta f_k \\; p^2 -
            \\delta f_kp + i\\gamma_{p,k}p} {\\gamma_{D,k}}\\right),

        where k indicates line dependent variables.  This function returns
        the line mixing parameters G, Y, and Delta-f.  The non-line
        mixing parameters are gamma_D as the Doppler broadening,  gamma_p
        as the pressure broadening, f as frequency,
        f_0 as the line frequency, delta-f as the first order pressure induced
        frequency shift, and p as pressure.  The function F(···) is the
        Faddeeva function and gives the line shape.  Many scaling factors are
        ignored in the equation above...

        Note 1: that for no line mixing, this function returns all zeroes

        Developer note: the internal variables used emulates the theory for
        each type of allowed line mixing.  Thus it should be easy to extend
        this for other types and for partial derivatives

        Input:
            temperature (float or ndarray) in Kelvin

        Output:
            G(temperature), Delta-f(temperature), Y(temperature)
        """
        if self._type is self._none:
            return np.zeros_like(temperature), np.zeros_like(temperature), \
                np.zeros_like(temperature)
        elif self._type is self._for_band:
            return np.zeros_like(temperature) * np.nan, \
                np.zeros_like(temperature) * np.nan, \
                np.zeros_like(temperature) * np.nan
        elif self._type is self._lblrtm:
            return self._g(temperature), np.zeros_like(temperature), \
                self._y(temperature)
        elif self._type is self._first_order:
            return np.zeros_like(temperature), np.zeros_like(temperature), \
                self._y0 * (self._t0/temperature) ** self._ey
        elif self._type is self._lblrtm_nonresonant:
            return np.full_like(temperature, self._data[0]), \
                np.zeros_like(temperature), np.zeros_like(temperature)
        elif self._type is self._second_order:
            th = self._t0 / temperature
            return (self._g0 + self._g1 * (th - 1)) * th ** self._eg, \
                (self._f0 + self._f1 * (th - 1)) * th ** self._ef, \
                (self._y0 + self._y1 * (th - 1)) * th ** self._ey


class PressureBroadening:
    """PressureBroadening data as in ARTS

    Not fully feature-complete

    Used by ARTSCAT5 to estimated ARTSCAT-5 style pressure broadening in
    cross_section
    """

    _none = None
    _air = "N2"
    _air_and_water = "WA"
    _all_planets = "AP"
    _sd_air = "SD-AIR"
    _possible_kinds = [_none, _air, _air_and_water, _all_planets, _sd_air]

    def __init__(self, data=None, kind=None):
        self.data = data
        if kind is not None:
            self.kind = kind
        self._manually_changed_data = False

        self._assert_sanity_()
        self._make_data_as_in_arts_()

    def _assert_sanity_(self):
        if self._type is self._none:
            assert len(self._data) == 0, "Data available for none-type"
        elif self._type is self._air:
            assert len(self._data) == 10, "mismatching air broadening "
        elif self._type is self._air_and_water:
            assert len(self._data) == 9, "mismatching air and water broadening"
        elif self._type is self._all_planets:
            assert len(self._data) == 20, "mismatching all planets data"
        elif self._type is self._sd_air:
            assert len(self._data) == 8, "mismatching speed dependent air data"
        else:
            assert False, "Cannot recognize data type at all"

    def _make_data_as_in_arts_(self):
        if self._type is self._none:
            return
        elif self._type is self._air:
            self._sgam = self._data[0]
            self._sn = self._data[1]
            self._sdel = 0

            self._agam = self._data[2]
            self._an = self._data[3]
            self._adel = self._data[4]

            self._dsgam = self._data[5]
            self._dnself = self._data[6]

            self._dagam = self._data[7]
            self._dnair = self._data[8]

            self._dadel = self._data[9]
        elif self._type is self._air_and_water:
            self._sgam = self._data[0]
            self._sn = self._data[1]
            self._sdel = self._data[2]

            self._agam = self._data[3]
            self._an = self._data[4]
            self._adel = self._data[5]

            self._wgam = self._data[6]
            self._wn = self._data[7]
            self._wdel = self._data[8]
        elif self._type is self._all_planets:
            self._sgam = self._data[0]
            self._sn = self._data[7]
            self._sdel = 0
            self._gam = {'N2': self._data[1], 'O2': self._data[2],
                         'H2O': self._data[3], 'CO2': self._data[4],
                         'H2': self._data[5], 'He': self._data[6]}
            self._n = {'N2': self._data[8], 'O2': self._data[9],
                       'H2O': self._data[10], 'CO2': self._data[11],
                       'H2': self._data[12], 'He': self._data[13]}
            self._delta_f = {'N2': self._data[14], 'O2': self._data[15],
                             'H2O': self._data[16], 'CO2': self._data[17],
                             'H2': self._data[18], 'He': self._data[19]}
        else:
            assert False, "Unknown data type"

    def _revert_from_arts_to_data_(self):
        if self._manually_changed_data:
            return

        if self._type is self._none:
            return
        elif self._type is self._air:
            self._data[0] = self._sgam
            self._data[1] = self._sn

            self._data[2] = self._agam
            self._data[3] = self._an
            self._data[4] = self._adel

            self._data[5] = self._dsgam
            self._data[6] = self._dnself

            self._data[7] = self._dagam
            self._data[8] = self._dnair

            self._data[9] = self._dadel
        elif self._type is self._air_and_water:
            self._data[0] = self._sgam
            self._data[1] = self._sn
            self._data[2] = self._sdel

            self._data[3] = self._agam
            self._data[4] = self._an
            self._data[5] = self._adel

            self._data[6] = self._wgam
            self._data[7] = self._wn
            self._data[8] = self._wdel
        elif self._type is self._all_planets:
            self._data[0] = self._sgam

            self._data[1] = self._gam['N2']
            self._data[2] = self._gam['O2']
            self._data[3] = self._gam['H2O']
            self._data[4] = self._gam['CO2']
            self._data[5] = self._gam['H2']
            self._data[6] = self._gam['He']

            self._data[7] = self._sn

            self._data[8] = self._n['N2']
            self._data[9] = self._n['O2']
            self._data[10] = self._n['H2O']
            self._data[11] = self._n['CO2']
            self._data[12] = self._n['H2']
            self._data[13] = self._n['He']

            self._data[14] = self._delta_f['N2']
            self._data[15] = self._delta_f['O2']
            self._data[16] = self._delta_f['H2O']
            self._data[17] = self._delta_f['CO2']
            self._data[18] = self._delta_f['H2']
            self._data[19] = self._delta_f['He']
        else:
            assert False, "Unknown data type"

    def __repr__(self):
        self._revert_from_arts_to_data_()
        out = ''
        if self._type is self._none:
            return "No Pressure-Broadening"
        elif self._type in self._possible_kinds:
            out += self._type
        else:
            assert False, "Cannot recognize kind"
        for i in self.data:
            out += ' ' + str(i)
        return out

    def __str__(self):
        self._revert_from_arts_to_data_()
        out = ''
        if self._type is self._none:
            return out
        elif self._type in self._possible_kinds:
            out += self._type
        else:
            assert False, "Cannot recognize kind"
        for i in self.data:
            out += ' ' + str(i)
        return out

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, val):
        self.data[index] = val
        self._make_data_as_in_arts_()

    @property
    def data(self):
        return self._data

    @property
    def kind(self):
        return self._type

    @kind.setter
    def kind(self, val):
        found = False
        for i in self._possible_kinds:
            if i == val:
                self._type = i
                found = True
                break
        assert found, "Cannot recognize kind"

    @data.setter
    def data(self, val):
        self._data = val
        if self._data is None:
            self._data = np.array([], dtype=float)
            self._type = self._none
        elif type(self._data) is dict:
            if self._data['Type'] is None:
                self._type = self._none
            else:
                self.kind = self._data['Type']
            self._data = self._data['Data']
        else:
            if len(self._data) == 10:
                self._type = self._air
            elif len(self._data) == 9:
                self._type = self._air_and_water
            elif len(self._data) == 20:
                self._type = self._all_planets
            elif len(self._data) == 0:
                self._type = self._none
            else:
                assert False, "Cannot recognize data type automatically"
        self._manually_changed_data = True

    def compute_pressurebroadening_params(self, temperature, line_temperature,
                                          pressure, vmrs):
        """Computes the pressure broadening parameters for the given atmosphere

        Cross-section is found from summing all lines

        .. math::
            \\sigma(f) \\propto \\sum_{k=1}^{k=n-1}
            F\\left(\\frac{f - f_{0,k} -  \\Delta f_k \\; p^2 -
            \\delta f_kp + i\\gamma_{p,k}p} {\\gamma_{D,k}}\\right),

        where k indicates line dependent variables.  This function returns
        the pressure broadening parameters p*gamma_p and p*delta-f.  The non-
        pressure broadening parameters are gamma_D as the Doppler broadening,
        f as frequency, f_0 as the line frequency, delta-f as the first order
        pressure induced frequency shift, and p as pressure.  The function
        F(···) is the Faddeeva function and gives the line shape.  Many scaling
        factors are ignored in the equation above...

        The pressure broadening parameters are summed from the contribution of
        each individual perturber so that for i perturbers

        .. math::
            \\gamma_pp = \\sum_i \\gamma_{p,i} p_i

        and

        .. math::
            \\delta f_pp = \\sum_i \\delta f_{p,i} p_i

        Parameters:
            temperature (float or ndarray): Temperature [Kelvin]

            line_temperature (float): Line temperature [Kelvin]

            pressure (float or like temperature): Total pressure [Pascal]

            vmrs (dict):  Volume mixing ratio of atmospheric species.
            dict should be {'self': self_vmr} for 'N2', {'self': self_vmr,
            'H2O': h2o_vmr} for kind 'WA', and each species of 'AP' should be
            represented in the same manner.  When 'self' is one of the list of
            species, then vmrs['self'] should not exist.  Missing data is
            treated as a missing species.  No data at all is assumed to mean
            1.0 VMR of self (len(vmrs) == 0 must evaluate as True).  The
            internal self_vmr, h2o_vmr, etc., variables must have sime size
            as pressure or be constants

        Returns:
            p · gamma0_p, p · delta-f0
        """
        theta = line_temperature / temperature
        if len(vmrs) == 0:
            return self._sgam * theta ** self._sn * pressure, \
                self._sdel * theta ** (0.25 + 1.5 * self._sn) * pressure

        sum_vmrs = 0.0
        gamma = np.zeros_like(temperature)
        delta_f = np.zeros_like(temperature)

        if self._type is self._none:
            return np.zeros_like(temperature), np.zeros_like(temperature)
        elif self._type is self._air:
            for species in vmrs:
                if species == 'self':
                    gamma += self._sgam * theta ** self._sn * \
                        pressure * vmrs[species]
                    delta_f += self._sdel * \
                        theta ** (0.25 + 1.5 * self._sn) * \
                        pressure * vmrs[species]
                    sum_vmrs += vmrs[species]
            gamma += self._agam * theta ** self._an * \
                pressure * (1 - sum_vmrs)
            delta_f += self._adel * theta ** (0.25 + 1.5 * self._an) * \
                pressure * (1 - sum_vmrs)
        elif self._type is self._air_and_water:
            for species in vmrs:
                if species == 'self':
                    gamma += self._sgam * theta ** self._sn * \
                        pressure * vmrs[species]
                    delta_f += self._sdel * \
                        theta ** (0.25 + 1.5 * self._sn) * \
                        pressure * vmrs[species]
                    sum_vmrs += vmrs[species]
                elif species == 'H2O':
                    gamma += self._wgam * theta ** self._wn * \
                        pressure * vmrs[species]
                    delta_f += self._wdel * \
                        theta ** (0.25 + 1.5 * self._wn) * \
                        pressure * vmrs[species]
                    sum_vmrs += vmrs[species]
            gamma += self._agam * theta ** self._an * \
                pressure * (1 - sum_vmrs)
            delta_f += self._adel * theta ** (0.25 + 1.5 * self._an) * \
                pressure * (1 - sum_vmrs)
        elif self._type is self._all_planets:
            for species in vmrs:
                if species == 'self':
                    gamma += self._sgam * theta ** self._sn * \
                        pressure * vmrs[species]
                    delta_f += self._sdel * \
                        theta ** (0.25 + 1.5 * self._sn) * \
                        pressure * vmrs[species]
                    sum_vmrs += vmrs[species]
                elif species in self._gam:
                    gamma += self._gam[species] * theta ** self._n[species] * \
                        pressure * vmrs[species]
                    delta_f += self._delta_f[species] * \
                        theta ** (0.25 + 1.5 * self._n[species]) * \
                        pressure * vmrs[species]
                    sum_vmrs += vmrs[species]
            gamma /= sum_vmrs
            delta_f /= sum_vmrs
        return gamma, delta_f


class PartitionFunctions:
    """Class to compute partition functions given ARTS-like partition functions
    """
    _default_test = 296.0

    def __init__(self, init_data=None):
        self._data = {}
        self.append(init_data)
        self._assert_sanity_()

    def append(self, data):
        if type(data) is SpeciesAuxData:
            self._from_species_aux_data_(data)
        elif type(data) is dict:
            self._from_dict_(data)
        elif type(data) is PartitionFunctions:
            self.data = PartitionFunctions.data
        elif data is not None:
            assert False, "Cannot recognize the initialization data type"

    def _from_species_aux_data_(self, sad):
        assert sad.version == 2, "Must be version 2 data"
        self._from_dict_(sad._data_dict)

    def _from_dict_(self, d):
        for k in d:
            assert type(d[k]) is list, "lowest level data must be list"
            self._from_list_(d[k], k)

    def _from_list_(self, l, k):
        if l[0] == 'PART_TFIELD':
            self.data[k] = _ip.interp1d(l[1][0].grids[0], l[1][0].data)
        elif l[0] == 'PART_COEFF':
            self.data[k] = _P(l[1][0].data)
        else:
            raise RuntimeError("Unknown or not implemented " +
                               "partition_functions type encountered")

    def _assert_sanity_(self):
        assert type(self.data) is dict, "Sanity check fail, calss is wrong"

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, data):
        if type(data) is list:
            self._from_list_(data, key)
        elif type(data) in [_P, _ip.interp1d]:
            self.data[key] = data
        else:
            try:
                data(self._default_test)
                self.data[key] = data
            except:
                raise RuntimeError("Cannot determine type")

    def __iter__(self):
        return iter(self.data)

    def __contains__(self, key):
        return key in self.data

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return "partition functions for " + str(len(self)) + " species"

    def keys(self):
        return self.data.keys()
    species = keys

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, val):
        assert type(val) is dict, "new values must be dictionary type"
        self._data = val

from .catalogues import SpeciesAuxData
from .catalogues import ArrayOfLineRecord
from .catalogues import QuantumNumberRecord
from .utils import as_quantumnumbers
