# -*- coding: utf-8 -*-
"""
Implementation of classes to handle various catalogue information.

"""


import numpy as np
import scipy.sparse

from pyarts.utils.arts import return_if_arts_type

__all__ = ['ArrayOfLineRecord',
           'CIARecord',
           'GasAbsLookup',
           'LineMixingRecord',
           'QuantumIdentifier',
           'QuantumNumberRecord',
           'QuantumNumbers',
           'Sparse',
           'SpeciesAuxData',
           'SpeciesTag',
           'PropagationMatrix',
           'StokesVector',
           'Ppath',
           'GridPos',
           ]


class GridPos:
    """Representation of ARTS GridPos"""
    
    def __init__(self, ind=None, n1=None, n2=None):
        self.ind = ind  # Index
        self.n1 = n1    # Numeric
        self.n2 = n2    # Numeric
        
    @classmethod
    def from_xml(cls, xmlelement):
        """Loads GridPos object from an existing file.
        """
        obj = cls()
        
        obj.ind = int(xmlelement[0].value())
        obj.n1 = xmlelement[1].value()
        obj.n2 = xmlelement[2].value()
        
        return obj
    
    def write_xml(self, xmlwriter, attr=None):
        """Write GridPos object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag("GridPos", attr)
        xmlwriter.write_xml(self.ind, {"name": "OriginalGridIndexBelowInterpolationPoint"})
        xmlwriter.write_xml(self.n1, {"name": "FractionalDistanceToNextPoint_1"})
        xmlwriter.write_xml(self.n2, {"name": "FractionalDistanceToNextPoint_2"})
        xmlwriter.close_tag()
    
    def __repr__(self):
        return "GridPos {}: n1={}; n2={}".format(self.ind, self.n1, self.n2)
    
    def __eq__(self, other):
        return self.ind == other.ind and self.n1 == other.n1 and self.n2 == other.n2
    
    def __ne__(self, other):
        return not (self==other)


class Ppath:
    """Represents the Ppath variable in ARTS"""
    
    def __init__(self):
        self.dim = None          # Index
        self.np = None           # Index
        self.constant = None     # Numeric
        self.background = None   # String
        self.start_pos = None    # Vector
        self.start_los = None    # Vector
        self.start_lstep = None  # Numeric
        self.pos = None          # Matrix
        self.los = None          # Matrix
        self.r = None            # Vector
        self.lstep = None        # Vector
        self.end_pos = None      # Vector
        self.end_los = None      # Vector
        self.end_lstep = None    # Vector
        self.nreal = None        # Vector
        self.ngroup = None       # Vector
        self.gp_p = None         # ArrayOfGridPos
        self.gp_lat = None       # ArrayOfGridPos
        self.gp_lon = None       # ArrayOfGridPos
    
    def __repr__(self):
        return "Ppath of {} steps in {}D Atmosphere".format(self.np, self.dim)
    
    @classmethod
    def from_xml(cls, xmlelement):
        """Loads Ppath object from an existing file.
        """
        obj = cls()
        
        obj.dim = int(xmlelement[0].value())
        obj.np = int(xmlelement[1].value())
        obj.constant = xmlelement[2].value()
        obj.background = xmlelement[3].value()
        obj.start_pos = xmlelement[4].value()
        obj.start_los = xmlelement[5].value()
        obj.start_lstep = xmlelement[6].value()
        obj.pos = xmlelement[7].value()
        obj.los = xmlelement[8].value()
        obj.r = xmlelement[9].value()
        obj.lstep = xmlelement[10].value()
        obj.end_pos = xmlelement[11].value()
        obj.end_los = xmlelement[12].value()
        obj.end_lstep = xmlelement[13].value()
        obj.nreal = xmlelement[14].value()
        obj.ngroup = xmlelement[15].value()
        obj.gp_p = xmlelement[16].value()
        obj.gp_lat = xmlelement[17].value()
        obj.gp_lon = xmlelement[18].value()
        return obj
    
    def write_xml(self, xmlwriter, attr=None):
        """Write Ppath object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag("Ppath", attr)
        xmlwriter.write_xml(self.dim, {"name": "AtmosphericDimensionality"})
        xmlwriter.write_xml(self.np, {"name": "NumberOfPositionInPropagationPath"})
        xmlwriter.write_xml(self.constant, {"name": "PropagationPathConstant"})
        xmlwriter.write_xml(self.background, {"name": "RadiativeBackground"})
        xmlwriter.write_xml(self.start_pos, {"name": "StartPositionOfPropagationPath"})
        xmlwriter.write_xml(self.start_los, {"name": "StartLOSOfPropagationPath"})
        xmlwriter.write_xml(self.start_lstep, {"name": "StartLstepOfPropagationPath"})
        xmlwriter.write_xml(self.pos, {"name": "PropagationPathPointPositions"})
        xmlwriter.write_xml(self.los, {"name": "LineOfSight"})
        xmlwriter.write_xml(self.r, {"name": "PropagationPathPointRadii"})
        xmlwriter.write_xml(self.lstep, {"name": "PropagationPathPositionLength"})
        xmlwriter.write_xml(self.end_pos, {"name": "EndPositionOfPropagationPath"})
        xmlwriter.write_xml(self.end_los, {"name": "EndLOSOfPropagationPath"})
        xmlwriter.write_xml(self.end_lstep, {"name": "EndLstepPropagationPath"})
        xmlwriter.write_xml(self.nreal, {"name": "RefractiveIndexRealPart"})
        xmlwriter.write_xml(self.ngroup, {"name": "GroupRefractiveIndex"})
        xmlwriter.write_xml(self.gp_p, {"name": "PressureGridIndexPosition"}, arraytype="GridPos")
        xmlwriter.write_xml(self.gp_lat, {"name": "LatitudeGridIndexPosition"}, arraytype="GridPos")
        xmlwriter.write_xml(self.gp_lon, {"name": "LongitudeGridIndexPosition"}, arraytype="GridPos")
        xmlwriter.close_tag()
        
    def alt_lat_lon_za_aa(self):
        alt = self.pos[:, 0]
        lat = self.pos[:, 1] if self.pos.shape[1] > 1 else np.zeros_like(alt)
        lon = self.pos[:, 2] if self.pos.shape[1] > 2 else np.zeros_like(alt)
        
        za = self.los[:, 0]
        aa = self.los[:, 1] if self.los.shape[1] > 1 else np.zeros_like(za)
        
        return alt, lat, lon, za, aa


class ArrayOfLineRecord:
    """Represents an :arts:`ArrayOfLineRecord` object."""

    def __init__(self, data=None, version=None):
        self.data = data
        self.version = version

    def __repr__(self):
        if len(self.data) > 1:
            return "ArrayOfLineRecord. " + self.version + ". " + \
                str(len(self.data)) + " lines."
        elif len(self.data) == 1:
            if '@' in self.data[0]:
                return "ArrayOfLineRecord. " + self.version + ". 1 line."
        return "ArrayOfLineRecord. " + self.version + ". No lines."

    def __getitem__(self, index):
        return self.data[index]

    def __len__(self):
        return len(self.data)

    def as_ARTSCAT5(self):
        """Returns manipulable ARTSCAT5 class of this linerecord array
        """
        assert self.version == 'ARTSCAT-5', "Only for ARTSCAT-5 data"
        return ARTSCAT5(self)

    @property
    def version(self):
        """ArrayOfRecord version number."""
        return self._version

    @property
    def data(self):
        """List of strings representing line records."""
        return self._data

    @version.setter
    def version(self, version):
        if version is None:
            self._version = None
            return

        if isinstance(version, str):
            self._version = version
        else:
            raise TypeError('version has to be String.')

    @data.setter
    def data(self, data):
        self._data = return_if_arts_type(data, 'ArrayOfString')

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads an ArrayOfLineRecord object from an existing file.
        """
        obj = cls()
        obj.version = xmlelement.attrib['version']
        obj.data = xmlelement.text.strip().split('\n')
        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write an ArrayOfLineRecord object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        attr['version'] = self.version
        attr['nelem'] = len(self.data)

        xmlwriter.open_tag("ArrayOfLineRecord", attr)
        xmlwriter.write('\n'.join(self.data) + '\n')
        xmlwriter.close_tag()


class CIARecord:
    """Represents a CIARecord object.

    See online ARTS documentation for object details.

    """

    def __init__(self, molecule1=None, molecule2=None, data=None):
        self.molecule1 = molecule1
        self.molecule2 = molecule2
        self.data = data

    @property
    def molecule1(self):
        """Name of the first molecule."""
        return self._molecule1

    @property
    def molecule2(self):
        """Name of the second molecule."""
        return self._molecule2

    @property
    def data(self):
        """Actual data stored in (list of) GriddedField2 objects."""
        return self._data

    @molecule1.setter
    def molecule1(self, molecule1):
        self._molecule1 = return_if_arts_type(molecule1, 'String')

    @molecule2.setter
    def molecule2(self, molecule2):
        self._molecule2 = return_if_arts_type(molecule2, 'String')

    @data.setter
    def data(self, data):
        self._data = return_if_arts_type(data, 'ArrayOfGriddedField2')

    def __repr__(self):
        return self._molecule1 + "-CIA-" + self.molecule2 + " " + \
            str(self.data)

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a CIARecord object from an existing file.
        """

        obj = cls()
        obj.molecule1 = xmlelement.attrib['molecule1']
        obj.molecule2 = xmlelement.attrib['molecule2']
        obj.data = xmlelement[0].value()

        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write a CIARecord object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        attr['molecule1'] = self.molecule1
        attr['molecule2'] = self.molecule2

        xmlwriter.open_tag("CIARecord", attr)
        xmlwriter.write_xml(self.data)
        xmlwriter.close_tag()


# TODO(LKL): consider splitting SpeciesAuxData into seperate classes for each
# version. SpeciesAuxData could be used as wrapper class.
class SpeciesAuxData:
    """Represents a SpeciesAuxData object.

    See online ARTS documentation for object details.

    """

    def __init__(self, data, version, nparam=None):
        self.version = version
        self.nparam = nparam
        self.data = data

    def __repr__(self):
        return "SpeciesAuxData Version " + str(self.version) + ' ' + \
            'for ' + str(len(self.species())) + ' species'

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        self._data = data
        if self.version == 1:
            self._data_dict = {}
            self._keys = {}
            for ii in range(len(data)):
                iso_data = data[ii]
                tmp = iso_data.split()
                self._keys[tmp[1]] = ii
                self._data_dict[tmp[1]] = float(tmp[2])
        elif self.version == 2:
            self._data_dict = {}
            self._keys = {}
            for ii in range(len(data)):
                tmp = data[ii]
                self._keys[tmp[0]] = ii
                self._data_dict[tmp[0]] = [tmp[1], tmp[2]]

    def __getitem__(self, key):
        return self._data_dict[key]

    def __setitem__(self, key, val):
        self._data_dict[key] = val
        if self.version == 1:
            self._data[(self._keys[key])] = '@ ' + key + ' ' + str(val)
        elif self.version == 2:
            self._data[(self._keys[key])] = val

    def __contains__(self, key):
        return key in self._data_dict

    def species(self):
        return list(self._data_dict.keys())

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a SpeciesAuxData object from an existing file.
        """

        version = int(xmlelement.attrib['version'])

        if version == 1:
            nparam = int(xmlelement.attrib['nparam'])
            data = [s for s in xmlelement.text.split('\n') if s != '']
        elif version == 2:
            nparam = None
            data = []
            sub_list = []
            for n, elem in enumerate(xmlelement):
                if n != 0 and n % 3 == 0:
                    data.append(sub_list)
                    sub_list = []
                sub_list.append(elem.value())
            data.append(sub_list)
        else:
            raise Exception(
                "Unknown SpeciesAuxData version {}.".format(version))

        obj = cls(data, version, nparam=nparam)
        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write a ScatterinMetaData object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        attr['version'] = self.version
        attr['nelem'] = len(self.data)

        if self.version == 1:
            attr['nparam'] = self.nparam

            xmlwriter.open_tag("SpeciesAuxData", attr)
            xmlwriter.write('\n'.join(self.data) + '\n')
            xmlwriter.close_tag()
        elif self.version == 2:
            xmlwriter.open_tag("SpeciesAuxData", attr)
            for sub_list in self.data:
                for element in sub_list:
                    xmlwriter.write_xml(element)
            xmlwriter.close_tag()

    def as_PartitionFunctions(self):
        return PartitionFunctions(self)


class GasAbsLookup:
    """Represents a GasAbsLookup object.

    See online ARTS documentation for object details.

    """

    def __init__(self,
                 speciestags=None,
                 nonlinearspecies=None,
                 frequencygrid=None,
                 pressuregrid=None,
                 referencevmrprofiles=None,
                 referencetemperatureprofile=None,
                 temperatureperturbations=None,
                 nonlinearspeciesvmrperturbations=None,
                 absorptioncrosssection=None):

        self.speciestags = speciestags
        self.nonlinearspecies = nonlinearspecies
        self.frequencygrid = frequencygrid
        self.pressuregrid = pressuregrid
        self.referencevmrprofiles = referencevmrprofiles
        self.referencetemperatureprofile = referencetemperatureprofile
        self.temperatureperturbations = temperatureperturbations
        self.nonlinearspeciesvmrperturbations = nonlinearspeciesvmrperturbations
        self.absorptioncrosssection = absorptioncrosssection

    @property
    def speciestags(self):
        """List of :class:`SpeciesTag`."""
        return self._speciestags

    @property
    def nonlinearspecies(self):
        """Indices to indentify nonlinear species."""
        return self._nonlinearspecies

    @property
    def frequencygrid(self):
        """Frequency vector."""
        return self._frequencygrid

    @property
    def pressuregrid(self):
        """Pressure level vector."""
        return self._pressuregrid

    @property
    def referencevmrprofiles(self):
        """Reference VMR profiles."""
        return self._referencevmrprofiles

    @property
    def referencetemperatureprofile(self):
        """Reference temperature profile."""
        return self._referencetemperatureprofile

    @property
    def temperatureperturbations(self):
        """Vector with temperature perturbations."""
        return self._temperatureperturbations

    @property
    def nonlinearspeciesvmrperturbations(self):
        """Vector with VMR perturbations for nonlinear species."""
        return self._nonlinearspeciesvmrperturbations

    @property
    def absorptioncrosssection(self):
        """Absorption crosssections."""
        return self._absorptioncrosssection

    @speciestags.setter
    def speciestags(self, speciestags):
        self._speciestags = return_if_arts_type(
            speciestags, 'ArrayOfArrayOfSpeciesTag')

    @nonlinearspecies.setter
    def nonlinearspecies(self, nonlinearspecies):
        self._nonlinearspecies = return_if_arts_type(
            nonlinearspecies, 'ArrayOfIndex')

    @frequencygrid.setter
    def frequencygrid(self, frequencygrid):
        self._frequencygrid = return_if_arts_type(
            frequencygrid, 'Vector')

    @pressuregrid.setter
    def pressuregrid(self, pressuregrid):
        self._pressuregrid = return_if_arts_type(
            pressuregrid, 'Vector')

    @referencevmrprofiles.setter
    def referencevmrprofiles(self, referencevmrprofiles):
        self._referencevmrprofiles = return_if_arts_type(
            referencevmrprofiles, 'Matrix')

    @referencetemperatureprofile.setter
    def referencetemperatureprofile(self, referencetemperatureprofile):
        self._referencetemperatureprofile = return_if_arts_type(
            referencetemperatureprofile, 'Vector')

    @temperatureperturbations.setter
    def temperatureperturbations(self, temperatureperturbations):
        self._temperatureperturbations = return_if_arts_type(
            temperatureperturbations, 'Vector')

    @nonlinearspeciesvmrperturbations.setter
    def nonlinearspeciesvmrperturbations(self, nonlinearspeciesvmrperturbations):
        self._nonlinearspeciesvmrperturbations = return_if_arts_type(
            nonlinearspeciesvmrperturbations, 'Vector')

    @absorptioncrosssection.setter
    def absorptioncrosssection(self, absorptioncrosssection):
        self._absorptioncrosssection = return_if_arts_type(
            absorptioncrosssection, 'Tensor4')

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a GasAbsLookup object from an existing file.
        """

        obj = cls()
        obj.speciestags = xmlelement[0].value()
        obj.nonlinearspecies = xmlelement[1].value()
        obj.frequencygrid = xmlelement[2].value()
        obj.pressuregrid = xmlelement[3].value()
        obj.referencevmrprofiles = xmlelement[4].value()
        obj.referencetemperatureprofile = xmlelement[5].value()
        obj.temperatureperturbations = xmlelement[6].value()
        obj.nonlinearspeciesvmrperturbations = xmlelement[7].value()
        obj.absorptioncrosssection = xmlelement[8].value()

        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write a ScatterinMetaData object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag("GasAbsLookup", attr)
        xmlwriter.write_xml(self.speciestags)
        if self.nonlinearspecies is None:
            nonlinearspecies = []
        else:
            nonlinearspecies = self.nonlinearspecies
        xmlwriter.write_xml(nonlinearspecies,
                            {'name': 'NonlinearSpecies'},
                            arraytype='Index')
        xmlwriter.write_xml(self.frequencygrid,
                            {'name': 'FrequencyGrid'})
        xmlwriter.write_xml(self.pressuregrid,
                            {'name': 'PressureGrid'})
        xmlwriter.write_xml(self.referencevmrprofiles,
                            {'name': 'ReferenceVmrProfiles'})
        xmlwriter.write_xml(self.referencetemperatureprofile,
                            {'name': 'ReferenceTemperatureProfile'})
        xmlwriter.write_xml(self.temperatureperturbations,
                            {'name': 'TemperaturePerturbations'})
        xmlwriter.write_xml(self.nonlinearspeciesvmrperturbations,
                            {'name': 'NonlinearSpeciesVmrPerturbations'})
        xmlwriter.write_xml(self.absorptioncrosssection,
                            {'name': 'AbsorptionsCrossSections'})
        xmlwriter.close_tag()


class SpeciesTag(str):
    """Represents a SpeciesTag object.

    See online ARTS documentation for object details.

    """

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a SpeciesTag object from an existing file.
        """
        if xmlelement.text is None:
            raise Exception('SpeciesTag must not be empty.')
        return cls(xmlelement.text.strip()[1:-1])

    def write_xml(self, xmlwriter, attr=None):
        """Write a SpeciesTag object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag('SpeciesTag', attr, newline=False)
        xmlwriter.write('"' + self + '"')
        xmlwriter.close_tag()


class Sparse(scipy.sparse.csc_matrix):
    """Wrapper around :class:`scipy.sparse.csc_matrix`.

    This class wraps around the SciPy Compressed Sparse Column matrix. The
    usage is exactly the same, but support for reading and writing XML files
    is added. Also additional attributes were added to map the ARTS
    implementation of :arts:`Sparse`.

    """
    @property
    def nrows(self):
        """Number of rows."""
        return self.shape[0]

    @property
    def ncols(self):
        """Number of columns."""
        return self.shape[0]

    @property
    def rowindex(self):
        """Row indices to locate data in matrix."""
        return self.tocoo().row

    @property
    def colindex(self):
        """Column indices to locate data in matrix."""
        return self.tocoo().col

    @property
    def sparsedata(self):
        """Data value at specified positions in matrix."""
        return self.tocoo().data

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a Sparse object from an existing file."""

        binaryfp = xmlelement.binaryfp
        nelem = int(xmlelement[0].attrib['nelem'])
        nrows = int(xmlelement.attrib['nrows'])
        ncols = int(xmlelement.attrib['ncols'])

        if binaryfp is None:
            rowindex = np.fromstring(xmlelement[0].text, sep=' ').astype(int)
            colindex = np.fromstring(xmlelement[1].text, sep=' ').astype(int)
            sparsedata = np.fromstring(xmlelement[2].text, sep=' ')
        else:
            rowindex = np.fromfile(binaryfp, dtype='<i4', count=nelem)
            colindex = np.fromfile(binaryfp, dtype='<i4', count=nelem)
            sparsedata = np.fromfile(binaryfp, dtype='<d', count=nelem)

        return cls((sparsedata, (rowindex, colindex)), [nrows, ncols])

    def write_xml(self, xmlwriter, attr=None):
        """Write a Sparse object to an ARTS XML file."""

        # Get ARTS-style information from CSC matrix.
        nrows = self.shape[0]
        ncols = self.shape[1]
        rowindex = self.tocoo().row
        colindex = self.tocoo().col
        sparsedata = self.tocoo().data

        precision = xmlwriter.precision

        if attr is None:
            attr = {}

        attr['nrows'] = nrows
        attr['ncols'] = ncols

        xmlwriter.open_tag('Sparse', attr)

        binaryfp = xmlwriter.binaryfilepointer

        if binaryfp is None:
            xmlwriter.open_tag('RowIndex', {'nelem': rowindex.size})
            for i in rowindex:
                xmlwriter.write('%d' % i + '\n')
            xmlwriter.close_tag()
            xmlwriter.open_tag('ColIndex', {'nelem': colindex.size})
            for i in colindex:
                xmlwriter.write('%d' % i + '\n')
            xmlwriter.close_tag()
            xmlwriter.open_tag('SparseData', {'nelem': sparsedata.size})
            for i in sparsedata:
                xmlwriter.write(('%' + precision) % i + '\n')
            xmlwriter.close_tag()
            xmlwriter.close_tag()
        else:
            xmlwriter.open_tag('RowIndex', {'nelem': rowindex.size})
            np.array(rowindex, dtype='i4').tofile(binaryfp)
            xmlwriter.close_tag()
            xmlwriter.open_tag('ColIndex', {'nelem': colindex.size})
            np.array(colindex, dtype='i4').tofile(binaryfp)
            xmlwriter.close_tag()
            xmlwriter.open_tag('SparseData', {'nelem': sparsedata.size})
            np.array(sparsedata, dtype='d').tofile(binaryfp)
            xmlwriter.close_tag()
            xmlwriter.close_tag()


class QuantumIdentifier:
    """Represents a QuantumIdentifier object.

    See online ARTS documentation for object details.

    """

    def __init__(self, qid):

        assert type(qid) is str, "Need String input"
        these = qid.split()
        assert len(these) > 0, "No QuantumIdentifier"
        spec = these[0].split('-')  # UGLY: What about negative charge?
        if len(spec) == 1:
            self._afgl = None
            if spec[0] == 'None':
                self._spec = None
            else:
                self._spec = spec[0]
        elif len(spec) == 2:
            self._spec = spec[0]
            self._afgl = int(spec[1])
        else:
            assert False, "Cannot recognize species"

        if len(these) == 1:
            self._transition = False
            self._level = False
            return

        if these[1] == 'TR':
            self._transition = True
            self._level = False
        elif these[1] == 'EN':
            self._transition = False
            self._level = True
        else:
            assert False, "Must be energy level [EN] or transition [TR] type"

        self._qns = as_quantumnumbers(" ".join(these[2:]))

        self._assert_sanity_()

    def __repr__(self):
        out = str(self._spec)
        if self._afgl is not None:
            out += '-' + str(self._afgl)
        if self._transition or self._level:
            if self._transition:
                out += ' TR '
            else:
                out += ' EN '
            out += str(self._qns)
        return out

    def _assert_sanity_(self):
        if self._transition:
            assert type(self._qns) is QuantumNumberRecord, "Mismatching types"
        elif self._level:
            assert type(self._qns) is QuantumNumbers, "Mismatching types"
        else:
            assert False, "Programmer error?"

    def __str__(self):
        assert self.afgl is not None or self.species is not None, \
            "Bad data cannot be converted to str.  Contains no species or iso"
        return self.__repr__()

    @property
    def qns(self):
        return self._qns

    @qns.setter
    def qns(self, qns):
        self._qns = as_quantumnumbers(qns)
        if type(self._qns) is QuantumNumberRecord:
            self._transition = True
            self._level = False
        elif type(self._qns) is QuantumNumbers:
            self._transition = False
            self._level = True
        else:
            assert False, "Programmer error?"

    @property
    def species(self):
        return self._spec

    @species.setter
    def species(self, value):
        self._spec = return_if_arts_type(value, 'String')

    @property
    def afgl(self):
        return self._afgl

    @afgl.setter
    def afgl(self, value):
        self._afgl = return_if_arts_type(value, 'Index')

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a QuantumIdentifier object from an existing file.
        """
        if xmlelement.text is None:
            raise Exception('QuantumIdentifier must not be empty.')
        return cls(xmlelement.text.strip())

    def write_xml(self, xmlwriter, attr=None):
        """Write a QuantumIdentifier object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag('QuantumIdentifier', attr, newline=False)
        xmlwriter.write(self.__str__())
        xmlwriter.close_tag()


class QuantumNumberRecord:
    """Represents a QuantumNumberRecord object.

    See online ARTS documentation for object details.

    """

    def __init__(self, upper=None, lower=None):
        self._qns = {'UP': QuantumNumbers(), 'LO': QuantumNumbers()}
        self._qns['UP'] = return_if_arts_type(upper, 'QuantumNumbers')
        self._qns['LO'] = return_if_arts_type(lower, 'QuantumNumbers')

    def __repr__(self):
        if len(self._qns['UP']) == 0 and len(self._qns['LO']) == 0:
            return 'No Quantum-Numbers'
        else:
            return "UP " + str(self._qns['UP']) + " LO " + str(self._qns['LO'])

    def __str__(self):
        if len(self._qns['UP']) == 0 and len(self._qns['LO']) == 0:
            return ''
        else:
            return self.__repr__()

    def __getitem__(self, key):
        return self._qns[key]

    def __setitem__(self, key, value):
        self._qns[key] = return_if_arts_type(as_quantumnumbers(value),
                                             'QuantumNumbers')

    def __iter__(self):
        return iter(self._qns)

    def __contains__(self, value):
        return value in ['UP', 'LO']

    @classmethod
    def from_dict(cls, d):
        """Creates a QuantumNumberRecord from dictionary
        """
        if len(d) == 0:
            return QuantumNumberRecord(upper=QuantumNumbers(),
                                       lower=QuantumNumbers())

        assert 'UP' in d and 'LO' in d, "Need UP and LO to create"
        qnr = QuantumNumberRecord(upper=QuantumNumbers(d['UP']),
                                  lower=QuantumNumbers(d['LO']))
        return qnr

    @classmethod
    def from_str(cls, s):
        """Creates a QuantumNumberRecord from string
        """
        s = s.strip()
        if len(s) == 0:
            return QuantumNumberRecord(upper=QuantumNumbers(),
                                       lower=QuantumNumbers())

        assert 'UP' in s and 'LO' in s, "Need UP and LO to create"
        _t1 = s.split('UP')
        assert len(_t1) == 2, "Unexpectedly many/few UP in s"
        if len(_t1[0]) == 0:
            _t2 = _t1[1].split('LO')
            assert len(_t2) == 2, "Unexpectedly many/few LO in s"
            lo = _t2[1]
            up = _t2[0]
        else:
            up = _t1[1]
            _t2 = _t1[0].split('LO')
            assert len(_t2) == 2, "Unexpectedly many/few LO in s"
            lo = _t2[1]

        qnr = QuantumNumberRecord(upper=QuantumNumbers(up),
                                  lower=QuantumNumbers(lo))
        return qnr

    @property
    def upper(self):
        """QuantumNumbers object representing the upper quantumnumber."""
        return self._qns['UP']

    @property
    def lower(self):
        """QuantumNumbers object representing the lower quantumnumber."""
        return self._qns['LO']

    @upper.setter
    def upper(self, upper):
        self._qns['UP'] = return_if_arts_type(upper, 'QuantumNumbers')

    @lower.setter
    def lower(self, lower):
        self._qns['LO'] = return_if_arts_type(lower, 'QuantumNumbers')

    @property
    def qns(self):
        return self._qns

    @qns.setter
    def qns(self, value):
        if 'LO' in value:
            self._qns['LO'] = QuantumNumbers(value['LO'])
        else:
            self._qns['LO'] = QuantumNumbers()

        if 'UP' in value:
            self._qns['UP'] = QuantumNumbers(value['UP'])
        else:
            self._qns['UP'] = QuantumNumbers()

    def zeeman_splitting(self, type=None, case=None, H=1):
        from ..physics.em import zeeman_splitting, landau_g_factor, \
            zeeman_transitions

        if case.lower() == 'a':
            gu = landau_g_factor(self.upper['Omega'], self.upper['J'],
                                 self.upper['S'], self.upper['Lambda'],
                                 gs=2, gl=1, case='a')
            gl = landau_g_factor(self.lower['Omega'], self.lower['J'],
                                 self.lower['S'], self.lower['Lambda'],
                                 gs=2, gl=1, case='a')
        elif case.lower() == 'b':
            gu = landau_g_factor(self.upper['N'], self.upper['J'],
                                 self.upper['S'], self.upper['Lambda'],
                                 gs=2, gl=1, case='b')
            gl = landau_g_factor(self.lower['N'], self.lower['J'],
                                 self.lower['S'], self.lower['Lambda'],
                                 gs=2, gl=1, case='b')
        else:
            raise RuntimeError("No unknown cases allowed")

        if type is not None:
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], type)
            return zeeman_splitting(gu, gl, mu, ml, H)
        else:
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], "PI")
            pi = zeeman_splitting(gu, gl, mu, ml, H)
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], "S+")
            sp = zeeman_splitting(gu, gl, mu, ml, H)
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], "S-")
            sm = zeeman_splitting(gu, gl, mu, ml, H)
            return {"pi": pi, "s+": sp, "s-": sm}

    def zeeman_transitions(self, type=None):
        from ..physics.em import zeeman_transitions, zeeman_strength
        if type is not None:
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], type)
            return zeeman_strength(self.upper['J'], self.lower['J'], mu, ml)
        else:
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], "PI")
            pi = zeeman_strength(self.upper['J'], self.lower['J'], mu, ml)
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], "S+")
            sp = zeeman_strength(self.upper['J'], self.lower['J'], mu, ml)
            mu, ml = zeeman_transitions(self.upper['J'], self.lower['J'], "S-")
            sm = zeeman_strength(self.upper['J'], self.lower['J'], mu, ml)
            return {"pi": pi, "s+": sp, "s-": sm}

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a QuantumNumberRecord object from an existing file.
        """

        obj = cls()
        obj.upper = xmlelement[0][0].value()
        obj.lower = xmlelement[1][0].value()

        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write a SpeciesTag object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag('QuantumNumberRecord', attr)
        xmlwriter.open_tag('Upper', attr, newline=False)
        xmlwriter.write_xml(self.upper)
        xmlwriter.close_tag()
        xmlwriter.open_tag('Lower', attr, newline=False)
        xmlwriter.write_xml(self.lower)
        xmlwriter.close_tag()
        xmlwriter.close_tag()

    def __iadd__(self, qnr):
        self._qns['UP'] += qnr['UP']
        self._qns['LO'] += qnr['LO']

    def __isub__(self, qnr):
        self._qns['UP'] -= qnr['UP']
        self._qns['LO'] -= qnr['LO']

    def __eq__(self, qns):
        if type(qns) is QuantumNumberRecord:
            return self['LO'] == qns['LO'] and self['UP'] == qns['UP']
        elif type(qns) is QuantumNumbers:
            return self['UP'] == qns, self['LO'] == qns
        else:
            return self == as_quantumnumbers(qns)

    def __ne__(self, qns):
        if type(qns) is QuantumNumberRecord:
            return self['LO'] != qns['LO'] and self['UP'] != qns['UP']
        elif type(qns) is QuantumNumbers:
            return self['UP'] != qns, self['LO'] != qns
        else:
            return self != as_quantumnumbers(qns)

    def __lt__(self, qns):
        if type(qns) is QuantumNumberRecord:
            return self['LO'] < qns['LO'] and self['UP'] < qns['UP']
        elif type(qns) is QuantumNumbers:
            return self['UP'] < qns, self['LO'] < qns
        else:
            return self < as_quantumnumbers(qns)

    def __gt__(self, qns):
        if type(qns) is QuantumNumberRecord:
            return self['LO'] > qns['LO'] and self['UP'] > qns['UP']
        elif type(qns) is QuantumNumbers:
            return self['UP'] > qns, self['LO'] > qns
        else:
            return self > as_quantumnumbers(qns)

    def __le__(self, qns):
        if type(qns) is QuantumNumberRecord:
            return self['LO'] <= qns['LO'] and self['UP'] <= qns['UP']
        elif type(qns) is QuantumNumbers:
            return self['UP'] <= qns, self['LO'] <= qns
        else:
            return self <= as_quantumnumbers(qns)

    def __ge__(self, qns):
        if type(qns) is QuantumNumberRecord:
            return self['LO'] >= qns['LO'] and self['UP'] >= qns['UP']
        elif type(qns) is QuantumNumbers:
            return self['UP'] >= qns, self['LO'] >= qns
        else:
            return self >= as_quantumnumbers(qns)


class QuantumNumbers:
    """Represents a QuantumNumbers object.

    See online ARTS documentation for object details.

    """

    def __init__(self, numbers=None, nelem=None):

        self.numbers = numbers
        if nelem is not None:
            self.nelem = nelem
        else:
            self.nelem = len(self.numbers)

        self._assert_sanity_()

    def _assert_sanity_(self):
        if self.nelem is None or self.numbers is None:
            return
        assert len(self.numbers) == self.nelem, "mismatching quantum numbers"

    def __repr__(self):
        out = ''
        for qn in self.numbers:
            out += qn + ' ' + str(self.numbers[qn]) + ' '
        return out[:-1]

    def __getitem__(self, key):
        """Returns the value.  Mimics ARTS' behavior for mismatched data
        """
        if key in self:
            return self.numbers[key]
        else:
            return None

    def __setitem__(self, key, value):
        """Sets a value and counts up the quantum numbers
        """
        if key in self.numbers:
            self.numbers[key] = Rational(value)
        else:
            self.numbers[key] = Rational(value)
            self.nelem += 1
        self._assert_sanity_()

    def __iadd__(self, qns):
        for qn in qns:
            assert qn not in self, "Addition means adding new QN. Access " + \
                "individual elements to change their values"
            self.numbers[qn] = qns[qn]
            self.nelem += 1
        return self

    def __isub__(self, qns):
        for qn in qns:
            assert qn in self, "Subtraction means removing QN. Access " + \
                "individual elements to change their values"
            del self.numbers[qn]
            self.nelem -= 1
        return self

    def __contains__(self, key):
        """Are these quantum numbers here?
        """
        return key in self.numbers

    def __iter__(self):
        return iter(self.numbers)

    def __eq__(self, qns):
        """Tests for complete equality ==
        """
        return self <= qns and len(qns) == self.nelem

    def __ne__(self, qns):
        """Tests for lacking complete equality !=
        """
        return not self == qns

    def __le__(self, qns):
        """Tests for all in self being in qns <=
        """
        try:
            for qn in self:
                if qns[qn] != self[qn]:
                    return False
            return True
        except:
            return False

    def __ge__(self, qns):
        """Tests for all in qns being in self >=
        """
        try:
            for qn in qns:
                if qns[qn] != self[qn]:
                    return False
            return True
        except:
            return False

    def __lt__(self, qns):
        """Tests for all in self being in qns and if there is more in qns <
        """
        return self <= qns and self.nelem < len(qns)

    def __gt__(self, qns):
        """Tests for all in self being in qns and if there is more in self >
        """
        return qns <= self and len(qns) < self.nelem

    def __len__(self):
        return self.nelem

    def array_of_M(self):
        """Returns all possible M in a list.  Requires presence of J
        """
        assert 'J' in self, "Need J to define M"
        assert self['J'] >= 0, "Negative J in this instance?"
        _t = []
        _s = -self['J']
        while _s <= self['J']:
            _t.append(_s)
            _s += 1
        return np.array(_t)

    @property
    def numbers(self):
        """Dict representing the quantumnumbers."""
        return self._numbers

    @property
    def nelem(self):
        """Number of quantumnumbers stored."""
        return self._nelem

    @numbers.setter
    def numbers(self, numbers):
        if type(numbers) is str:
            _t = numbers.split()
            nums = {}
            i = 0
            assert len(_t) % 2 == 0, "Not of form 'key1 value1 key2 value2'"
            while i < len(_t):
                nums[_t[i]] = Rational(_t[i+1])
                i += 2
            self._numbers = nums
        elif type(numbers) is dict:
            for i in numbers:
                numbers[i] = Rational(numbers[i])
            self._numbers = numbers
        elif type(numbers) is QuantumNumbers:
            self._numbers = numbers.numbers
        elif numbers is None:
            self._numbers = {}
        else:
            assert False, "Expected dict or String for QuantumNumbers"
        # OLD: self._numbers = return_if_arts_type(numbers, 'String')

    @nelem.setter
    def nelem(self, nelem):
        if nelem is None:
            self._nelem = None
            return

        self._nelem = nelem

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a QuantumNumbers object from an existing file.
        """

        obj = cls()
        obj.numbers = xmlelement.text
        obj.nelem = int(xmlelement.attrib['nelem'])

        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write a SpeciesTag object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        attr['nelem'] = self.nelem

        xmlwriter.open_tag('QuantumNumbers', attr, newline=False)
        xmlwriter.write(self.__str__())
        xmlwriter.close_tag(newline=False)


class LineMixingRecord:
    """Represents a LineMixingRecord object.

    See online ARTS documentation for object details.

    """

    def __init__(self, tag=None, quantumnumberrecord=None, data=None):

        self.tag = tag
        self.quantumnumberrecord = quantumnumberrecord
        self.data = LineMixing(data)

    @property
    def tag(self):
        """:class:`SpeciesTag`"""
        return self._tag

    @property
    def quantumnumberrecord(self):
        """:class:`QuantumNumberRecord`"""
        return self._quantumnumberrecord

    @property
    def data(self):
        """Lineshape parameters."""
        return self._data

    @tag.setter
    def tag(self, tag):
        if tag is None:
            self._tag = None
            return

        self._tag = SpeciesTag(tag)

    def __repr__(self):
        return self.tag + ' ' + str(self.quantumnumberrecord) + ' ' + \
            str(self.data)

    @quantumnumberrecord.setter
    def quantumnumberrecord(self, quantumnumberrecord):
        self._quantumnumberrecord = return_if_arts_type(
            quantumnumberrecord, 'QuantumNumberRecord')

    @data.setter
    def data(self, data):
        self._data = return_if_arts_type(data, 'LineMixing')

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a LineMixingRecord object from an existing file.
        """

        obj = cls()
        obj.tag = xmlelement[0].value()
        obj.quantumnumberrecord = xmlelement[1].value()
        obj.data = LineMixing(xmlelement[2].value())

        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write a LineMixingRecord object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag("LineMixingRecord", attr)
        xmlwriter.write_xml(self.tag)
        xmlwriter.write_xml(self.quantumnumberrecord)
        xmlwriter.write_xml(self.data.data)
        xmlwriter.close_tag()


class PropagationMatrix:
    """Represents a PropagationMatrix object.

    See online ARTS documentation for object details.

    """

    def __init__(self, data=np.zeros((1))):
        n = len(data.shape)
        if n == 4:
            self.aa = data.shape[0]
            self.za = data.shape[1]
            self.nf = data.shape[2]
        elif n == 3:
            self.aa = 1
            self.za = data.shape[0]
            self.nf = data.shape[1]
        elif n == 2:
            self.aa = 1
            self.za = 1
            self.nf = data.shape[0]
        elif n != 1:
            raise RuntimeError("Bad input")
        else:
            self.aa = 1
            self.za = 1
            self.nf = 1

        if data.shape[-1] == 7:
            self.stokes = 4
        elif data.shape[-1] == 4:
            self.stokes = 3
        elif data.shape[-1] == 1 or data.shape[-1] == 2:
            self.stokes = data.shape[-1]
        else:
            raise RuntimeError("Bad input")

        self.data = data.reshape(self.aa, self.za, self.nf, data.shape[-1])
    
    def Kjj(self, aa=0, za=0):
        return self.data[aa, za, :, 0]
    
    def K12(self, aa=0, za=0):
        return self.data[aa, za, :, 1]
    
    def K13(self, aa=0, za=0):
        return self.data[aa, za, :, 2]
    
    def K14(self, aa=0, za=0):
        return self.data[aa, za, :, 3]
    
    def K23(self, aa=0, za=0):
        return self.data[aa, za, :, 4]
    
    def K24(self, aa=0, za=0):
        return self.data[aa, za, :, 5]
    
    def K34(self, aa=0, za=0):
        return self.data[aa, za, :, 6]
    
    def __add__(self, other):
        if isinstance(other, PropagationMatrix):
            return PropagationMatrix(self.data + other.data)
        else:
            return PropagationMatrix(self.data + other)
    __radd__ = __add__
    
    def __repr__(self):
        size = "Stokes Dim: {}, Freqs: {}, Zeniths: {}, Azimuths: {}".format(self.stokes, self.nf, self.za, self.aa)
        return "PropagationMatrix of size <{}>".format(size)
    __str__=__repr__

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a PropagationMatrix object from an existing file.
        """
        return cls(xmlelement[0].value())

    def write_xml(self, xmlwriter, attr=None):
        """Write a PropagationMatrix object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag("PropagationMatrix", attr)
        xmlwriter.write_xml(self.data)
        xmlwriter.close_tag()


class StokesVector:
    """Represents a StokesVector object.

    See online ARTS documentation for object details.

    """

    def __init__(self, data=np.zeros((1))):
        n = len(data.shape)
        if n == 4:
            self.aa = data.shape[0]
            self.za = data.shape[1]
            self.nf = data.shape[2]
        elif n == 3:
            self.aa = 1
            self.za = data.shape[0]
            self.nf = data.shape[1]
        elif n == 2:
            self.aa = 1
            self.za = 1
            self.nf = data.shape[0]
        elif n != 1:
            raise RuntimeError("Bad input")
        else:
            self.aa = 1
            self.za = 1
            self.nf = 1

        if data.shape[-1] == 1 or data.shape[-1] == 2 or \
           data.shape[-1] == 3 or data.shape[-1] == 4:
            self.stokes = data.shape[-1]
        else:
            raise RuntimeError("Bad input")

        self.data = data.reshape(self.aa, self.za, self.nf, data.shape[-1])

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a StokesVector object from an existing file.
        """
        return cls(xmlelement[0].value())

    def write_xml(self, xmlwriter, attr=None):
        """Write a StokesVector object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag("StokesVector", attr)
        xmlwriter.write_xml(self.data)
        xmlwriter.close_tag()


    from pyarts.utils.arts import return_if_arts_type, as_quantumnumbers
    from pyarts.internals import PartitionFunctions, ARTSCAT5, Rational, \
        LineMixing
