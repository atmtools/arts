# -*- coding: utf-8 -*-

"""Read ARTS XML types

This packages contains the internal implementation for reading ARTS XML files.
"""

from xml.etree import ElementTree

import numpy as np

from .names import dimension_names
from .. import types

__all__ = ['parse']


class ARTSTypesLoadMultiplexer:
    """Used by the xml.etree.ElementTree to parse ARTS variables.

    Tag names in the XML file are mapped to the corresponding parsing method.

    """

    @staticmethod
    def arts(elem):
        if (elem.attrib['format'] not in ('ascii', 'binary')):
            raise RuntimeError('Unknown format in <arts> tag: {}'.format(
                elem.attrib['format']))

        ret = elem[0].value()

        # Try next element, if return value is None (comment tags).
        n = 1
        while ret is None:
            try:
                ret = elem[n].value()
                n += 1
            except:
                break

        return ret

    @staticmethod
    def comment(elem):
        return

    @staticmethod
    def Array(elem):
        arr = [t.value() for t in elem]
        if len(arr) != int(elem.attrib['nelem']):
            raise RuntimeError('Expected {:s} elements in Array, found {:d}'
                               ' elements!'.format(elem.attrib['nelem'],
                                                   len(arr)))
        return arr

    @staticmethod
    def String(elem):
        if elem.text is None:
            return ''
        return elem.text.strip()[1:-1]

    @staticmethod
    def Index(elem):
        if elem.binaryfp is not None:
            return np.fromfile(elem.binaryfp, dtype='<i4', count=1)[0]
        else:
            return int(elem.text)

    @staticmethod
    def Numeric(elem):
        if elem.binaryfp is not None:
            return np.fromfile(elem.binaryfp, dtype='<d', count=1)[0]
        else:
            return float(elem.text)

    @staticmethod
    def Vector(elem):
        nelem = int(elem.attrib['nelem'])
        if nelem == 0:
            arr = np.ndarray((0,))
        else:
            # sep=' ' seems to work even when separated by newlines, see
            # http://stackoverflow.com/q/31882167/974555
            if elem.binaryfp is not None:
                arr = np.fromfile(elem.binaryfp, dtype='<d', count=nelem)
            else:
                arr = np.fromstring(elem.text, sep=' ')
            if arr.size != nelem:
                raise RuntimeError(
                    'Expected {:s} elements in Vector, found {:d}'
                    ' elements!'.format(elem.attrib['nelem'],
                                        arr.size))
        return arr

    @staticmethod
    def ComplexVector(elem):
        nelem = int(elem.attrib['nelem'])
        if nelem == 0:
            arr = np.ndarray((0,), dtype=np.complex128)
        else:
            # sep=' ' seems to work even when separated by newlines, see
            # http://stackoverflow.com/q/31882167/974555
            if elem.binaryfp is not None:
                arr = np.fromfile(elem.binaryfp, dtype=np.complex128,
                                  count=nelem)
            else:
                arr = np.fromstring(elem.text, sep=' ', dtype=np.float64)
                arr.dtype = np.complex128
            if arr.size != nelem:
                raise RuntimeError(
                    'Expected {:s} elements in Vector, found {:d}'
                    ' elements!'.format(elem.attrib['nelem'],
                                        arr.size))
        return arr

    @staticmethod
    def Matrix(elem):
        # turn dims around: in ARTS, [10 x 1 x 1] means 10 pages, 1 row, 1 col
        dimnames = [dim for dim in dimension_names
                    if dim in elem.attrib.keys()][::-1]
        dims = [int(elem.attrib[dim]) for dim in dimnames]
        if np.prod(dims) == 0:
            flatarr = np.ndarray(dims)
        elif elem.binaryfp is not None:
            flatarr = np.fromfile(elem.binaryfp, dtype=np.float64,
                                  count=np.prod(np.array(dims)).item())
            flatarr = flatarr.reshape(dims)
        else:
            flatarr = np.fromstring(elem.text, sep=' ')
            flatarr = flatarr.reshape(dims)
        return flatarr

    @staticmethod
    def ComplexMatrix(elem):
        # turn dims around: in ARTS, [10 x 1 x 1] means 10 pages, 1 row, 1 col
        dimnames = [dim for dim in dimension_names
                    if dim in elem.attrib.keys()][::-1]
        dims = [int(elem.attrib[dim]) for dim in dimnames]
        if np.prod(dims) == 0:
            flatarr = np.ndarray(dims, dtype=np.complex128)
        elif elem.binaryfp is not None:
            flatarr = np.fromfile(elem.binaryfp, dtype=np.complex128,
                                  count=np.prod(np.array(dims)).item())
            flatarr = flatarr.reshape(dims)
        else:
            flatarr = np.fromstring(elem.text, sep=' ', dtype=np.float64)
            flatarr.dtype = np.complex128
            flatarr = flatarr.reshape(dims)
        return flatarr

    Tensor3 = Tensor4 = Tensor5 = Tensor6 = Tensor7 = Matrix
    ComplexTensor3 = ComplexTensor4 = ComplexTensor5 = ComplexTensor6 = ComplexTensor7 = ComplexMatrix


class ARTSElement(ElementTree.Element):
    """Element with value interpretation."""
    binaryfp = None

    def value(self):
        if hasattr(types, self.tag):
            try:
                return types.classes[self.tag].from_xml(self)
            except AttributeError:
                raise RuntimeError('Type {} exists, but has no XML parsing '
                                   'support.'.format(self.tag))
        else:
            try:
                return getattr(ARTSTypesLoadMultiplexer, self.tag)(self)
            except AttributeError:
                raise RuntimeError('Unknown ARTS type {}'.format(self.tag))


def parse(source, binaryfp=None):
    """Parse ArtsXML file from source.

    Args:
        source (str): Filename or file pointer.

    Returns:
        xml.etree.ElementTree: XML Tree of the ARTS data file.

    """
    arts_element = type('ARTSElementBinaryFP',
                        ARTSElement.__bases__,
                        dict(ARTSElement.__dict__))
    arts_element.binaryfp = binaryfp
    return ElementTree.parse(source,
                             parser=ElementTree.XMLParser(
                                 target=ElementTree.TreeBuilder(
                                     element_factory=arts_element)))
