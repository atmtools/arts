# -*- coding: utf-8 -*-
"""
Implementation of RetrievalQuantity.

"""

from arts.utils.arts import return_if_arts_type

__all__ = ['RetrievalQuantity',
           ]

class RetrievalQuantity:
    """Represents a RetrievalQuantity object.

    See online ARTS documentation for object details.

    """

    def __init__(self, maintag=None, subtag=None, subsubtag=None, mode=None,
                 analytical=None, perturbation=None, grids=None):

        self.maintag = maintag
        self.subtag = subtag
        self.subsubtag = subsubtag
        self.mode = mode
        self.analytical = analytical
        self.perturbation = perturbation
        self.grids = grids

    @property
    def maintag(self):
        """MainTag of retrieval species."""
        return self._maintag

    @property
    def subtag(self):
        """Subtag of retrieval species."""
        return self._subtag

    @property
    def subsubtag(self):
        """Subsubtag of retrieval species."""
        return self._subsubtag

    @property
    def mode(self):
        """Retrieval mode."""
        return self._mode

    @property
    def analytical(self):
        """Flag to determine whether the retrieval was done analytically."""
        return self._analytical

    @property
    def perturbation(self):
        """Amplitude of the perturbation."""
        return self._perturbation

    @property
    def grids(self):
        """Pressure grid."""
        return self._grids

    @maintag.setter
    def maintag(self, maintag):
        self._maintag = return_if_arts_type(maintag, 'String')

    @subtag.setter
    def subtag(self, subtag):
        self._subtag = return_if_arts_type(subtag, 'String')

    @subsubtag.setter
    def subsubtag(self, subsubtag):
        self._subsubtag = return_if_arts_type(subsubtag, 'String')

    @mode.setter
    def mode(self, mode):
        self._mode = return_if_arts_type(mode, 'String')

    @analytical.setter
    def analytical(self, analytical):
        self._analytical = return_if_arts_type(analytical, 'Index')

    @perturbation.setter
    def perturbation(self, perturbation):
        self._perturbation = return_if_arts_type(perturbation, 'Numeric')

    @grids.setter
    def grids(self, grids):
        self._grids = return_if_arts_type(grids, 'ArrayOfVector')

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a RetrievalQuantity object from an existing file.

        """
        obj = cls()
        obj.maintag = xmlelement[0].value()
        obj.subtag = xmlelement[1].value()
        obj.subsubtag = xmlelement[2].value()
        obj.mode = xmlelement[3].value()
        obj.analytical = xmlelement[4].value()
        obj.perturbation = xmlelement[5].value()
        obj.grids = xmlelement[6].value()

        return obj

    def write_xml(self, xmlwriter, attr=None):
        """Write a RetrievalQuantity object to an ARTS XML file.
        """
        if attr is None:
            attr = {}

        xmlwriter.open_tag("RetrievalQuantity", attr)
        xmlwriter.write_xml(self.maintag, {'name': 'MainTag'})
        xmlwriter.write_xml(self.subtag, {'name': 'Subtag'})
        xmlwriter.write_xml(self.subsubtag, {'name': 'SubSubtag'})
        xmlwriter.write_xml(self.mode, {'name': 'Mode'})
        xmlwriter.write_xml(self.analytical, {'name': 'Analytical'})
        xmlwriter.write_xml(self.perturbation, {'name': 'Perturbation'})
        xmlwriter.write_xml(self.grids, {'name': 'Grids'})
        xmlwriter.close_tag()
