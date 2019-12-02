__all__ = ['XsecRecord']


class XsecRecord:
    """:class:`XsecRecord` implements the same-named ARTS datatype.

    Contains the reference cross section data at low pressure and
    the coefficients for the broadening formula.
    """

    def __init__(self, species=None, coeffs=None, fmin=None, fmax=None,
                 refpressure=None, reftemperature=None, xsec=None,
                 tfit_slope=None, tfit_intersect=None):
        """Initialize XsecRecord object.
        """
        self.version = 1
        self.species = species
        self.coeffs = coeffs
        self.fmin = fmin
        self.fmax = fmax
        self.refpressure = refpressure
        self.reftemperature = reftemperature
        self.xsec = xsec
        self.tfit_slope = tfit_slope
        self.tfit_intersect = tfit_intersect

    def write_xml(self, xmlwriter, attr=None):
        """Write a XsecRecord object to an ARTS XML file.
        """
        # self.checksize()
        if attr is None:
            attr = {}
        attr['version'] = self.version
        xmlwriter.open_tag("XsecRecord", attr)
        xmlwriter.write_xml(self.species, {'name': 'Species'})
        xmlwriter.write_xml(self.coeffs, {'name': 'Broadening Coefficients'})
        xmlwriter.write_xml(self.fmin, {'name': 'fmin'})
        xmlwriter.write_xml(self.fmax, {'name': 'fmax'})
        xmlwriter.write_xml(self.refpressure,
                            {'name': 'Reference Pressure'})
        xmlwriter.write_xml(self.reftemperature,
                            {'name': 'Reference Temperature'})
        xmlwriter.write_xml(self.xsec, {'name': 'Cross Sections'})
        xmlwriter.write_xml(self.tfit_slope, {'name': 'Temperature Fit Slope'})
        xmlwriter.write_xml(self.tfit_intersect,
                            {'name': 'Temperature Fit Intersect'})
        xmlwriter.close_tag()

    @classmethod
    def from_xml(cls, xmlelement):
        """Loads a XsecRecord object from an xml.ElementTree.Element.
        """

        obj = cls()
        if 'version' in xmlelement.attrib.keys():
            obj.version = int(xmlelement.attrib['version'])
        else:
            obj.version = 1

        if obj.version != 1:
            raise RuntimeError(f'Unknown XsecRecord version {obj.version}')

        obj.species = xmlelement[0].value()
        obj.coeffs = xmlelement[1].value()
        obj.fmin = xmlelement[2].value()
        obj.fmax = xmlelement[3].value()
        obj.refpressure = xmlelement[4].value()
        obj.reftemperature = xmlelement[5].value()
        obj.xsec = xmlelement[6].value()
        obj.tfit_slope = xmlelement[7].value()
        obj.tfit_intersect = xmlelement[8].value()

        return obj
