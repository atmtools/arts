import os
from tempfile import mkstemp
from copy import deepcopy
import pyarts


class TestXsecRecord:
    ref_dir = os.path.join(os.path.dirname(__file__), "reference", "")

    def setup_method(self):
        """Create a temporary file."""
        fd, self.f = mkstemp()
        os.close(fd)

        self.xr = pyarts.xml.load(os.path.join(self.ref_dir, "xsec.xml"))

    def teardown_method(self):
        """Delete temporary file."""
        os.remove(self.f)

    def test_eq(self):
        xr2 = deepcopy(self.xr)
        assert self.xr == xr2

    def test_ne_fitcoeffs(self):
        xr2 = deepcopy(self.xr)
        xr2.fitcoeffs[0].data[0, 0] = 11.11
        assert self.xr != xr2

    def test_ne_species(self):
        xr2 = deepcopy(self.xr)
        xr2.species = "CFC12"
        assert self.xr != xr2

    def test_ne_fitminpressures(self):
        xr2 = deepcopy(self.xr)
        xr2.fitminpressures = [0]
        assert self.xr != xr2

    def test_ne_fitmaxpressures(self):
        xr2 = deepcopy(self.xr)
        xr2.fitmaxpressures = [0]
        assert self.xr != xr2

    def test_ne_fitmintemperatures(self):
        xr2 = deepcopy(self.xr)
        xr2.fitmintemperatures = [0]
        assert self.xr != xr2

    def test_ne_fitmaxtemperatures(self):
        xr2 = deepcopy(self.xr)
        xr2.fitmaxtemperatures = [0]
        assert self.xr != xr2

    def test_species(self):
        xr = pyarts.classes.XsecRecord()
        xr.species = "CFC11"
        assert str(xr.species) == "CFC11"

    def test_xarray(self):
        xa = self.xr.to_xarray()
        xr2 = pyarts.classes.XsecRecord.from_xarray(xa)
        assert self.xr == xr2

    def test_netcdf(self):
        self.xr.to_netcdf(self.f)
        xr2 = pyarts.classes.XsecRecord.from_netcdf(self.f)
        assert self.xr == xr2
