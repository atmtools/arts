# -*- coding: utf-8 -*-
"""Testing high-level functionality in arts.xml.
"""
from os import environ
from os.path import (dirname, join)

import numpy as np
import pytest

from pyarts import xml


class TestXML:
    """Testing high-level functionality in arts.xml."""
    ref_dir = join(dirname(__file__), "reference")

    def test_load_directory(self):
        """Test loading all XML files in a directory."""
        t = xml.load_directory(self.ref_dir)
        ref = xml.load(join(self.ref_dir, 'vector.xml'))

        assert np.allclose(t['vector'], ref)

    def test_load_directory_exclude(self):
        """Test excluding files when loading directory content."""
        t = xml.load_directory(self.ref_dir, exclude=['vector.xml'])

        with pytest.raises(KeyError):
            t['vector']

    def test_load_xml_search_arts_path(self):
        """Test loading file from ARTS_DATA_PATH"""
        backup_path = environ.get("ARTS_DATA_PATH")
        environ["ARTS_DATA_PATH"] = self.ref_dir

        xml.load("vector.xml")

        if backup_path:
            environ["ARTS_DATA_PATH"] = backup_path
        else:
            environ.pop("ARTS_DATA_PATH")
