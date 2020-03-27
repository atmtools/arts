"""Special test case which tries to load all files in ARTS XML data."""
import os

import pytest

from pyarts import xml
from pyarts import environment


def collect_xml_files():
    """Collect all XML files in the ARTS XML data tree."""
    for d in environment.ARTS_DATA_PATH.split(os.path.pathsep):
        for root, _, filenames in os.walk(d):
            for filename in filenames:
                if filename.endswith(('.xml', '.xml.gz')):
                    yield os.path.join(root, filename)


@pytest.mark.slow
@pytest.mark.skipif(environment.ARTS_DATA_PATH is None,
                    reason='ARTS_DATA_PATH not set.')
@pytest.mark.parametrize('xmlfile', collect_xml_files())
def test_load_arts_xml_data(xmlfile):
    """Try to load all XML files in ARTS_DATA_PATH.

    Search for XML files in ARTS_DATA_PATH. If files are found, try to load
    them. It is just checked, if xml.load runs without exception.
    """
    xml.load(xmlfile)
    pass
