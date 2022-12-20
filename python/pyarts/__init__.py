# -*- coding: utf-8 -*-

"""This module contains functions to interact with ARTS.
"""
from pathlib import Path

from pyarts import sensor  # noqa
from pyarts import xml  # noqa
from pyarts import arts  # noqa
from pyarts import arts_ext  # noqa
from pyarts import plots  # noqa
from pyarts import workspace  # noqa
from pyarts.common import *  # noqa
from pyarts import hitran  # noqa
from pyarts import cat  # noqa

__all__ = [s for s in dir() if not s.startswith('_')]
__version__ = "@ARTS_VERSION@"
version = __version__
ARTS_BINARY_DIR = Path("@CMAKE_BINARY_DIR@")
