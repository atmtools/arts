# -*- coding: utf-8 -*-

"""This module contains functions to interact with ARTS.
"""

from arts import sensor  # noqa
from arts import xml  # noqa
from arts.common import *  # noqa

__all__ = [s for s in dir() if not s.startswith('_')]
__version__ = "@ARTS_VERSION@"
