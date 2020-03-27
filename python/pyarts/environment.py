# -*- coding: utf-8 -*-
"""Handling of environment variables.

Implements a mapping object for environment variables.

The handler is used like ``os.environ``:

>>> environ['ARTS_DATA_PATH']
'path/to/data'
>>> environ.get('ARTS_DATA_PATH')
'path/to/data'

In addition to the user's environment, variables can be set in the in the
configuration file (:mod:`config`) in the ``environment`` section::

    [environment]
    ARTS_BUILD_PATH: /path/to/arts/build/

Note:
    If the environment variable is set explicitly, the value set in
    the configuration file is ignored.
"""
import os


__all__ = [
    'environ',
]


class _EnvironmentHandler:
    """A mapping object for environment variables.

    See module docstring of :mod:`environment` for more information.
    """
    def __getitem__(self, key):
        try:
            # First, try to return the value from the user's environment...
            return os.environ[key]
        except :
            # TODO: Error handling.
            pass

    def __setitem__(self, key, value):
        # If an environment variable is set, pass it to the actual user's
        # environment. This ensures consistent environments for subprocesses.
        os.environ[key] = value

    def __contains__(self, item):
        """Implement membership test operators.

        Returns:
            `True` if item is in self, `False` otherwise.
        """
        try:
            self[item]
            return True
        except KeyError:
            return False

    def get(self, key, default=None):
        """D.get(k[, d]) -> D[k] if k in D, else d.d defaults to None."""
        # Try to return the value from the user's environment.
        # If the key is not set, try to find it in the config.
        # If this also fails, return a default value.
        return os.environ.get(key)


environ = _EnvironmentHandler()
