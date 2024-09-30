# -*- coding: utf-8 -*-
"""Miscellaneous convenience functions.
"""

import os
from functools import partial, wraps
from tempfile import mkstemp
from warnings import warn

__all__ = [
    'deprecated',
    'unique',
    'path_append',
    'path_prepend',
    'path_remove',
    'TempFileHandler',
]


def deprecated(func=None, new_name=None, message=None):
    """Decorator which can be used to mark functions as deprecated.

    Examples:
        Calling ``foo()`` will raise a ``DeprecationWarning``.

        >>> @deprecated
        ... def deprecated_function():
        ...     pass

        Display message with additional information:

        >>> @deprecated(message='Additional information message.')
        ... def deprecated_function():
        ...     pass
    """
    # Return partial when no arguments are passed.
    # This allows a plain call of the decorator.
    if func is None:
        return partial(deprecated, new_name=new_name, message=message)

    # Build warning message (with optional information).
    msg = f'\nCall to deprecated function `{func.__name__}`.'

    if new_name is not None:
        msg += f' Use `{new_name}` instead.'

    if message is not None:
        msg += f'\n{message}'

    # Wrapper that prints the warning before calling the deprecated function.
    @wraps(func)
    def wrapper(*args, **kwargs):
        warn(msg, category=DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)

    if wrapper.__doc__ is None:
        wrapper.__doc__ = ''

    # Lines added to the docstring need to be indented with four spaces as this
    # is technically the case for all lines except the first one.
    # If this is not done, the Sphinx build will produce wrong results as the
    # relative indentation of the added lines and the original docstring does
    # not match.
    wrapper.__doc__ += ('\n\n    .. warning::\n       Function is deprecated'
                        ' and will be removed in a future version.')

    if new_name is not None:
        wrapper.__doc__ += f' Use :func:`{new_name}` instead.'

    if message is not None:
        wrapper.__doc__ += f'\n\n       {message}'

    return wrapper


def unique(seq):
    """Remove duplicates from list whilst keeping the original order

    Notes:
        If you do not care about keeping the order, use this code:
        >>>> list(set([0, 5, 1, 2, 0, 3, 1,]))
        [0, 1, 2, 3, 5]

    This code is taken from https://stackoverflow.com/a/480227.

    Args:
        seq: A sequence (list, etc.) of elements.

    Returns:
        A list with unique items with original order.

    Examples:
        >>>> unique([0, 5, 1, 2, 0, 3, 1,])
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def path_append(dirname, path='PATH'):
    """Append a directory to environment path variable.

    Append entries to colon-separated variables (e.g. the system path).
    If the entry is already in the list, it is moved to the end.
    A path variable is set, if not existing at function call.

    Parameters:
        dirname (str): Directory to add to the path.
        path (str): Name of the path variable to append to.
            Defaults to the system path 'PATH'.
    """
    if path in os.environ:
        dir_list = os.environ[path].split(os.pathsep)
        if dirname in dir_list:
            dir_list.remove(dirname)
        dir_list.append(dirname)
        os.environ[path] = os.pathsep.join(dir_list)
    else:
        os.environ[path] = dirname


def path_prepend(dirname, path='PATH'):
    """Prepend a directory to environment path variable.

    Append entries to colon-separated variables (e.g. the system path).
    If the entry is already in the list, it is moved to the end.
    A path variable is set, if not existing at function call.

    Parameters:
        dirname (str): Directory to add to the path.
        path (str): Name of the path variable to append to.
            Defaults to the system path 'PATH'.
    """
    if path in os.environ:
        dir_list = os.environ[path].split(os.pathsep)
        if dirname in dir_list:
            dir_list.remove(dirname)
        dir_list.insert(0, dirname)
        os.environ[path] = os.pathsep.join(dir_list)
    else:
        os.environ[path] = dirname


def path_remove(dirname, path='PATH'):
    """Remove a directory from environment path variable.

    Remove entries from colon-separated variables (e.g. the system path).
    If the path variable is not set, nothing is done.

    Parameters:
        dirname (str): Directory to add to the path.
        path (str): Name of the path variable to append to.
            Defaults to the system path 'PATH'.
    """
    if path in os.environ:
        dir_list = os.environ[path].split(os.pathsep)
        dir_list.remove(dirname)
        os.environ[path] = os.pathsep.join(dir_list)


class TempFileHandler:
    """Context manager for creating and deleting temporary files.

    This class is a context manager that creates a temporary file and
    deletes it when the context is exited. The file is automatically
    removed if an exception occurs.

    Parameters:
        prefix (str): Optional prefix for the temporary file name.
        suffix (str): Optional suffix for the temporary file name.
        dir (str): Optional directory for the temporary file.
    """
    def __init__(self, prefix="pyarts", suffix=".tmp", dir=None):
        self.prefix = prefix
        self.suffix = suffix
        self.dir = dir
        self.filename = None

    def __enter__(self):
        # Generate a unique filename
        (fd, self.filename) = mkstemp(
            suffix=self.suffix, prefix=self.prefix, dir=self.dir
        )
        os.close(fd)
        return self.filename

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Remove the file if it exists
        if self.filename and os.path.exists(self.filename):
            os.remove(self.filename)
        return False  # False means to propagate exceptions, if any occurred
