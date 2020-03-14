# -*- coding: utf-8 -*-

"""Miscellaneous convenience functions.
"""

import ast
import functools
import operator
import os
import shutil
import subprocess
import collections
import itertools
from numbers import Number
import traceback

from warnings import warn
from functools import (partial, wraps)

import numpy as np


__all__ = [
    'deprecated',
    'extract_block_diag',
    'safe_eval',
    'unique',
    'path_append',
    'path_prepend',
    'path_remove',
    'get_time_dimensions',
    'get_time_coordinates',
    'concat_each_time_coordinate',
    'image2mpeg',
    'split_units',
    'reraise_with_stack'
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
    wrapper.__doc__ += (
        '\n\n    .. warning::\n       Function is deprecated'
        ' and will be removed in a future version.'
    )

    if new_name is not None:
        wrapper.__doc__ += f' Use :func:`{new_name}` instead.'

    if message is not None:
        wrapper.__doc__ += f'\n\n       {message}'

    return wrapper


def extract_block_diag(M, n):
    """Extract diagonal blocks from square Matrix.

    Args:
        M (np.array): Square matrix.
        n (int): Number of blocks to extract.

    Example:
        >>> foo = np.array([[ 1.,  1.,  0.,  0.],
        ... [ 1.,  1.,  0.,  0.],
        ... [ 0.,  0.,  2.,  2.],
        ... [ 0.,  0.,  2.,  2.]])
        >>> extract_block_diag(foo, 2)
        [array([[ 1.,  1.],
                [ 1.,  1.]]), array([[ 2.,  2.],
                [ 2.,  2.]])]

    """
    return [np.split(m, n, axis=1)[i] for i, m in enumerate(np.split(M, n))]

# This code, or a previous version thereof, was posted by user 'J. F.
# Sebastian' on http://stackoverflow.com/a/9558001/974555 on 2012-03-04
# and is dual-licensed under CC-BY-SA 3.0 and MIT, as confirmed at
# https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string/9558001?noredirect=1#comment76927447_9558001
# on 2017-07-07

operators = {ast.Add: operator.add,
             ast.Sub: operator.sub,
             ast.Mult: operator.mul,
             ast.Div: operator.truediv,
             ast.Pow: operator.pow,
             ast.BitXor: operator.xor,
             ast.USub: operator.neg}


def safe_eval(expr):
    """Safely evaluate string that may contain basic arithmetic
    """

    return _safe_eval_node(ast.parse(expr, mode="eval").body)


def _safe_eval_node(node):
    if isinstance(node, ast.Num):  # <number>
        return node.n
    elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
        return operators[type(node.op)](
            _safe_eval_node(node.left), _safe_eval_node(node.right))
    elif isinstance(node, ast.UnaryOp):  # <operator> <operand> e.g., -1
        return operators[type(node.op)](_safe_eval_node(node.operand))
    else:
        raise TypeError(node)
# End of snippet derived from http://stackoverflow.com/a/9558001/974555


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


def get_time_dimensions(ds):
    """From a xarray dataset or dataarray, get dimensions corresponding to time coordinates

    """

    return {k for (k, v) in ds.coords.items() if k in ds.dims and v.dtype.kind == "M"}


def get_time_coordinates(ds):
    """From a xarray dataset or dataarray, get coordinates with at least 1 time dimension

    """

    time_dims = get_time_dimensions(ds)
    return {k for (k, v) in ds.coords.items() if set(v.dims)&time_dims}


# Any commits made to this module between 2015-05-01 and 2017-03-01
# by Gerrit Holl are developed for the EC project “Fidelity and
# Uncertainty in Climate Data Records from Earth Observations (FIDUCEO)”.
# Grant agreement: 638822.  This specifically applies to the function
# concat_each_time_coordinate.
#
# All those contributions are dual-licensed under the MIT license for use
# in typhon, and the GNU General Public License version 3.

def concat_each_time_coordinate(*datasets):
    """Concatenate xarray datasets along each time coordinate

    Given two or more xarray datasets, concatenate seperately data
    variables with different time coordinates.  For example, one might
    have dimensions 'scanline' and 'calibration_cycle' that are each along
    time coordinates.  Data variables may have dimension either scanline
    or calibration_cycle or neither, but not both.  Both correspond to
    coordinates a datetime index.  Ordinary xarray.concat along either
    dimension will broadcast the other one in a way similar to repmat,
    thus exploding memory usage (test case for one FCDR HIRS granule: 89
    MB to 81 GB).  Instead, here, for each data variable, we will
    concatenate only along at most one time coordinate.

    Arguments:

        *datasets: xarray.Dataset objects to be concatenated
    """

    time_coords = get_time_coordinates(datasets[0])
    time_dims = get_time_dimensions(datasets[0])
    # ensure each data-variable has zero or one of those time coordinates
    # as dimensions
    for ds in datasets:
        if not all([len(set(v.dims) & time_coords) <= 1
                    for (k, v) in ds.data_vars.items()]):
            raise ValueError("Found vars with multiple time coords")

    new_sizes = {k: sum(g.dims[k] for g in datasets)
                     if k in time_coords
                     else datasets[0].dims[k]
                 for k in datasets[0].dims.keys()}
    # note data vars per time coordinate
    time_vars = {k: (set(v.dims)&time_coords).pop() for (k, v) in datasets[0].variables.items() if set(v.dims)&time_coords}
    time_vars_per_time_dim = {k: {vn for (vn, dn) in time_vars.items() if dn==k} for k in time_coords}
    untimed_vars = datasets[0].data_vars.keys() - time_vars.keys()

    # allocate new
    new = xarray.Dataset(
        {k: (v.dims,
             np.zeros(shape=[new_sizes[d] for d in v.dims],
                         dtype=v.dtype))
                for (k, v) in datasets[0].data_vars.items()})
    # coordinates cannot be set in the same way so need to be allocated
    # separately
    new_coords = {k: xarray.DataArray(
                        np.zeros(shape=[new_sizes[d] for d in v.dims],
                                 dtype=datasets[0][k].dtype),
                        dims=v.dims)
                    for (k, v) in datasets[0].coords.items()}

    # copy over untimed vars
    for v in untimed_vars:
        new[v].values[...] = datasets[0][v].values

    # and untimed coords
    for c in datasets[0].coords.keys() - time_coords:
        new_coords[c][...] = datasets[0].coords[c]

    # keep track of progress per time dimension
    n_per_dim = dict.fromkeys(time_coords, 0)
    # copy timed vars dataset by dataset
    for ds in datasets:
        for (v, timedim) in time_vars.items():
            ncur = n_per_dim[timedim]
            nnew_cur = ds.dims[timedim]
            if nnew_cur == 0:
                # nothing to fill, but prevent
                # https://github.com/pydata/xarray/issues/1329
                continue
            slc = {dim: slice(ncur, ncur+nnew_cur)
                        if dim==timedim else slice(None)
                   for dim in ds[v].dims}
            if v in time_coords: # time coordinate
                new_coords[v][slc] = ds[v]
            else:
                new[v].loc[slc] = ds[v]
        for timedim in time_dims:
            n_per_dim[timedim] += ds.dims[timedim]
    # copy attributes
    new.attrs.update(**datasets[0].attrs)
    for k in new.variables.keys():
        new[k].attrs.update(**datasets[0][k].attrs)
        new[k].encoding.update(**datasets[0][k].encoding)
    return new.assign_coords(**new_coords)

def undo_xarray_floatification(ds, fields=None):
    """convert floats back to ints in xarray dataset where appropriate

    When xarray opens a NetCDF file with the default decode_cf=True,
    any integer values that have a _FillValue set are converted to float,
    such that any _FillValue-set values can be set to nan.  Some datasets
    may have such _FillValue set even though they are never used.
    In this case, it may be desirable to convert those values back to
    the original dtype (which is preserved in the .encoding attribute),
    for example, when those integers are intended to be used as indices.
    This function takes an xarray Dataset, checks all the variables which
    originally have an integer dtype and a fillvalue set, and converts
    those back to int.  Optionally only a subset of those is converted.

    Use this function only when those fill values are not used.  Behaviour
    when fill values are actually used is undefined.

    Parameters:

        ds (xarray.Dataset): xarray dataset to be converted.  Will be
        copied.

        fields (Collection or None): Describes what fields shall be
            converted.  If not given or None (default), convert all fields
            that were originally ints but converted to float due to having a
            _FillValue set.  Even when given, only fields meeting those
            criteria will be converted.

    Returns:
        The same dataset but with changes as described above.
    """

    to_correct = {k for (k, v) in ds.data_vars.items()
        if v.encoding.get("dtype", np.dtype("O")).kind[0] in "ui" and
        not v.dtype.kind in "uiMm"} # don't convert datetime/deltas

    if fields is not None:
        to_correct &= fields

    ds2 = ds.copy()
    for k in to_correct:
        ds2[k] = ds[k].astype(ds[k].encoding["dtype"])
        ds2[k].encoding.update(ds[k].encoding)

    return ds2

def image2mpeg(glob, outfile, framerate=12, resolution='1920x1080'):
    """Combine image files to a video using ``ffmpeg``.

    Notes:
        The function is tested for ``ffmpeg`` versions 2.8.6 and 3.2.2.

    Parameters:
        glob (str): Glob pattern for input files.
        outfile (str): Path to output file.
            The file fileformat is determined by the extension.
        framerate (int or str): Number of frames per second.
        resolution (str or tuple): Video resolution given in width and height
            (``"WxH"`` or ``(W, H)``).

    Raises:
        Exception: The function raises an exception if the
            underlying ``ffmpeg`` process returns a non-zero exit code.

    Example:
        >>> image2mpeg('foo_*.png', 'foo.mp4')
    """
    if not shutil.which('ffmpeg'):
        raise Exception('``ffmpeg`` not found.')

    # If video resolution is given as tuple, convert it into string format
    # to directly pass it to ffmpeg later.
    if isinstance(resolution, tuple) and len(resolution) == 2:
        resolution = '{width}x{height}'.format(width=resolution[0],
                                               height=resolution[1])

    p = subprocess.run(
        ['ffmpeg',
         '-framerate', str(framerate),
         '-pattern_type', 'glob', '-i', glob,
         '-s:v', resolution,
         '-c:v', 'libx264',
         '-profile:v', 'high',
         '-crf', '20',
         '-pix_fmt', 'yuv420p',
         '-y', outfile
         ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
        )

    # If the subprocess fails, raise exception including error message.
    if p.returncode != 0:
        raise Exception(p.stderr)

def split_units(value):
    """Splits a string into float number and potential unit

    References
        Taken from https://stackoverflow.com/a/30087094

    Args:
        value: String with number and unit.

    Returns
        A tuple of a float and unit string.

    Examples:

        >>> split_units("2GB")
        (2.0, 'GB')
        >>> split_units("17 ft")
        (17.0, 'ft')
        >>> split_units("   3.4e-27 frobnitzem ")
        (3.4e-27, 'frobnitzem')
        >>> split_units("9001")
        (9001.0, '')
        >>> split_units("spam sandwhiches")
        (0, 'spam sandwhiches')
        >>> split_units("")
        (0, '')
    """
    units = ""
    number = 0
    while value:
        try:
            number = float(value)
            break
        except ValueError:
            units = value[-1:] + units
            value = value[:-1]
    return number, units.strip()


def reraise_with_stack(func):
    """Make functions include the whole stack in raised exceptions

    Notes:
        This is a decorator function.

    When using the concurrent.futures module, the original traceback message
    gets lost, which makes it difficult to debug. This decorator solves the
    problem.

    References:
        Taken from https://stackoverflow.com/a/29357032.
    """

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            traceback_str = traceback.format_exc()
            raise Exception(
                "Error occurred. Original traceback is\n%s\n" % traceback_str
            )

    return wrapped

