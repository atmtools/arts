# -*- coding: utf-8 -*-
"""Common functions for the pyarts package.
"""
import collections
import os
import shutil
import subprocess

from pyarts.environment import environ
from pyarts.utils import path_append


__all__ = [
    'run_arts',
    ]


def run_arts(controlfile=None, arts=None, writetxt=False,
             ignore_error=False, **kwargs):
    """Start an ARTS Simulation.

    Parameters:
        controlfile (str): Path to the ARTS controlfile.
        arts (str): Path to the arts executable.
        writetxt (bool): Write stdout and stderr to ASCII files.
        ignore_error (bool): If set to True, erros during the ARTS run do not
            result in an exception (default is False).
        **kwargs: Additional command line arguments passed as keyword.
            See `arts --help` for more details.

    Returns:
        Named tuple containing the fields stdout, stderr and retcode.

    Examples:
        Run a simple ARTS job and set the output directory
        and the report level:

        >>> run_arts('foo.arts', outdir='bar', reporting='020')

        If a keyword is set to True it is added as flag.
        Show the ARTS help message:

        >>> run_arts(help=True)

    """
    # If path to the ARTS exectuable is not passed explicitly, construct it
    # from the ARTS_BUILD_PATH. Its actual existence is checked later.
    if arts is None and environ.get('ARTS_BUILD_PATH') is not None:
        arts = os.path.join(environ['ARTS_BUILD_PATH'], 'src', 'arts')
    # Try 'arts' as a fallback, maybe it is in the user's PATH.
    elif arts is None:
        arts = 'arts'

    # Append ARTS_INCLUDE_PATH and ARTS_DATA_PATH to the user's environment.
    if environ.get('ARTS_INCLUDE_PATH') is not None:
        path_append(environ.get('ARTS_INCLUDE_PATH'), path='ARTS_INCLUDE_PATH')

    if environ.get('ARTS_DATA_PATH') is not None:
        path_append(environ.get('ARTS_DATA_PATH'), path='ARTS_DATA_PATH')

    if not shutil.which(arts):
        raise Exception('ARTS executable not found at: {}'.format(arts))

    if controlfile is None:
        controlfile = ''
    elif not os.path.exists(controlfile):
        err_msg = 'Controlfile not found at: {}'.format(controlfile)
        raise FileNotFoundError(err_msg)

    opts = []
    for kw, arg in kwargs.items():
        if isinstance(arg, bool) and arg is True:
            if len(kw) == 1:
                opts.append('-{}'.format(kw))
            else:
                opts.append('--{}'.format(kw))
        elif len(kw) == 1:
            opts.append('-{}{}'.format(kw, arg))
        else:
            opts.append('--{}={}'.format(kw, arg))

    # Run ARTS job and redirect output.
    p = subprocess.run([arts, *opts, controlfile],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       universal_newlines=True
                       )

    # Write ARTS output and error to ASCII file.
    if writetxt:
        if controlfile.endswith('.arts'):
            outfile = controlfile.replace('.arts', '.out')
            errfile = controlfile.replace('.arts', '.err')
        else:
            outfile = 'arts.out'
            errfile = 'arts.err'

        for key in ['outdir', 'o']:
            if key in kwargs:
                outfile = os.path.join(kwargs[key], outfile)
                errfile = os.path.join(kwargs[key], errfile)

        with open(outfile, 'w') as out, open(errfile, 'w') as err:
            out.write(p.stdout)
            err.write(p.stderr)

    # Throw exception if ARTS run failed.
    if p.returncode != 0 and ignore_error is not True:
        raise Exception('ARTS run failed:\n{}'.format(p.stderr))

    # Store ARTS output in namedtuple.
    arts_out = collections.namedtuple(
        'ARTS_output',
        ['stdout', 'stderr', 'retcode']
    )

    return arts_out(stdout=p.stdout, stderr=p.stderr, retcode=p.returncode)

