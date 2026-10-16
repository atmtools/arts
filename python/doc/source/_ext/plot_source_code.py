"""Make ``inspect.getsource`` work inside ``plot`` directive code.

matplotlib's plot directive ``exec``'s the example code from a string, so the
resulting code objects carry the filename ``<string>`` and ``inspect`` cannot
retrieve the source of functions defined therein.  ARTS's ``arts_agenda`` and
``callback_operator`` decorators parse the source of the decorated function
with ``inspect``, so examples using them fail under the directive with
"could not get source code".

This extension replaces ``plot_directive._run_code`` with a version that
compiles the code with the real source filename.  The file exists on disk, so
``inspect`` finds its source through ``linecache`` as usual.
"""

import contextlib
import io
import os
import sys
import traceback

from matplotlib import cbook
from matplotlib.backend_bases import FigureManagerBase
from matplotlib.sphinxext import plot_directive


def _run_code(code, code_path, ns=None, function_name=None):
    """Same as ``plot_directive._run_code``, but compiling with the real filename."""
    pwd = os.getcwd()
    if plot_directive.setup.config.plot_working_directory is not None:
        os.chdir(plot_directive.setup.config.plot_working_directory)
    elif code_path is not None:
        os.chdir(os.path.abspath(os.path.dirname(code_path)))

    with cbook._setattr_cm(sys, argv=[code_path], path=[os.getcwd(), *sys.path]), \
            contextlib.redirect_stdout(io.StringIO()):
        try:
            if ns is None:
                ns = {}
            if not ns:
                if plot_directive.setup.config.plot_pre_code is None:
                    exec("import numpy as np\n"
                         "from matplotlib import pyplot as plt\n", ns)
                else:
                    exec(str(plot_directive.setup.config.plot_pre_code), ns)
            if "__main__" in code:
                ns['__name__'] = '__main__'

            # Patch out non-interactive show() to avoid triggering a warning.
            with cbook._setattr_cm(FigureManagerBase, show=lambda self: None):
                filename = code_path if code_path else "<string>"
                exec(compile(code, filename, "exec"), ns)
                if function_name is not None:
                    exec(function_name + "()", ns)
        except (Exception, SystemExit) as err:
            raise plot_directive.PlotError(traceback.format_exc()) from err
        finally:
            os.chdir(pwd)

    return ns


def setup(app):
    plot_directive._run_code = _run_code
    return {"parallel_read_safe": True, "parallel_write_safe": True}
