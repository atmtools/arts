"""
C stdout Redirection
====================

This module contains the CoutCapture class used to capture
the output from C sub

Attributes:
    cout_file: Temporary file to which the stdout is redirected if Python's
    `sys.stdout` does not have a file number, i.e. Python's output is not
    directed to the stdout file. This is the case the code is run in a
    jupyter notebook.

    If the python output is already directed to stdout, `cout_file` is none
    since no redirection of output is required.
"""
import io
import os
import sys
import tempfile

class CoutCapture():
    """
    CoutCapture class to capture output from stdout file. Implements
    the context management protocol.

    The CoutCapture class empties the file containing the redirected
    stdout upon entering and prints the file content using `print`
    upon exit.

    Usage:

    >>> with CoutCapture(ws):
    >>>     ws.Print(ws.Verbosity)

    """
    def __init__(self, ws, silent = False):
        """
        Create CoutCapture for given workspace object.

        Args:
            ws: The `Workspace` object from which the output will be
            captured.
        """
        self.ws = ws
        self.silent = silent

    def __enter__(self):
        if cout_file:
            cout_file.seek(0)
            cout_file.truncate()

    def __exit__(self, type, value, traceback):
        if cout_file and not self.silent:
            #cout_file.flush()
            cout_file.seek(0, io.SEEK_SET)
            lines = [l.decode("UTF8") for l in cout_file.readlines()]
            if lines:
                print("".join(["ARTS[{0}]: {1}".format(self.ws.ptr , l)
                               for l in lines]))

cout_file = None
try:
    sys.stdout.fileno()
except:
    cout_file = tempfile.TemporaryFile(mode='w+b')
    os.dup2(cout_file.fileno(), 1)
