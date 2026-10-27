"""
Utility functions.

Contains miscellaneous functions that do complicated things but
are note really relevant to other code.
"""

import ast
import inspect
import os
import sys
import tokenize
import types


def _used_names(code):
    """The names a code object refers to, its free variables included."""
    return frozenset(code.co_names) | frozenset(code.co_freevars)


def _consts_match(a, b):
    """Whether two code objects compiled independently use the same constants.

    The bytecode pins the position of the literals but not their values,
    which the constants do.  Constants can be code objects themselves, and
    those compare by identity across compilations, so the comparison is
    skipped when any is present.
    """
    for consts in (a.co_consts, b.co_consts):
        if any(isinstance(c, types.CodeType) for c in consts):
            return True
    return [repr(c) for c in a.co_consts] == [repr(c) for c in b.co_consts]


def _block_compiles_to(block, code, filename):
    """Return whether a source block is the definition of a code object.

    Reading a file can find a function of the same name that is not the one
    that was executed.  The line numbers of code executed from a string only
    match the original file if that string is the unmodified file content, and
    ``__file__`` or ``sys.argv[0]`` may point at an unrelated file.  A wrong
    block would silently produce a wrong agenda, so the block is compiled on
    its own and compared to the code object that is really running.

    The block keeps its original indentation, wrapped in a dummy block, so
    that string literals inside it are not rewritten.  Free variables of
    ``code`` resolve as globals in the stand-alone block, so the names are
    compared with :func:`_used_names` and the bytecode is only compared when
    no free variables are involved.  The bytecode pins the position of the
    literals but not their values, which :func:`_consts_match` compares
    instead.  For a function with free variables, or with code objects among
    its constants, only names, arguments and local variables are compared,
    so an unrelated function at the same line of a candidate file can still
    be accepted there.

    Args:
        block: The source lines of the located block.
        code: The code object of the function that is being inspected.
        filename: The file the block was read from.

    Returns:
        bool: True if the block is the definition of ``code``.
    """

    source = "".join(block)
    if source[:1] in (" ", "\t"):
        source = "if 1:\n" + source

    try:
        module = compile(source, filename, "exec")
    except (SyntaxError, ValueError):
        return False

    for block_code in module.co_consts:
        if not isinstance(block_code, types.CodeType):
            continue
        if block_code.co_name != code.co_name:
            continue

        if (
            block_code.co_argcount == code.co_argcount
            and block_code.co_posonlyargcount == code.co_posonlyargcount
            and block_code.co_kwonlyargcount == code.co_kwonlyargcount
            and block_code.co_varnames == code.co_varnames
            and _used_names(block_code) == _used_names(code)
            and (code.co_freevars or block_code.co_code == code.co_code)
            and _consts_match(block_code, code)
        ):
            return True

    return False


def source_and_location(func):
    """Return the source, filename and first line number of a function.

    ``inspect`` can only find source code through :mod:`linecache`.  Code
    that was executed from a string has the filename ``<string>`` and no
    linecache entry, so ``inspect.getsource`` fails.  Matplotlib's plot
    directive executes example files this way, but keeps the original path
    in ``sys.argv[0]``.  Fall back to that path (or ``__file__``), locate the
    function with the code object's first line number, and accept the block
    only if it is the definition of that code object.

    Args:
        func: The function whose source is requested.  It is unwrapped with
            :func:`inspect.unwrap` first, as :func:`inspect.getsourcelines`
            does, or the fallback would return the source of a wrapper
            applied with :func:`functools.wraps` instead of the wrapped
            function's.  Comparing against the wrapper's code object would
            not notice, the wrapper's own source matches it exactly.

    Raises:
        OSError: If the source cannot be found or no located block is the
            requested function.

    Returns:
        tuple: ``(source, filename, first_line)``.
    """
    func = inspect.unwrap(func)

    try:
        lines, first_line = inspect.getsourcelines(func)
        return "".join(lines), inspect.getfile(func), first_line
    except OSError as err:
        original_error = err

    code = func.__code__
    candidates = [code.co_filename, func.__globals__.get("__file__")]
    if sys.argv:
        candidates.append(sys.argv[0])

    for filename in candidates:
        if not filename or not os.path.isfile(filename):
            continue

        try:
            with tokenize.open(filename) as file:
                file_lines = file.readlines()

            first_line = code.co_firstlineno - 1
            if not 0 <= first_line < len(file_lines):
                continue

            block = inspect.getblock(file_lines[first_line:])
            tree = ast.parse(unindent("".join(block)))
        except Exception:
            # Any failure here only means that this candidate cannot be read
            # as this function, which is reported together with the original
            # error once every candidate is spent.
            continue

        if len(tree.body) != 1 or not isinstance(
            tree.body[0], (ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            continue

        if _block_compiles_to(block, code, os.fspath(filename)):
            return "".join(block), os.fspath(filename), first_line + 1

    raise OSError(
        f"could not get source code for {func.__qualname__!r}"
    ) from original_error


def unindent(source):
    """Unindent source code.

    Determines indent level of the first line and unindents
    all lines by the found indent level.

    Args:
        source: The source code to unindent as a string as
            obtained for example from inspect.getsource.

    Raises:
        Exception: If the non-whitespace characters are detected
            in the characters that are stripped off the code lines.

    Returns:
        str: The unindented source code.
    """

    if not type(source) == str:
        raise Exception("Argument must be a string.")

    lines = source.splitlines()

    if len(lines) < 1:
        return ""

    n_indent = len(lines[0]) - len(lines[0].lstrip())

    lines_new = []
    for i, l in enumerate(lines):
        stripped = l[:n_indent]
        if len(stripped.lstrip()) > 0:
            err = "Error when unindenting source code. Stripped characters" \
                  + stripped + " in line " + str(i) + " are non-whitespace " \
                  + " characters."
            raise Exception(err)
        lines_new += [l[n_indent:]]
    return "\n".join(lines_new)
