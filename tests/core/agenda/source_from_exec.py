"""Source inspection for functions defined by exec'ing a file.

Matplotlib's plot directive executes example files as strings.  The code
objects then have the filename ``<string>``, so ``inspect.getsource`` fails
and the decorators have to fall back to the path kept in ``sys.argv[0]``.
"""

import contextlib
import functools
import inspect
import pathlib
import sys
import tempfile

import pyarts3 as pyarts
from pyarts3.workspace.utility import source_and_location

SOURCE = """
import functools
import pyarts3 as pyarts

ws = pyarts.Workspace()


def build_agenda():
    @pyarts.arts_agenda
    def from_exec(ws):
        ws.abs_bandsFromModelState()

    return from_exec


def wrap():
    def nested(ws):
        ws.abs_bandsFromModelState()

    return nested


@pyarts.workspace.callback_operator(inputs=[], outputs=[])
def callback_from_exec():
    pass


@pyarts.in_parallel(ws=ws)
def in_parallel_from_exec(ws):
    ws.abs_speciesSet(species=["H2O"])


def undecorated(ws):
    ws.abs_bandsFromModelState()


def wraps_agenda(func):
    @functools.wraps(func)
    def inner(ws):
        ws.abs_speciesSet(species=["WRONG"])

    return inner


from_exec = build_agenda()
nested = wrap()
wrapped = wraps_agenda(undecorated)
"""


def agenda_source(species, decorated):
    """A single agenda function, at the same line in every variant."""
    decorator = "@pyarts.arts_agenda\n" if decorated else ""
    return (
        "\n"
        "import pyarts3 as pyarts\n"
        "\n"
        "\n"
        f"{decorator}def agenda(ws):\n"
        f"    ws.abs_speciesSet(species={species})\n"
    )


@contextlib.contextmanager
def argv_is(path):
    """Run with sys.argv[0] set, as matplotlib's plot directive does.

    The fallback can only read the executed file while that path is known,
    which is the case while the code is executed.
    """
    old_argv = sys.argv
    sys.argv = [str(path)]
    try:
        yield
    finally:
        sys.argv = old_argv


def exec_with_argv(source, path):
    """exec() source the way the plot directive does.

    Returns:
        dict: The namespace the source was executed in.
    """
    with argv_is(path):
        namespace = {}
        exec(source, namespace)
    return namespace


def without_fast_path(func):
    """source_and_location with the inspect lookup disabled."""

    def fails(*args, **kwargs):
        raise OSError("fast path disabled")

    real_getsourcelines = inspect.getsourcelines
    inspect.getsourcelines = fails
    try:
        return source_and_location(func)
    finally:
        inspect.getsourcelines = real_getsourcelines


def passthrough(func):
    return func


@passthrough
def decorated_here(ws):
    ws.abs_bandsFromModelState()


def plain_here(ws):
    ws.abs_bandsFromModelState()


def wrapping(func):
    @functools.wraps(func)
    def inner(ws):
        ws.abs_speciesSet(species=["WRONG"])

    return inner


@wrapping
def wrapped_here(ws):
    ws.abs_bandsFromModelState()


with tempfile.TemporaryDirectory() as directory:
    directory = pathlib.Path(directory)

    source_path = directory / "from_exec.py"
    source_path.write_text(SOURCE, encoding="utf-8")

    agenda_path = directory / "agenda.py"
    agenda_path.write_text(agenda_source('["H2O"]', True), encoding="utf-8")

    decoy_path = directory / "decoy.py"
    decoy_path.write_text(agenda_source('["O2", "N2"]', False), encoding="utf-8")

    decoy2_path = directory / "decoy2.py"
    decoy2_path.write_text(agenda_source('["O2"]', True), encoding="utf-8")

    namespace = exec_with_argv(SOURCE, source_path)

    agenda = namespace["from_exec"]
    assert len(agenda.methods) == 1
    assert agenda.methods[0].name == "abs_bandsFromModelState"

    assert isinstance(namespace["callback_from_exec"], pyarts.arts.CallbackOperator)

    assert namespace["ws"].abs_species == ["H2O"]

    # The fallback finds what inspect cannot, both for a function at the top
    # level of the executed source and for an indented one.
    undecorated = namespace["undecorated"]
    nested = namespace["nested"]
    with argv_is(source_path):
        assert without_fast_path(undecorated) == (
            "def undecorated(ws):\n    ws.abs_bandsFromModelState()\n",
            str(source_path),
            undecorated.__code__.co_firstlineno,
        )

        nested_source, nested_file, nested_line = without_fast_path(nested)
    assert nested_source.startswith("    def nested(ws):")
    assert nested_file == str(source_path)
    assert nested_line == nested.__code__.co_firstlineno

    # A wrapper applied with functools.wraps must not make the fallback read
    # the source of the wrapper instead of the wrapped function, as both the
    # wrapper and the file it is defined in pass the code object comparison.
    wrapped = namespace["wrapped"]
    with argv_is(source_path):
        wrapped_source, wrapped_file, wrapped_line = without_fast_path(wrapped)
        wrapped_agenda = pyarts.arts_agenda(wrapped)
    assert wrapped_source == (
        "def undecorated(ws):\n    ws.abs_bandsFromModelState()\n"
    )
    assert wrapped_file == str(source_path)
    assert wrapped_line == undecorated.__code__.co_firstlineno
    assert [m.name for m in wrapped_agenda.methods] == ["abs_bandsFromModelState"]

    # A same-named function at the same line of an unrelated file is not
    # accepted, it would silently build the wrong agenda, ...
    try:
        exec_with_argv(agenda_source('["H2O"]', True), decoy_path)
        raise AssertionError("agenda parsed from the source of another file")
    except OSError:
        pass

    # ... neither is one on the very same lines whose body differs only in a
    # literal value, which the bytecode alone cannot tell apart, ...
    try:
        exec_with_argv(agenda_source('["H2O"]', True), decoy2_path)
        raise AssertionError("agenda parsed from a decoy differing in literals")
    except OSError:
        pass

    # ... whereas the identical source is accepted as its own file.
    agenda = exec_with_argv(agenda_source('["H2O"]', True), agenda_path)["agenda"]
    assert [m.name for m in agenda.methods][-1] == "abs_speciesSet"

# For a function in a readable file the fallback agrees with inspect, including
# for a decorator that moves the first line of the definition and for a wrapper
# applied with functools.wraps.
for local in (decorated_here, plain_here, wrapped_here):
    fast = source_and_location(local)
    slow = without_fast_path(local)
    assert fast[0] == slow[0]
    assert fast[2] == slow[2]
