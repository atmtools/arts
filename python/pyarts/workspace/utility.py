"""
Utility functions.

Contains miscellaneous functions that do complicated things but
are note really relevant to other code.
"""

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
        new_source: The unindented source code.
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
