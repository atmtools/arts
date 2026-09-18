"""Keep NumPy docstring qualifiers separate from Python type expressions."""

import re


def python_array_type(declaration):
    # Nanobind adds dtype/shape constraints that are not Python type syntax.
    # Keep the real ndarray type link instead of treating the entire constraint
    # expression as the name of a class.
    return re.sub(r'numpy\.ndarray\[dtype=[^\]]*\]', 'numpy.ndarray', declaration)


def format_signature(app, what, name, obj, options, signature, return_annotation):
    if not name.startswith('pyarts3.arts.'):
        return None
    return (python_array_type(signature) if signature else signature,
            python_array_type(return_annotation) if return_annotation else return_annotation)


def format_docstring_types(app, what, name, obj, options, lines):
    from pyarts3 import arts

    aliases = {
        'ndarray': 'numpy.ndarray',
        'sequence': 'typing.Sequence',
        'mapping': 'typing.Mapping',
        'Tuple': 'typing.Tuple',
        'Container': 'typing.Container',
    }
    for i, line in enumerate(lines):
        match = re.match(r'(:type [^:]+:|:rtype:)\s*(.*)', line)
        if not match:
            continue
        declaration = python_array_type(match[2])
        # Sphinx understands Python unions/generics; Napoleon's optional marker
        # is prose and must not become a reference to a class named optional.
        declaration = re.sub(r',\s*optional\b', ', *optional*', declaration)
        declaration = declaration.replace('array-like', ':term:`array_like`')
        declaration = declaration.replace('keyword arguments', '``keyword arguments``')
        declaration = re.sub(
            r'(?<![\w.`~])([A-Za-z_]\w*)(?![\w.])',
            lambda m: f'~pyarts3.arts.{m[0]}'
            if isinstance(getattr(arts, m[0], None), type)
            else aliases.get(m[0], m[0]),
            declaration,
        )
        lines[i] = f'{match[1]} {declaration}'


def setup(app):
    app.setup_extension("sphinx.ext.autodoc")
    app.connect('autodoc-process-signature', format_signature)
    app.connect('autodoc-process-docstring', format_docstring_types, priority=600)
    return {'version': '1', 'parallel_read_safe': True, 'parallel_write_safe': True}
