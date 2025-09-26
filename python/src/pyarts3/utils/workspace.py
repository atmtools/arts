import pyarts3.arts as cxx


__all__ = [
    'builtin_groups',
    'builtin_group_methods',
]


def is_nb_module(obj):
    from inspect import ismodule
    return ismodule(obj)


def is_nb_class(obj):
    from inspect import isclass
    return isclass(obj)


def is_nb_function(obj):
    return isinstance(obj, type(cxx.igrf))


def is_nb_member_method(obj):
    return isinstance(obj, type(cxx.Index.conj))


def builtin_groups(mod=cxx):
    assert is_nb_module(mod), f"{mod} is not a module"
    out = []

    for x in dir(mod):
        v = getattr(mod, x)
        if is_nb_module(v):
            out.extend(builtin_groups(v))
        elif is_nb_class(v):
            out.append(v)

    return out


def builtin_group_methods(groups=None, dunders=None):
    if groups is None:
        groups = builtin_groups()
    if dunders is None:
        dunders = []

    out = {}

    for x in groups:
        for y in dir(x):
            v = getattr(x, y)
            if is_nb_member_method(v) and (not y.startswith("__") or y in dunders):
                if x not in out:
                    out[x] = {}
                out[x][y] = v

    return out
