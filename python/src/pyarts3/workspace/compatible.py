import pyarts3 as pyarts

wsgs = list(pyarts.arts.globals.workspace_groups().keys())
wsgs.sort()


def compat(*args, **kwargs):
    """A list of the compatible ARTS groups

    Parameters
    ----------
    args : object
        Any python object
    kwargs : object
        Any python objects

    Returns
    -------
    out : list
        All possible ARTS workspace groups that can be constructed from var

    """
    out = []
    for group in wsgs:
        if "Any" == group:
            continue
        try:
            g = eval(f"pyarts.arts.{group}")
            g(*args, **kwargs)
            out.append(g)
        except TypeError:
            pass
        except RuntimeError:
            pass
    return out
