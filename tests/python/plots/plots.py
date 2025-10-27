import pyarts3
import inspect

builtin = [v.__name__ for v in pyarts3.utils.builtin_groups()]
routines = [r for r in dir(pyarts3.plots) if not r.startswith("_")]
ignore = ["plot", "common"]

for ig in ignore:
    if ig in routines:
        routines.remove(ig)

errors = []
for plotting_routine in routines:
    if plotting_routine not in builtin:
        errors.append(
            f"The plotting routine {plotting_routine} is not a builtin group, it should not be where it is in the pyarts3 structure")
        continue

    v = getattr(pyarts3.plots, plotting_routine)

    if not hasattr(v, "plot"):
        errors.append(
            f"The plotting routine {plotting_routine} does not have a plot() function")
        continue

    inspection = inspect.signature(v.plot)

    for i, param in enumerate(inspection.parameters):
        par = inspection.parameters[param]

        if i == 0:
            if par.kind != inspect.Parameter.POSITIONAL_OR_KEYWORD:
                errors.append(
                    f"The parameter '{param}' of pyarts3.plots.{plotting_routine}.plot must be positional and keyword accessible")
            continue

        if i == 1 and par.kind == inspect.Parameter.VAR_POSITIONAL:
            # fine, this is also allowed
            continue

        if i == len(inspection.parameters) - 1 and par.kind == inspect.Parameter.VAR_KEYWORD:
            # fine, this is also allowed
            continue

        if par.kind != inspect.Parameter.KEYWORD_ONLY or isinstance(par.default, inspect.Parameter.empty):
            print(par.default, par.kind)
            errors.append(
                f"All but the first parameter of pyarts3.plots.{plotting_routine}.plot must be keyword only and have default values.  The parameter '{param}' breaks this rule.  If you intend to have optional unnamed parameters, use '*args*' to separate the first parameter.  If you intend to pass along optional keyword information, use '**kwargs' as the last argument.")
            continue

if errors:
    raise RuntimeError(f"Errors found in plotting routines:\n{'\n'.join(errors)}")

print(f"The plotting routines these types are at first glance valid:\n\n{'\n'.join(routines)}")
