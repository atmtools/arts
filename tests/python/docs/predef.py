import pyarts3 as pyarts

specs = [str(r.name) for r in pyarts.arts.globals.all_isotopologues() if r.predef]

fulldoc =  pyarts.Workspace().propagation_matrixAddPredefined.__doc__
docstr = fulldoc.split('\n')

order_error = False
errors = []
i = 0
for x in specs:
  found_ordered = False
  for line in docstr[i:]:
    i += 1
    if x in line:
      found_ordered = True
      break

  if not found_ordered:
    if x in fulldoc:
      if order_error == False:
        order_error = True
        errors.append("The order of the predefined species are:\n" + "\n".join(specs))
        errors.append(f"Misordered {x} in docstring of propagation_matrixAddPredefined)")
    else:
      errors.append(f"Missing {x} in docstring of propagation_matrixAddPredefined, please add your new docstring")

if errors:
  raise Exception("\n" + "\n".join(errors))
