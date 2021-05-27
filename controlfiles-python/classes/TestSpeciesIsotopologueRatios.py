import pyarts

arts = pyarts.workspace.Workspace()
arts.isotopologue_ratiosInitFromBuiltin()

ratios = arts.isotopologue_ratios.value

get_index = pyarts.classes.SpeciesIsotopeRecord.get_index

# Weird numeric checks based on known compile-time values of the ARTS builtin
# ratios
assert ratios[get_index("H2O", "161")] == .997317E+00
assert ratios[get_index("CO" , "28" )] == 1.97822E-03
