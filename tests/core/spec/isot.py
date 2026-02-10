import pyarts3 as pyarts

x = pyarts.arts.SpeciesIsotopologueRatios.builtin().valueless_isotopes()

for v in x:
    print("Missing builtin isotopologue ratio for", v)

assert len(x) == 0, "Missing values"
