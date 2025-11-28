import pyarts3 as pyarts


def quantum_level_identifier():
    return pyarts.arts.QuantumLevelIdentifier("H2O-161 J 1 ElecStateLabel X")


def quantum_identifier():
    return pyarts.arts.QuantumIdentifier("H2O-161 J 1 2 ElecStateLabel X A")


def quantum_level():
    return pyarts.arts.QuantumLevel("J 1 ElecStateLabel X")


def quantum_state():
    return pyarts.arts.QuantumState("J 1 2 ElecStateLabel X A")


def string():
    return pyarts.arts.String("J 1 2 ElecStateLabel X A")


def species_isotope():
    return pyarts.arts.SpeciesIsotope("H2O-161")


def species_tag():
    return pyarts.arts.SpeciesTag("H2O-161")


def rational():
    return pyarts.arts.Rational("1/2")


def surface_property_tag():
    return pyarts.arts.SurfacePropertyTag("x")


def subsurface_property_tag():
    return pyarts.arts.SubsurfacePropertyTag("x")


def scattering_species_property_tag():
    return pyarts.arts.ScatteringSpeciesProperty("xyz_n")


calls = [
    quantum_level_identifier,
    quantum_identifier,
    quantum_level,
    quantum_state,
    string,
    species_isotope,
    species_tag,
    rational,
    surface_property_tag,
    subsurface_property_tag,
    scattering_species_property_tag,
]

for call in calls:
    obj = call()
    str_obj = str(obj)
    print(type(obj))
    print(f"obj      = {obj}")
    print(f"str(obj) = {str_obj}")
    print("-" * 80)
    assert str_obj == obj, f"str({obj}) != {obj}"
