"""Runtime Python conversion and workspace dispatch for generic arguments."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts
ws = pyarts.Workspace()


def rejects(call):
    try:
        call()
    except (RuntimeError, TypeError, ValueError):
        return
    raise AssertionError("Invalid generic argument accepted")


# Inputs can convert native Python objects. Outputs must retain object identity.
for value in (1, 1.5, "test", arts.Vector([1.0, 2.0]), arts.AtmKey.temperature):
    ws.Ignore(input=value)
rejects(lambda: ws.Ignore(input=None))
rejects(lambda: ws.Touch(input=1.5))

with TemporaryDirectory() as directory:
    for source, destination in (
        (arts.Numeric(3.5), arts.Numeric(-1.0)),
        (arts.Index(7), arts.Index(-1)),
        (arts.String("shared output"), arts.String("before")),
        (arts.Vector([2.0, 3.0]), arts.Vector()),
        (arts.Matrix([[1.0, 2.0], [3.0, 4.0]]), arts.Matrix()),
    ):
        path = str(Path(directory) / "value.xml")
        ws.WriteXML(output_file_format="ascii", input=source, filename=path)
        ws.ReadXML(output=destination, filename=path)
        assert str(destination) == str(source), (source, destination)
        ws.Touch(input=destination)

    source = arts.Vector([5.0, 6.0])
    destination = arts.Vector()
    prefix = str(Path(directory) / "indexed.xml")
    ws.WriteXMLIndexed(
        output_file_format="ascii",
        file_index=3,
        input=source,
        filename=prefix,
        digits=3,
    )
    ws.ReadXMLIndexed(output=destination, file_index=3, filename=prefix, digits=3)
    np.testing.assert_array_equal(destination, source)

    # Execute the same Any-output path through the generated workspace adapter.
    ws.freq_grid = [1e9, 2e9]
    path = str(Path(directory) / "frequency.xml")
    ws.WriteXML(output_file_format="ascii", input=ws.freq_grid, filename=path)
    agenda = arts.Agenda("read_generic")
    agenda.add(arts.Method("xml_filename", arts.String(path)))
    agenda.add(
        arts.Method("ReadXML", [], {"output": "freq_grid", "filename": "xml_filename"})
    )
    ws.freq_grid = [4e9]
    agenda.execute(ws)
    np.testing.assert_array_equal(ws.freq_grid, [1e9, 2e9])

# A shared input type list needs no Python variant caster. String conversions
# must try the declared alternatives rather than settle on the generic String.
ws.jac_targetsInit()
ws.jac_targetsAddAtmosphere(target="temperature")
ws.jac_targetsAddAtmosphere(target="H2O")
ws.jac_targetsAddSurface(target="t")
assert len(ws.jac_targets.atm) == 2
assert len(ws.jac_targets.surf) == 1
rejects(lambda: ws.jac_targetsAddAtmosphere(target="not a target"))
rejects(lambda: ws.jac_targetsAddAtmosphere(target=arts.Matrix()))
rejects(lambda: ws.jac_targetsAddAtmosphere(target=None))
assert len(ws.jac_targets.atm) == 2

# Retrieval wrappers forward the same variant to the Jacobian method.
ws.RetrievalInit()
ws.RetrievalAddAtmosphere(target="temperature", matrix=np.eye(1))
ws.RetrievalAddAtmosphere(target="H2O", matrix=np.eye(1))
assert len(ws.jac_targets.atm) == 2
