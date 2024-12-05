import pyarts

DIR = "lines"
local_fn = "single_line.par"

ws = pyarts.Workspace()
ws.absorption_bandsReadHITRAN(file=local_fn)
ws.absorption_bandsSaveSplit(dir=DIR)
