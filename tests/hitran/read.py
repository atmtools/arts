import pyarts3 as pyarts

DIR = "lines"
local_fn = "single_line.par"

ws = pyarts.Workspace()
ws.abs_bandsReadHITRAN(file=local_fn, file_formatter=['par', 'statep', 'statepp'])
ws.abs_bandsSaveSplit(dir=DIR)
