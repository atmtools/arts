import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.StringCreate("s")
ws.StringSet(ws.iy_unit, "  Global value of iy_unit")
# Create a custom ArrayOfAgenda
ws.ArrayOfAgendaCreate("myagendaarray")
# Test appending to the array
@arts_agenda
def myagendaarray(ws):
    ws.StringSet(ws.s, "  First agenda in array")
    ws.Print(ws.s, 0)
    ws.Print(ws.iy_unit, 0)
    ws.Print(ws.agenda_array_index, 0)


ws.Append(ws.myagendaarray, myagendaarray)


@arts_agenda
def myagendaarray(ws):
    ws.StringSet(ws.s, "  Second agenda in array")
    ws.Print(ws.s, 0)
    ws.Print(ws.iy_unit, 0)
    ws.Print(ws.agenda_array_index, 0)


ws.Append(ws.myagendaarray, myagendaarray)

# Create new agenda
ws.AgendaCreate("myagenda")


@arts_agenda
def myagenda(ws):
    ws.StringSet(ws.s, "  Third agenda in array")
    ws.Print(ws.s, 0)
    ws.Print(ws.iy_unit, 0)
    ws.Print(ws.agenda_array_index, 0)


ws.myagenda = myagenda

# Append new agenda to custom array
ws.Append(ws.myagendaarray, ws.myagenda)
# Copy agenda array to WSV
ws.Copy(ws.test_agenda_array, ws.myagendaarray)
# Append another agenda directly to WSV
@arts_agenda
def myagenda(ws):
    ws.StringSet(ws.s, "  Forth agenda in array")
    ws.Print(ws.s, 0)
    ws.Print(ws.iy_unit, 0)
    ws.Print(ws.agenda_array_index, 0)


ws.myagenda = myagenda

ws.Append(ws.test_agenda_array, ws.myagenda)
# Execute agenda array manually
ws.IndexSet(ws.agenda_array_index, 0)
ws.ArrayOfAgendaExecute(agendas=ws.test_agenda_array)
ws.IndexSet(ws.agenda_array_index, 1)
ws.ArrayOfAgendaExecute(agendas=ws.test_agenda_array)
ws.IndexSet(ws.agenda_array_index, 2)
ws.ArrayOfAgendaExecute(agendas=ws.test_agenda_array)
ws.IndexSet(ws.agenda_array_index, 3)
ws.ArrayOfAgendaExecute(agendas=ws.test_agenda_array)
# Execute test_agenda_array through WSM
ws.TestArrayOfAgenda(index=0)
ws.TestArrayOfAgenda(index=1)
ws.TestArrayOfAgenda(index=2)
ws.TestArrayOfAgenda(index=3)
