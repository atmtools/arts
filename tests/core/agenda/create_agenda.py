import pyarts3 as pyarts


@pyarts.arts_agenda
def create_agenda(ws: pyarts.Workspace):
    ws.abs_bandsFromModelState()


assert [method.name for method in create_agenda.methods] == ["abs_bandsFromModelState"]
