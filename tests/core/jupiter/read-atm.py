import pyarts3 as pyarts

ws = pyarts.Workspace()

ws.atmospheric_fieldRead(toa=500e3,
                         basename="planets/Jupiter/MPS/",
                         deal_with_field_component="Ignore")
