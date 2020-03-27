from pyarts.workspace import Workspace

ws = Workspace()

ws.StringCreate("mystring")
ws.StringSet(ws.mystring, "Hello World!")
ws.Print(ws.mystring, 0)
