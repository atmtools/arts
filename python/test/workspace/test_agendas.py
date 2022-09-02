"""
Test handling of agendas of the Python interface.
"""
import os
import numpy as np
import pytest
import scipy as sp
import pyarts
from pyarts.workspace import Workspace, arts_agenda
from pyarts.xml import load, save


@arts_agenda
def ppath_agenda(ws):
      ws.Ignore(ws.rte_pos2)
      ws.ppathStepByStep()


class TestAgendas:
    """
    Tests the calling of ARTS workspace methods.
    """
    def setup_method(self):
        """
        This ensures a new Workspace for every test.
        """
        self.setup_workspace()

    def setup_workspace(self):
        self.ws = Workspace()
        self.ws.execute_controlfile("artscomponents/clearsky/TestClearSky.arts")

    def test_assignment(self):
        """
        Test assignment of agendas.
        """
        ws = self.ws
        ws.ppath_agenda = ppath_agenda

    def test_include(self):
        ws = self.ws

        @arts_agenda(ws=self.ws)
        def ppath_agenda_inc(ws):
            INCLUDE(ppath_agenda)

        ws.ppath_agenda = ppath_agenda_inc

    def test_execution(self):
        """
        Test definition and execution of agendas.
        """

        self.ws.atmosphere_dim = 1

        @arts_agenda(ws=self.ws)
        def add_1(ws):
            ws.IndexAdd(ws.atmosphere_dim,
                        ws.atmosphere_dim,
                        1)
        add_1.execute(self.ws)

        assert self.ws.atmosphere_dim.value == 2

        add_1.append_agenda_methods(add_1)
        add_1.execute(self.ws)

        assert self.ws.atmosphere_dim.value == 4

        args = [self.ws.atmosphere_dim, self.ws.atmosphere_dim, 1]

        @arts_agenda(ws=self.ws)
        def add_2(ws):
            ws.IndexAdd(*args)

        add_2.execute(self.ws)

        assert self.ws.atmosphere_dim.value == 5

    # This test gets stuck only on the CI build servers
    # Reason is not clear and it's nearly impossible to debug
    # Therefore it is commented out for now
    #
    # def test_callback(self):
        # """
        # Test callbacks by re-implementing iy_space_agenda in Python and
        # comparing results of yCalc.
        # """
        # # z_ppath = []

        # self.ws.yCalc()
        # y_old = np.array(self.ws.y.value)

        # print("import scipy")
        # import scipy.constants as c

        # print("Define agenda")
        # @arts_agenda(ws=self.ws, allow_callbacks=True)
        # def space_agenda(ws):
        #     # Since everything happens in Python we need
        #     # to tell ARTS that we are using all in and outputs.
        #     ws.Ignore(ws.f_grid)
        #     ws.Ignore(ws.rtp_pos)
        #     ws.Ignore(ws.rtp_los)
        #     ws.Touch(ws.iy)

        #     # Temperatures and frequency
        #     t = 2.735
        #     f = ws.f_grid.value

        #     # Compute radiances
        #     c1 = 2.0 * c.h / c.c ** 2
        #     c2 = c.h / c.k
        #     b = c1 * f ** 3 / (np.exp(c2 * f / t) - 1.0)

        #     # Put into iy vector.
        #     ws.iy = np.zeros((f.size, ws.stokes_dim.value.val))
        #     ws.iy.value.value[:, 0] = b

        # print("Assign agenda")
        # # Copy ppath_agenda into workspace.
        # self.ws.iy_space_agenda = space_agenda
        # print("Call yCalc again")
        # self.ws.yCalc()

        # print("Assign y_new")
        # y_new = np.array(self.ws.y.value)

        # print("Return allclose")
        # assert(np.allclose(y_new, y_old))


    def test_callback_2(self):
        """
        Test a very complicated Python callback.
        """

        @arts_agenda(ws=self.ws, allow_callbacks=True)
        def agenda(ws):
            """
            This agenda sets a workspace variable in a very
            obscure way.
            """

            class Foo:
                def __init__(self, ws):
                    self.ws = ws

                def ooo(self):
                    self.ws.IndexSet(ws.stokes_dim, 42)

            foo = Foo(ws)
            ws.IndexSet(ws.stokes_dim, 21)
            foo.ooo()

        agenda.execute(self.ws)

    def test_unknown_wsv(self):
        """
        Ensure that an exception is thrown when an unknown WSV is
        used inside an agenda.

        This covers https://github.com/atmtools/arts/issues/368
        """
        with pytest.raises(ValueError):
            @arts_agenda(ws=self.ws)
            def my_agenda(ws):
                  ws.UnknownMethod()

    def test_starred(self):
        """
        Test expansion of starred expression.
        """
        @arts_agenda(ws=self.ws)
        def agenda(ws):
            """
            This agenda uses a starred expression.
            """
            ws.IndexSet(*[ws.stokes_dim, 42])

        self.ws.stokes_dim = 0
        agenda.execute(self.ws)
        assert self.ws.stokes_dim.value == 42

    def test_double_starred(self):
        """
        Test expansion of starred expression.
        """
        @arts_agenda(ws=self.ws)
        def agenda(ws):
            """
            This agenda uses a starred expression.
            """
            ws.IndexSet(**{"out" : ws.stokes_dim,
                           "value" : 42})

        self.ws.stokes_dim = 0
        agenda.execute(self.ws)
        assert self.ws.stokes_dim.value == 42

    def test_exception(self):
        """
        Ensure that exception is thrown when a agenda
        variable is set to an invalid value.
        """
        @arts_agenda(ws=self.ws, allow_callbacks=True)
        def propmat_clearsky_agenda(ws):
              pass
        self.ws = pyarts.workspace.Workspace()
        with pytest.raises(Exception):
              self.ws.propmat_clearsky_agenda = propmat_clearsky_agenda
              
    def test_multiple_workspace_defaults(self):
        """
        Tests that we can reliably use multiple workspaces with default GINs
        """
        @pyarts.workspace.arts_agenda
        def myagenda(ws):
            """
            A simple delayed agenda that has defaults
            """
            ws.WriteXML("ascii", [1,2,3])

        def do_something():
            """
            A simple function that has its own workspace
            """
            ws = pyarts.workspace.Workspace()
            ws.test_agenda = myagenda
            ws.AgendaExecute(ws.test_agenda)
        
        # Call twice to test multiple defaults
        do_something()
        do_something()
    
    def test_agenda_set(self):
        def get_agendas():
            wsvdata = pyarts.arts.get_wsv_data()
            out = []
            for wsv in wsvdata:
                if wsv.groupname == "Agenda":
                    out.append(str(wsv.name))
            return out

        def get_options(enum_object):
            out = []
            for thing in dir(enum_object):
                if thing.startswith("__") or thing == "value": continue
                out.append(thing)
            return out

        def set_agendas(ws, agenda_string):
            options = get_options(eval(f"pyarts.arts.options.{agenda_string}DefaultOptions"))
            for enum_option in options:
                if (enum_option + ":" )not in eval(f"ws.{agenda_string}Set").__doc__:
                    raise RuntimeError(f"The option {enum_option} is not documented")
                try:
                    eval(f"ws.{agenda_string}Set(option=enum_option)")
                except RuntimeError as err:
                    print(f"Failed to parse {enum_option} of {agenda_string} with error:\n\n{err}")

        ws = pyarts.workspace.Workspace()

        agendas = get_agendas()

        for agenda in agendas: 
            if f"{agenda}Set" in dir(ws):
                set_agendas(ws, agenda)
            else:
                raise RuntimeError(f"""THERE ARE MISSING AGENDA DEFAULTS.
If you want a no-defaults version, copy-paste the following (fixing obvious hints in methods.cc):

To agenda_set.cc:
Agenda get_{agenda}(Workspace& ws, const String& option) {"{"}
  AgendaCreator agenda(ws, "{agenda}");

  using enum Options::{agenda}DefaultOptions;
  switch (Options::to{agenda}DefaultOptionsOrThrow(option)) {"{"}
    case FINAL:
      break;
  {"}"}
  
  return agenda.finalize();
{"}"}

To agenda_set.h:
Agenda get_{agenda}(Workspace& ws, const String& option);

To m_agenda_set.cc:
void {agenda}Set(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {"{"}
  out = get_{agenda}(ws, option);
{"}"}

To methods.cc:
  md_data_raw.push_back(
      create_mdrecord(NAME("{agenda}Set"),
                      DESCRIPTION(R"--(Sets *{agenda}* to a default value

Options are:
    There are currently no options, calling this function is an error.
    It only exist to enforce defaultable options for future agendas
    If you are copy-pasting this into methods.cc, dear author,
    please add one default to help us use your agenda :)
    If you do not foresee adding other options in the near-future, make this
    default the GIN_DEFAULT
)--"),
                      AUTHORS("Automatic Nonsense Name, Please Fix"),
                      OUT("{agenda}"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));
    
To arts_options.h:
/** Options for setting {agenda} --- CHANGE TO ENUMCLASS WHEN ADDING ANY OPTIONS AND REMOVE THIS PART OF THE COMMENT */
ENUMCLASS_EMPTY({agenda}DefaultOptions, char)

To py_options.cpp:
DeclareOption(Options, {agenda}DefaultOptions)

NOTE THAT THIS DESCRIPTION WAS VALID ON 1-Sep 2022, and if anythin has changed,
please update the description.  Also note that the commented out code is useful
to have only, as it allows us to easily add some agenda defaults if they do
exist in the future
""")
    def test_planet_set(self):
        def get_options(enum_object):
            out = []
            for thing in dir(enum_object):
                if thing.startswith("__") or thing == "value": continue
                out.append(thing)
            return out
        options = get_options(pyarts.arts.options.planetDefaultOptions)

        ws = pyarts.workspace.Workspace()
        for opt in options:
            assert (opt+':') in ws.PlanetSet.__doc__, f"The {opt}-option is not documented correctly"
            ws.PlanetSet(option=opt)


if __name__ == "__main__":
    ta = TestAgendas()
    ta.setup_method()
    ta.test_planet_set()

