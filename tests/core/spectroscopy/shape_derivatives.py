import numpy as np
import pyarts


DO_PLOT = False
DO_CMP = True


def temp_model(t, x1, x2, x3, x4):
    return pyarts.arts.LineShapeModelParameters(t, x1, x2, x3, x4)


rtol = 1e-2
def deriv(temp, d):
    scale = rtol
    if d == 0:
        return temp, 0
    elif d == 1:
        d = temp.X0*scale
        temp.X0 += d
        return temp, d
    elif d == 2:
        d = temp.X1*scale
        temp.X1 += d
        return temp, d
    elif d == 3:
        d = temp.X2*scale
        temp.X2 += d
        return temp, d
    elif d == 4:
        d = temp.X3*scale
        temp.X3 += d
        return temp, d
    assert False, "d > 4 or d<0"


def model_var(var, scale, d):
    if "G0" == str(var):
        return deriv(temp_model("T1", scale*10000, 0.8, 999, 999), d)
    elif "D0" == str(var):
        return deriv(temp_model("T0", scale*1000, 999, 999, 999), d)
    elif "G2" == str(var):
        return deriv(temp_model("T2", scale*100, 0.5, 1, 999), d)
    elif "D2" == str(var):
        return deriv(temp_model("T3", scale*100, scale*50, 999, 999), d)
    elif "ETA" == str(var):
        return deriv(temp_model("T4", scale*0.1, 0.1, 0.5, 999), d)
    elif "FVC" == str(var):
        return deriv(temp_model("T5", scale*100, scale*1, 999, 999), d)
    elif "Y" == str(var):
        return deriv(temp_model("LM_AER", scale*1e-8, scale*1e-9, scale*5e-10, scale*3e-10), d)
    elif "G" == str(var):
        return deriv(temp_model("DPL", scale*1e-16, 0.5, scale*1e-14, -.7), d)
    elif "DV" == str(var):
        return deriv(temp_model("POLY", scale*1e-7, scale*1e-10, scale*1e-13, scale*1e-16), d)
    assert False, "Cannot find variable"


def species_model(scale, varderiv=None, d=None):
    if varderiv is None or d is None:
        assert varderiv is None and d is None, "Both or neither"

    vars = pyarts.arts.options.LineShapeVariable.get_options_as_strings()
    sels = []
    for var in vars:
        if var == str(varderiv):
            sels.append(model_var(var, scale, d))
        else:
            sels.append(model_var(var, scale, 0))
    ins = [a[0] for a in sels]
    der = [a[1] for a in sels]
    t = min(der)
    
    return (pyarts.arts.LineShapeSingleSpeciesModel(*ins),
            t if t < 0 else max(der))



scales = 1.0, 0.9
bspec = "O2", "AIR"
assert len(bspec) == len(scales)

def lsm_model(specderiv=None, varderiv=None, d=None):
    if varderiv is None or d is None or specderiv is None:
        assert varderiv is None and d is None and specderiv is None, "All or nothing"
    
    sels=[]
    dx = 0
    for i in range(len(scales)):
        if specderiv == i:
            s, dx = species_model(scales[i], varderiv, d)
            sels.append(s)
        else:
            sels.append(species_model(scales[i])[0])
    return pyarts.arts.LineShapeModel(sels), dx

def set_jacobian(ws, qid, specderiv, varderiv, d):
    ws.jacobianInit()
    ws.jacobianAddShapeCatalogParameter(line_identity=qid, variable=str(varderiv), coefficient=f"X{d-1}", species=bspec[spec])
    ws.jacobianClose()

lqid = "J 0 1 N 1 1"
locqn = pyarts.arts.QuantumNumberLocalState(lqid)
glow = locqn.state.get("J").low * 2 + 1
gupp = locqn.state.get("J").upp * 2 + 1
zee = pyarts.arts.ZeemanModel([1e-3, 1e-4])
qid = "O2-66 v 0 0 S 1 1"
full_qid = qid + ' ' + lqid
f0=100e9
i0=1e-18

for lineshapetype in ["DP", "LP", "VP", "HTP"]:
    for spec in range(len(scales)):
        for var in pyarts.arts.options.LineShapeVariable.get_options_as_strings():
            for d in range(1, 5):
                lsm, _ = lsm_model()
                lsm2, dx = lsm_model(spec, var, d)
                
                # Original valyes
                line = pyarts.arts.AbsorptionSingleLine(F0=f0, I0=i0, glow=glow, gupp=gupp,
                                                        zeeman=zee, lineshape=lsm,
                                                        localquanta=locqn)
                lines = pyarts.arts.AbsorptionLines(selfbroadening=True, bathbroadening=True,
                                                    lineshapetype=lineshapetype, quantumidentity=qid,
                                                    broadeningspecies=bspec, lines=[line])
                
                # Altered valyes
                dline = pyarts.arts.AbsorptionSingleLine(F0=f0, I0=i0, glow=glow, gupp=gupp,
                                                        zeeman=zee, lineshape=lsm2,
                                                        localquanta=locqn)
                dlines = pyarts.arts.AbsorptionLines(selfbroadening=True, bathbroadening=True,
                                                     lineshapetype=lineshapetype, quantumidentity=qid,
                                                     broadeningspecies=bspec, lines=[dline])
                
                ws = pyarts.workspace.Workspace()
                
                ws.iy_main_agendaSet(option="Emission")
                ws.ppath_agendaSet(option="FollowSensorLosPath")
                ws.PlanetSet(option="Earth")
                ws.ppath_step_agendaSet(option="GeometricPath")
                ws.water_p_eq_agendaSet()
                ws.iy_space_agendaSet()
                ws.iy_surface_agendaSet()
                ws.surface_rtprop_agendaSet(option="Blackbody_SurfTFromt_field")
                
                set_jacobian(ws, full_qid, spec, var, d)
                ws.abs_speciesSet(species=["O2"])
                ws.abs_lines_per_species = [[lines]]
                ws.propmat_clearsky_agendaAuto()
                
                ws.stokes_dim = 1
                ws.f_grid = np.linspace(-5e9, 5e9, 50) + f0
                ws.sensor_pos = [[300_000]]
                ws.sensor_los = [[180]]
                ws.z_surface = [[300]]
                ws.p_grid = np.logspace(5, -1)
                
                ws.AtmRawRead( basename = "planets/Earth/Fascod/tropical/tropical" )
                ws.AtmosphereSet1D()
                ws.AtmFieldsCalc()
                
                ws.cloudboxOff()
                ws.sensorOff()
                
                ws.atmfields_checkedCalc()
                ws.atmgeom_checkedCalc()
                ws.cloudbox_checkedCalc()
                ws.sensor_checkedCalc()
                ws.lbl_checkedCalc()
                ws.yCalc()
                
                jac = ws.jacobian.value.flatten() * 1.0
                y0 = ws.y.value * 1.0
                
                ws.jacobianOff()
                ws.abs_lines_per_species = [[dlines]]
                ws.yCalc()
                
                y1 = ws.y.value * 1.0
                jac_perturbed = (y1 - y0) / dx
                f = (ws.f_grid.value - f0) / 1e9
                
                if DO_PLOT:
                    import matplotlib.pyplot as plt
                    plt.clf()
                    plt.plot(f, jac)
                    plt.plot(f, jac_perturbed)
                    plt.title(f"Derivative {var} X{d-1} for {bspec[spec]} using {lineshapetype}")
                    plt.xlabel("Line Center Offset [GHz]")
                    plt.ylabel("Derivative [?]")
                    plt.show()
                
                if DO_CMP:
                    if (jac == jac_perturbed).all(): continue  # skip zeroes
                    
                    ### Test that things are within a percentage or so of the maximum value
                    test = np.isclose(jac, jac_perturbed, atol=rtol*max([max(jac), abs(min(jac))]), rtol=0)
                    
                    if not test.all() and DO_PLOT:
                        plt.clf()
                        plt.plot(f, jac-jac_perturbed)
                        plt.title(f"Error in Derivative for {var} X{d-1} for {bspec[spec]} using {lineshapetype}")
                        plt.ylabel("Difference error in the jacobian [?]")
                        plt.xlabel("Line Center Offset [GHz]")
                        plt.show()
                        plt.clf()
                        plt.plot(f, test)
                        plt.title(f"Error in Derivative for {var} X{d-1} for {bspec[spec]} using {lineshapetype}")
                        plt.xlabel("Good bins [-]")
                        plt.xlabel("Line Center Offset [GHz]")
                        plt.show()
                    assert test.all(), f"Error in Derivative for {var} X{d-1} for {bspec[spec]} using {lineshapetype}"
