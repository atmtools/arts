import matplotlib.pyplot as plt
import numpy as np
import pyarts


DO_PLOT = False
DO_CMP = True


def temp_model(t, x1, x2=0, x3=0, x4=0):
    return pyarts.arts.LineShapeModelParameters(t, x1, x2, x3, x4)


def spec_model(spec, scale, offset):
    return pyarts.arts.LineShapeSingleSpeciesModel(
        temp_model("T1", 20000 * scale + offset, 0.7*scale))


lsm = pyarts.arts.LineShapeModel([spec_model("H2O", 1, 0), spec_model(
    "O2", 0.7, 200), spec_model("BATH", 1.2, 2000)])
F0 = 100e9
I0 = 1e-16
E0 = 0
line = pyarts.arts.AbsorptionSingleLine(F0=F0, I0=I0, E0=E0, lineshape=lsm)
lines = pyarts.arts.AbsorptionLines(selfbroadening=True, bathbroadening=True, broadeningspecies=["H2O", "O2", "N2"], lines=[line], quantumidentity="H2O-161")

specs = ["H2O", "O2"]
vmrs = [0.01, 0.21]
t = 296
p = 1e5
f = np.linspace(-10e9, 10e9) + F0
perturbation = 1e-4

unpert_data = 0
perturbed = np.array([0])
for lineshapetype in ["DP", "LP", "VP", "SDVP", "HTP", "SplitLP", "SplitVP", "SplitSDVP", "SplitHTP"]:
    for spec_ind in [0, 1]:
        for perturb in [False, True]:
            lines.lineshapetype = lineshapetype
            ws = pyarts.workspace.Workspace()
            
            ws.Touch(ws.rtp_nlte)
            ws.p_grid = [p]
            ws.lat_grid = [0]
            ws.lon_grid = [0]
            ws.atmosphere_dim = 1
            ws.rtp_temperature = t
            ws.rtp_pressure = p
            ws.rtp_vmr = vmrs
            if perturb:
                ws.rtp_vmr.value[spec_ind] += perturbation
                
            ws.abs_lines_per_species = [[lines], []]
            
            ws.abs_speciesSet(species=specs)
            if perturb:
                ws.jacobianOff()
            else:
                ws.jacobianInit()
                ws.jacobianAddAbsSpecies(species="H2O", g1=[p], g2=[], g3=[], for_species_tag=0)
                ws.jacobianAddAbsSpecies(species="O2", g1=[p], g2=[], g3=[], for_species_tag=0)
                
            ws.stokes_dim = 1
            ws.f_grid = f
            ws.propmat_clearsky_agenda_checked = 1
            ws.lbl_checked = 1
            
            ws.propmat_clearskyInit()
            ws.propmat_clearskyAddLines()
            
            if perturb:
                perturbed = (1.0 * ws.propmat_clearsky.value.data.flatten() - unpert_data) / perturbation
            else:
                unpert_data = 1.0 * ws.propmat_clearsky.value.data.flatten()
                analytical = 1.0 * ws.dpropmat_clearsky_dx.value[spec_ind].data.flatten()
            
            
            if perturb:
                if (perturbed == 0).all(): continue
                x = np.allclose(perturbed/analytical, 1, rtol=1e-3)
                    
                if not x:
                    msg = f"Fail for {lineshapetype} for perturbation of {specs[spec_ind]}"
                    plt.plot(f, perturbed/analytical)
                    plt.title(msg)
                    plt.show()
                    
                if DO_CMP:
                    assert x, msg

