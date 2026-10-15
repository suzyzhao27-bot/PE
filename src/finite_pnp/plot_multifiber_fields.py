"""Plot calculated five-row aqueous field arrays, not an SEM reconstruction."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def make_plot(datafile,pngfile,equilibrium=False):
    report=json.loads((OUT/"five_rows_chemical_refinement.json").read_text())
    case=report["cases"][-1]
    data=np.load(OUT/datafile)
    mask=data["mask"];h=float(data["grid_nm"])
    H,W=mask.shape[0]*h,mask.shape[1]*h
    fields=[]
    for key in ("potential_mv","c_plus_mol_m3","c_minus_mol_m3"):
        a=np.full(mask.shape,np.nan);a[mask]=data[key];fields.append(a)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
    fig,axes=plt.subplots(3,1,figsize=(8.8,9.2),sharex=True,sharey=True)
    fig.subplots_adjust(left=.11,right=.90,bottom=.25,top=.85,hspace=.22)
    both=np.r_[data["c_plus_mol_m3"],data["c_minus_mol_m3"]]
    for i,(ax,field,title) in enumerate(zip(axes,fields,("Aqueous potential (mV)","Na+ (mol/m3 = mM)","Cl- (mol/m3 = mM)"))):
        cmap=plt.get_cmap("viridis").copy();cmap.set_bad("#c6c6c6")
        opts=dict(vmin=float(np.min(both)),vmax=float(np.max(both))) if i else {}
        im=ax.imshow(field.T,origin="lower",extent=[0,H,0,W],interpolation="nearest",aspect="equal",cmap=cmap,**opts)
        ax.set_title(title,fontsize=11,loc="left")
        ax.set_ylabel("Periodic x (nm)")
        fig.colorbar(im,ax=ax,fraction=.025,pad=.025)
    axes[-1].set_xlabel("Axial z (nm); left reservoir to right reservoir")
    heading="Static 1 mM aqueous equilibrium" if equilibrium else "30 mM to 1 mM open-circuit aqueous PNP"
    fig.text(.11,.975,"Finite 2D PNP: five-row unsupported fiber model",fontsize=14,va="top")
    fig.text(.11,.94,heading,fontsize=11,va="top")
    fig.text(.11,.91,"Assumed geometry and chemistry; not PE measurements",fontsize=11,fontweight="bold",va="top")
    sigma=case["chemical_conductance"]["wall_sigma_mean_equilibrium_c_m2"] if equilibrium else case["chemical_transport"]["wall_sigma_mean_c_m2"]
    captions=[
        "5 staggered fibers: side 26.6 nm; span 239.4 nm; period 79.8 nm; external water 53.2 nm/side.",
        "Grey solid = model fibers; local fluid velocity = 0. Actual 3D/SEM network is not reconstructed.",
        "Assumed sites = 0.05/nm2; K_Na = 0; K_Cl = 10000/M.",
        f"Calculated mean wall sigma = {sigma:.7f} C/m2; chemical Gauss water-side boundary.",
        "PE internal dielectric field, H/OH, DI and electrode chemistry are omitted."
    ]
    for i,line in enumerate(captions):
        fig.text(.11,.205-i*.027,line,fontsize=9,va="top")
    path=OUT/pngfile;fig.savefig(path,dpi=250);plt.close(fig)
    return dict(field_npz=datafile,field_sha256=sha(OUT/datafile),png=pngfile,png_sha256=sha(path),
                caption=captions,source_case=case["chemical_conductance"] if equilibrium else case["chemical_transport"],
                data_is_actual_solver_output=True,visual_qa="pending image inspection")


def main():
    plots=[make_plot("five_rows_chemical_30_1_fields.npz","five_rows_PNP_phi_Na_Cl.png"),
           make_plot("five_rows_chemical_static_1mM_fields.npz","five_rows_static_phi_Na_Cl.png",True)]
    receipt=dict(input_result_sha256=sha(OUT/"five_rows_chemical_refinement.json"),
                 script_sha256=sha(Path(__file__)),solver_sha256=sha(ROOT/"solver.py"),plots=plots)
    (OUT/"five_rows_plot_receipt.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"plots":[p["png"] for p in plots]}))


if __name__=="__main__":
    main()
