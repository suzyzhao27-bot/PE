"""Source-backed finite-water PNP field maps; no real-film image is invented."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from solver import Geometry,solve,RG,F

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    manifest=json.loads((OUT/"run_manifest.json").read_text())
    g=Geometry(**manifest["geometry"])
    case=next(r for r in manifest["transport"] if r["scenario"]=="HO_reference_sigma_frozen" and r["c_left_mol_m3"]==30.)
    result,y,e,m,h=solve(case["c_left_mol_m3"],case["c_right_mol_m3"],case["sigma_c_m2"],geometry=g)
    if abs(result["common_salt_flux_mol_m2_s"]-case["common_salt_flux_mol_m2_s"])>1e-12:
        raise ValueError("Field solve does not reproduce frozen output case")
    data=OUT/"HO_reference_sigma_frozen_30_1_fields.npz"
    np.savez_compressed(data,mask=m.mask,ids=m.ids,potential_mv=e["u"]*RG*298.15/F*1000,
                        c_plus_mol_m3=e["cp"],c_minus_mol_m3=e["cm"],
                        local_sigma_c_m2=e["sigma_local"],grid_nm=m.h)
    fields=[]
    for values in (e["u"]*RG*298.15/F*1000,e["cp"],e["cm"]):
        grid=np.full(m.mask.shape,np.nan);grid[m.mask]=values;fields.append(grid)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
    fig,axes=plt.subplots(1,3,figsize=(10.8,6.6),sharex=True,sharey=True)
    fig.subplots_adjust(left=.075,right=.97,bottom=.28,top=.82,wspace=.35)
    cmin=min(float(np.min(e["cp"])),float(np.min(e["cm"])))
    cmax=max(float(np.max(e["cp"])),float(np.max(e["cm"])))
    titles=["Aqueous potential (mV)","Na+ (mol/m3 = mM)","Cl- (mol/m3 = mM)"]
    for i,(ax,field,title) in enumerate(zip(axes,fields,titles)):
        kw={"vmin":cmin,"vmax":cmax} if i else {}
        im=ax.imshow(field,origin="lower",extent=[0,m.width_m/1e-9,0,m.height_nm],
                     cmap="viridis",interpolation="nearest",aspect="equal",**kw)
        xs=(m.width_m/1e-9-g.fiber_side_nm)/2
        zs=(m.height_nm-g.fiber_side_nm)/2
        ax.add_patch(Rectangle((xs,zs),g.fiber_side_nm,g.fiber_side_nm,
                               facecolor="#d0d0d0",edgecolor="#222",linewidth=.7))
        ax.set_title(title,fontsize=11)
        ax.set_xlabel("Periodic transverse x (nm)")
        fig.colorbar(im,ax=ax,fraction=.055,pad=.035)
    axes[0].set_ylabel("Axial z (nm): 30 mM to 1 mM")
    fig.text(.075,.96,"Finite 2D aqueous PNP around a representative unsupported fiber layer",fontsize=14,va="top")
    fig.text(.075,.915,"Hypothetical geometry and frozen charge; not real PE measurements",fontsize=12,fontweight="bold",va="top")
    captions=[
        f"Square fiber side = {g.fiber_side_nm:g} nm; period = {m.width_m/1e-9:g} nm; gap = {(g.period_fiber_sides-1)*g.fiber_side_nm:g} nm; external water = {g.external_fiber_sides_each*g.fiber_side_nm:g} nm per side.",
        f"Given sigma = {case['sigma_c_m2']:.8f} C/m2; local u(r,z) = 0. The 230 nm full network is not reconstructed.",
        f"Open ionic current; V(right)-V(left) = {result['right_minus_left_theoretical_aqueous_voltage_mv']:.4f} mV. This excludes electrode/liquid-junction readings.",
        "Origin reference sigma is frozen into a different geometry; no coupled H/OH chemistry, actual wet gap, or DI calculation."
    ]
    for i,line in enumerate(captions):
        fig.text(.075,.20-i*.037,line,fontsize=9,va="top")
    png=OUT/"finite_PNP_fields.png"
    fig.savefig(png,dpi=250);plt.close(fig)
    receipt=dict(input_manifest_sha256=sha(OUT/"run_manifest.json"),
                 solver_sha256=sha(ROOT/"solver.py"),script_sha256=sha(Path(__file__)),
                 fields_npz_sha256=sha(data),output_png_sha256=sha(png),
                 source_case=result,caption=captions,
                 heatmap_data_source="actual numerical field arrays; masked rectangle is representative impermeable fiber, not a membrane photograph",
                 visual_qa="pending actual PNG inspection")
    (OUT/"field_plot_receipt.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"output":str(png),"sha256":receipt["output_png_sha256"]}))


if __name__=="__main__":
    main()
