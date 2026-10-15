from pathlib import Path
import csv,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from fibre_pb import EPS0,E,R,F
O=Path(__file__).resolve().parent.parent/'outputs'
rows=list(csv.DictReader((O/'surface_scenarios.csv').open(encoding='utf-8-sig')))
profiles=list(csv.DictReader((O/'surface_radial_profiles.csv').open(encoding='utf-8-sig')))
fig,ax=plt.subplots(1,3,figsize=(13.2,4.1),constrained_layout=True)
colors=['#975330','#237e8a','#6e59a0']
for pH,col in zip([3,7,10],colors):
 p=[x for x in profiles if float(x['pH'])==pH]
 ax[0].plot([float(x['r_nm'])-13.3 for x in p],[float(x['psi_mV']) for x in p],label=f'pH {pH}',color=col)
ax[0].axhline(0,color='gray',lw=.6);ax[0].set(xlabel='Distance from solid fibre / nm',ylabel='Potential relative to reservoir / mV',title='A. Charge arises from declared adsorption')
ax[0].legend(frameon=False)
p=[x for x in profiles if float(x['pH'])==7]
ax[1].semilogy([float(x['r_nm'])-13.3 for x in p],[float(x['Na_mol_m3'])/1 for x in p],label='Na+',color='#237e8a')
ax[1].semilogy([float(x['r_nm'])-13.3 for x in p],[float(x['Cl_mol_m3'])/1 for x in p],label='Cl-',color='#bb874b')
ax[1].axhline(1,color='gray',lw=.6);ax[1].set(xlabel='Distance from solid fibre / nm',ylabel='c(r) / reservoir concentration',title='B. Mobile ions around the fibre, pH 7');ax[1].legend(frameon=False)
x=[r for r in rows if r['model']=='regulated_H_OH_hypothesis' and float(r['pH'])==7 and float(r['assumed_dG_OH_kj'])==-50]
salts=[float(r['salt_mM']) for r in x]
ax[2].semilogx(salts,[float(r['t_cation']) for r in x],'o-',color='#237e8a',label='Regulated fibre cell')
ax[2].semilogx(salts,[float(r['same_water_uncharged_t_cation']) for r in x],'--',color='gray',label='Uncharged, same water')
ax[2].set(xlabel='Background NaCl / mM',ylabel='Cation fraction of linear axial current',ylim=(0,1),title='C. Screening weakens transport preference');ax[2].legend(frameon=False,fontsize=9)
fig.suptitle('Hypothetical exterior-fibre cell: a=13.3 nm, b=39.5 nm, Gamma=0.05 / nm^2\nH/OH binding energies -20/-50 kJ/mol are assumptions, not measured PE parameters',fontsize=11)
fig.savefig(O/'surface_regulation_and_selectivity.png',dpi=220);plt.close(fig)
# Image self-energy is an analytical symmetric-valence control, not atomistic PMF.
z=np.linspace(.2,2,160);epsw,epspe,T=78.3,2.3,298.15;kb=R/6.02214076e23
lb=E**2/(4*np.pi*EPS0*epsw*kb*T)*1e9;eta=(epsw-epspe)/(epsw+epspe);dg=lb*eta/(4*z)
with (O/'dielectric_image_control.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.writer(f);w.writerow(['z_nm','dG_plus_kBT','dG_minus_kBT','K_plus','K_minus']);w.writerows(zip(z,dg,dg,np.exp(-dg),np.exp(-dg)))
fig,ax=plt.subplots(figsize=(6.4,4.1),constrained_layout=True);ax.plot(z,dg,color='#237e8a',label='Both monovalent signs (identical q^2)');ax.set(xlabel='Distance from planar dielectric boundary / nm',ylabel='Image self-energy / kBT',title='Dielectric contrast alone has no charge-sign preference');ax.legend(frameon=False);fig.text(.1,.01,'Assumed eps_water=78.3, eps_PE=2.3; z << fibre radius; no specific affinity.',fontsize=8);fig.savefig(O/'dielectric_symmetric_control.png',dpi=220);plt.close(fig)
(O/'plot_provenance.json').write_text(json.dumps({'surface_data_sha256':hashlib.sha256((O/'surface_scenarios.csv').read_bytes()).hexdigest(),'radial_data_sha256':hashlib.sha256((O/'surface_radial_profiles.csv').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'image_self_energy_scope':'analytical planar local approximation, hypothetical dielectric constants, exactly sign symmetric; not measured PE PMF','bjerrum_length_nm':lb,'dielectric_factor':eta},indent=2))
print('2 PNG, analytical control CSV written')
