"""
Plots MI of the Landau calculation for several beta values of the same params
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from newfig import newfig
import pickle
from utils import read_json
from calc_MI_Landau_from_energy import calc_energy
from scipy.sparse import csr_matrix
#

setup_file = '../Data/Raw/Landau_nc1000_steadystate_sweeph/setup_000.json'
ising_x, ising_y, setup_data = read_json(setup_file)
ising_x.h = 0
ising_y.h = 0
yind = ising_x.nc + 1

## For reading the MC simulation distributions - not needed since they perfectly match the numeric
fig, ax = newfig(width_cm=4.3, height_cm=4.3, font_size=9)
mc_betax10 = pickle.load(open(f'../Data/Collected/out_040__Landau_nc1000_steadystate_sweeph_betax10.pkl', 'rb'))
mc_betax1 = pickle.load(open(f'../Data/Collected/out_040__Landau_nc1000_steadystate_sweeph.pkl', 'rb'))
mc_betaover10 = pickle.load(open(f'../Data/Collected/out_040__Landau_nc1000_steadystate_sweeph_betaover10.pkl', 'rb'))
val_betax1 = mc_betax1['Pnm'][:, yind].toarray().squeeze()
val_betax10 = mc_betax10['Pnm'][:, yind].toarray().squeeze()
val_betaover10 = mc_betaover10['Pnm'][:, yind].toarray().squeeze()
m = np.asarray(np.arange(len(val_betax1))/ mc_betax1['params']['nc'] - 1)
# dm = m[1] - m[0]
val_betax1 = val_betax1 / np.trapz(val_betax1, m)
val_betax10 = val_betax10 / np.trapz(val_betax10, m)
val_betaover10 = val_betaover10 / np.trapz(val_betaover10, m)

# cmap = plt.cm.summer
# colors = cmap(np.linspace(-1, 1, 5))
colors = plt.cm.tab10(np.linspace(0, 1, 10))

plt.plot(m, val_betax10, label='$k_BT/10$', color=colors[0])
plt.plot(m, val_betax1, label='$k_BT$', color=colors[1])
plt.plot(m, val_betaover10, label='$10k_BT$', color=colors[4])

# Analytic calculation from numerical values
m = np.asarray(np.arange(ising_x.nc*2+1)/ ising_x.nc - 1)
E = calc_energy(m, m[yind], ising_x, ising_y)
beta = 2 * ising_x.nc
Pn_betax10 = np.exp(-beta * 10 * E)
Pn_betax10 = Pn_betax10 / np.trapz(Pn_betax10, m)
Pn_betax1 = np.exp(-beta * E)
Pn_betax1 = Pn_betax1 / np.trapz(Pn_betax1, m)
Pn_betaover10 = np.exp(-beta / 10 * E)
Pn_betaover10 = Pn_betaover10 / np.trapz(Pn_betaover10, m)
ax.plot(m, Pn_betax10, 'k:', label='Analytic')
ax.plot(m, Pn_betax1, 'k:', label='Analytic')
ax.plot(m, Pn_betaover10, 'k:', label='Analytic')

# ax.set_yscale('log')
ax.set_xlim([-0.1, 0.1])
plt.xlabel('$m$')
plt.ylabel('$P(m_X | m_Y=0 )$')
# plt.legend()
plt.tight_layout()

# Landau_numeric = pd.read_csv('../Data/Collected/collected_Landau_nc1000_analytic.tsv', sep="\t")
# Landau_numeric['MI'] = Landau_numeric['MI']/np.log(2)
# plt.plot(Landau_numeric['hx']+Landau_numeric['hy'], Landau_numeric['MI'], 'k:', label='Analytic')

# Landau_betax10['MI'] = Landau_betax10['MI']
# plt.plot(Landau_betax10['hx']+Landau_betax10['hy'], Landau_betax10['MI'], 'c:', label='kBT/10')

plt.savefig('../Figures/fig_compare_Pn_betas.png', dpi=600)
plt.savefig('../Figures/fig_compare_Pn_betas.svg')
plt.show()


# Landau_numeric = pd.read_csv('../Data/Collected/collected_Landau_nc1000_analytic.tsv', sep="\t")
# Landau_numeric['MI'] = Landau_numeric['MI']/np.log(2)
# plt.plot(Landau_numeric['hx']+Landau_numeric['hy'], Landau_numeric['MI'], 'k:', label='Analytic')

fig, ax = newfig(width_cm=4.3, height_cm=4.3, font_size=9)
# betas = np.array([0.1, 1, 10]) * ising_x.nc * 2
# numeric_dfs = []
# for beta in betas:
#     numeric_dfs.append(pd.read_csv(f'../Data/Collected/collected_Landau_nc1000_analytic_beta_{beta}.tsv', sep="\t"))
#     ax.plot(numeric_dfs[-1]['hx']+numeric_dfs[-1]['hy'], numeric_dfs[-1]['MI'], label=f'kBT/{beta/(2*ising_x.nc)}')

Landau_betaover2 = pd.read_csv('../Data/Collected/collected_Landau_nc1000_steadystate_sweeph_betaover2.tsv', sep="\t").sort_values(by="hx")
Landau_betaover3 = pd.read_csv('../Data/Collected/collected_Landau_nc1000_steadystate_sweeph_betaover3.tsv', sep="\t").sort_values(by="hx")
Landau_betaover10 = pd.read_csv('../Data/Collected/collected_Landau_nc1000_steadystate_sweeph_betaover10.tsv', sep="\t").sort_values(by="hx")
Landau_betax1 = pd.read_csv('../Data/Collected/collected_Landau_nc1000_steadystate_sweeph.tsv', sep="\t").sort_values(by="hx")
Landau_betax10 = pd.read_csv('../Data/Collected/collected_Landau_nc1000_steadystate_sweeph_betax10.tsv', sep="\t").sort_values(by="hx")

ax.set_xlim([-0.1, 0.1])
ax.set_ylim([0, 3])
ax.set_xticks([-0.1, 0, 0.1])
ax.set_yticks([0, 2])

plt.plot(Landau_betaover10['hx']+Landau_betaover10['hy'], Landau_betaover10['MI'], label='$10 k_BT_{MC}$', color=colors[4])
plt.plot(Landau_betaover3['hx']+Landau_betaover3['hy'], Landau_betaover3['MI'], label='$3k_BT_{MC}$', color=colors[3])
plt.plot(Landau_betaover2['hx']+Landau_betaover2['hy'], Landau_betaover2['MI'], label='$2 k_BT_{MC}$', color=colors[2])
plt.plot(Landau_betax1['hx']+Landau_betax1['hy'], Landau_betax1['MI'], label='$k_BT_{MC}$', color=colors[1])
plt.plot(Landau_betax10['hx']+Landau_betax10['hy'], Landau_betax10['MI'], label='$k_BT_{MC}/10$', color=colors[0])





analytic_df = pd.read_csv('../Mathematica/Fig1_LNA_Landau_Same_Pts.csv', header=None, sep=",", names=['H', 'I'])
analytic_df['I'] = analytic_df['I'] / np.log(2)
ax.plot(analytic_df['H'], analytic_df['I'], 'k:', label='Gauss', linewidth=2)

ax.set_xlim([-0.1, 0.1])
ax.set_ylim([0, 4])
ax.set_xticks([-0.1, 0, 0.1])
ax.set_yticks([0, 2, 4])
ax.set_xlabel(r'$H$')
ax.set_ylabel(r'$I(X;Y)$')
#plt.legend()
ax.grid()
ax.set_position([0.35, 0.35, 0.6, 0.55])

plt.savefig('../Figures/fig_MI_nc_1000_sweeph_compare_betas.png', dpi=600)
plt.legend()
# plt.tight_layout()
plt.savefig('../Figures/fig_MI_nc_1000_sweeph_compare_betas.svg')
plt.show()

# Define colormap as copper and use it for the following plots
fig, ax = newfig(width_cm=4.3, height_cm=4.3, font_size=9)
# Set colormap
plt.plot(Landau_betaover10['hx']+Landau_betaover10['hy'], Landau_betaover10['varmx'], label='$10 k_BT_{MC}$', color=colors[4])
plt.plot(Landau_betaover3['hx']+Landau_betaover3['hy'], Landau_betaover3['varmx'], label='$3k_BT_{MC}$', color=colors[3])
plt.plot(Landau_betaover2['hx']+Landau_betaover2['hy'], Landau_betaover2['varmx'], label='$2 k_BT_{MC}$', color=colors[2])
plt.plot(Landau_betax1['hx']+Landau_betax1['hy'], Landau_betax1['varmx'], label='$k_BT_{MC}$', color=colors[1])
plt.plot(Landau_betax10['hx']+Landau_betax10['hy'], Landau_betax10['varmx'], label='$k_BT_{MC}/10$', color=colors[0])



df_schlogl = pd.read_csv('../Data/Collected/collected_Schlogl_nc1000_steadystate_sweeph.tsv', sep="\t").sort_values(by="hx")
plt.plot(df_schlogl['hx']+df_schlogl['hy'], df_schlogl['varmx'], 'k:', linewidth=2, label='Schlogl')
# plt.yscale('log')
plt.ylabel(r'$Var(m_X)$')
plt.xlim([-0.1, 0.1])
ax.set_xticks([-0.1, 0, 0.1])
ax.set_xlabel(r'$H$')
# Show legend with latex interpreter for the labels
# plt.legend()
plt.tight_layout()
plt.savefig('../Figures/fig_varmx_compare_betas.png', dpi=600)
plt.legend()
plt.savefig('../Figures/fig_varmx_compare_betas.svg')
plt.show()

print('Finished')

