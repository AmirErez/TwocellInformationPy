"""
Plots MI of the Landau calculation
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from newfig import newfig
#

collected_files = ('../Data/Collected/collected_Schlogl_nc1000_steadystate_sweeph.tsv', '../Data/Collected/collected_Landau_nc1000_steadystate_sweeph.tsv')
fig, ax = newfig(width_cm=4.3, height_cm=4.3, font_size=9)

for collected_file in collected_files:
    df = pd.read_csv(collected_file, sep="\t")
    # Plot
    df = df.sort_values(by="hx")
    ax.set_xlabel(r"$H$")
    ax.set_ylabel(r"$I(X;Y)$")
    simmethod = df['simulation_method'].unique()[0]
    # simnc = df['nc'].unique()[0]
    H = df["hx"] + df["hy"]
    ax.plot(H, df["MI"], "-", label=f"{simmethod}")

ax.set_xlim([-0.1, 0.1])
ax.set_ylim([0, 3])
ax.set_xticks([-0.1, 0, 0.1])
ax.set_yticks([0, 2])
# ax.set_xlabel(r'$H$')
# ax.set_ylabel(r'$I(X;Y)$')
# plt.legend()
# plt.plot(df['mean_hx'], df['corr'], '.-', label='KZ')
ax.grid()
ax.set_position([0.35, 0.35, 0.6, 0.55])

analytic_df = pd.read_csv('../Mathematica/Fig1_LNA_Landau_Same_Pts.csv', header=None, sep=",", names=['H', 'I'])
analytic_df['I'] = analytic_df['I'] / np.log(2)
ax.plot(analytic_df['H'], analytic_df['I'], 'k:', label='Gauss', linewidth=2)

analytic_df = pd.read_csv('../Mathematica/Fig1_LNA_Schlogl_Same_Pts.csv', header=None, sep=",", names=['H', 'I'])
analytic_df['I'] = analytic_df['I'] / np.log(2)
ax.plot(analytic_df['H'], analytic_df['I'], 'k:', label='Gauss', linewidth=2)


# Landau_numeric = pd.read_csv('../Data/Collected/collected_Landau_nc1000_analytic.tsv', sep="\t")
# Landau_numeric['MI'] = Landau_numeric['MI']/np.log(2)
# plt.plot(Landau_numeric['hx']+Landau_numeric['hy'], Landau_numeric['MI'], 'k:', label='Analytic')

plt.savefig('../Figures/fig_MI_nc_1000_sweeph.png', dpi=600)
plt.savefig('../Figures/fig_MI_nc_1000_sweeph.svg')
plt.show()


print('Finished')

