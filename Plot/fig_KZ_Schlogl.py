import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from newfig import newfig



files = ('../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate1.tsv',
         '../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate2.tsv',
         '../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate3.tsv',
         '../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate4.tsv',
         '../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate1_return.tsv',
         '../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate2_return.tsv',
         '../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate3_return.tsv',
         '../Data/Collected/collected_Schlogl_ncvar_KZ_prot1_rate4_return.tsv')

tau_rs = (2000, 1000, 500, 250, 2000, 1000, 500, 250)

cols = np.asarray((1, 2, 3, 4, 1, 2, 3, 4))/4

# Define colormap for tau_rs
cmap = mpl.colormaps.get_cmap('PuRd')
fig, ax = newfig(width_cm=4.3, height_cm=4.3, font_size=9)

for ff, file in enumerate(files):
    df = pd.read_csv(file, sep="\t")
    tau_r = tau_rs[ff]
    H_g1 = (df['mean_hx'] + df['mean_hy']) + 0.5 * (df['mean_hx'] * df['mean_theta_y'] + df['mean_hy'] * df['mean_theta_x'])
    m = (df['mean_mx'] + df['mean_my']) / 2
    ax.plot(H_g1[2:], m[2:], '-', label=rf'$\tau_d={tau_r}$', color=cmap(cols[ff]), markersize=1.5)

ax.set_xlim([-0.1, 0.1])
ax.set_ylim([-0.6, 0.6])
ax.set_xticks([-0.1, 0, 0.1])
ax.set_yticks([-0.6, 0, 0.6])
ax.set_xlabel(r'$H$')
ax.set_ylabel(r'$M$')
# plt.legend()
ax.set_position([0.35, 0.35, 0.6, 0.55])
plt.savefig('../Figures/fig_KZ_Schlogl__ncvar.png', dpi=600)
plt.savefig('../Figures/fig_KZ_Schlogl__ncvar.svg')
plt.show()

delta = 3
beta = 0.5
nu = 0.5
z = 2

fig, ax = newfig(width_cm=4.3, height_cm=4.3, font_size=9)
for ff, file in enumerate(files):
    df = pd.read_csv(file, sep="\t")
    tau_r = tau_rs[ff]
    H_g1 = (df['mean_hx'] + df['mean_hy']) + 0.5 * (df['mean_hx'] * df['mean_theta_y'] + df['mean_hy'] * df['mean_theta_x'])
    m = (df['mean_mx'] + df['mean_my']) / 2
    scale_H = H_g1 * (tau_r ** (beta * delta / (nu * z + beta * delta)))
    scale_m = m * (tau_r ** (beta / (nu * z + beta * delta)))
    ax.plot(scale_H[10:], scale_m[10:], '-', label=fr'$\tau_d={tau_r}$', color=cmap(cols[ff]), markersize=1.5)

ax.set_xlim([-10, 10])
ax.set_ylim([-3, 3])
ax.set_xticks([-10, 0, 10])
ax.set_yticks([-2, 0, 2])
ax.set_xlabel(r'$H\,\tau_d^{\beta \delta/(\nu z + \beta \delta)}$')
ax.set_ylabel(r'$M\,\tau_d^{\beta/(\nu z + \beta\delta)}$')
# ax.grid()
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles[:4], labels[:4], loc='upper left', ncol=1, frameon=False)
# ax.grid()
ax.set_position([0.35, 0.35, 0.6, 0.55])
plt.savefig('../Figures/fig_KZ_Schlogl__ncvar_collapse.png', dpi=600)
plt.savefig('../Figures/fig_KZ_Schlogl__ncvar_collapse.svg')
plt.show()


print('Done')