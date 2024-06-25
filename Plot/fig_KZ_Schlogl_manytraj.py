import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from newfig import newfig

files = ['../Data/Collected/collected_Schlogl_nc1000_KZ_prot1_rate4.tsv',
         '../Data/Collected/collected_Schlogl_nc1000_KZ_prot2_rate4.tsv',
         '../Data/Collected/collected_Schlogl_nc1000_KZ_prot3_rate4.tsv']
ss_df = pd.read_csv('../Data/Collected/collected_Schlogl_nc1000_steadystate_sweeph.tsv', sep="\t")

ss_df = ss_df.sort_values(by="hx")

tau_rs = [250, 250, 250]
cols = np.asarray([1, 2, 3])/3

colors = ['c', 'y', 'm']
linestyles = ['-', '--', ':']
fig, ax = newfig(width_cm=4.3, height_cm=4.3, font_size=9)
for ff, file in enumerate(files):
    df = pd.read_csv(file, sep="\t")
    tau_r = tau_rs[ff]
    H_g1 = (df['mean_hx'] + df['mean_hy']) + 0.5 * (df['mean_hx'] * df['mean_theta_y'] + df['mean_hy'] * df['mean_theta_x'])
    m = (df['mean_mx'] + df['mean_my']) / 2
    ax.plot(H_g1[1:], df['MI'][1:], label=f'P{ff + 1}', color=colors[ff], linestyle=linestyles[ff], markersize=1)
    # ax.plot(H_g1, df['MI'], '.', label=f'P{ff+1}', color=cmap(cols[ff]), alpha=0.2, markersize=1)
ax.plot(ss_df['hx']+ss_df['hy'], ss_df['MI']/np.log(2), '--', label=f'SS', color='k', markersize=2)

ax.set_xlabel(r'$H$')
ax.set_ylabel(r'$I(X;Y)$')
# ax.set_xlim([-0.1, 0.1])
ax.set_xlim([-0.1, 0.1])
# ax.set_xticks([-0.1, 0, 0.1])
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', ncol=1, frameon=False)
# plt.grid()
ax.set_ylim([0, 3])
ax.set_yticks([0, 1, 2, 3])
ax.set_position([0.25, 0.25, 0.7, 0.65])
plt.savefig('../Figures/fig_KZ_Schlogl_MI_manytraj.png', dpi=600)
plt.savefig('../Figures/fig_KZ_Schlogl_MI_manytraj.svg')
plt.show()
# plt.close()


# Plot H vs t
# cols = np.asarray([1, 1, 2, 2, 3])/3
# cmap = mpl.colormaps.get_cmap('Dark2')
fig, axs = newfig(width_cm=8.6, height_cm=4.3, font_size=9,nrows=1, ncols=3)
for ff, file in enumerate(files):
    ax = axs[ff]
    df = pd.read_csv(file, sep="\t")
    ax.plot(df['t'], df['mean_hx'], '-', label=r'$h_x$', color=plt.cm.Dark2.colors[0], markersize=1)
    ax.plot(df['t'], df['mean_hy'], '--', label=r'$h_y$', color=plt.cm.Dark2.colors[0], markersize=1)
    ax.plot(df['t'], df['mean_theta_x'], '-', label=r'$\theta_x$', color=plt.cm.Dark2.colors[1], markersize=1)
    ax.plot(df['t'], df['mean_theta_y'], '--', label=r'$\theta_x$', color=plt.cm.Dark2.colors[1], markersize=1)
    ax.plot(df['t'], df['H'], '-k', label=f'H', markersize=1)
    ax.set_xlabel(r'Time, $t$')
    ax.set_title(f'Protocol {ff+1}', color=colors[ff])
    if ff > 0:
        ax.set_yticklabels([])
    # ax.set_xticks([0, 225, 450])
    # ax.legend()
    # ax.grid()
plt.tight_layout()
fig.subplots_adjust(top=0.73, bottom=0.25, wspace=0.08)
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=5, frameon=False)

plt.savefig(f'../Figures/fig_KZ_Schlogl_params.png', dpi=600)
plt.savefig(f'../Figures/fig_KZ_Schlogl_params.svg')
plt.show()

print('Done')