"""Calculat the mutual information between X and Y in the 2-cell Landay model using the energy values
"""
import numpy as np
import utils
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import xlogy


def calc_energy(mx, my, n_ising, m_ising):
    Ex = -n_ising.h * mx + 0.5 * n_ising.theta * (mx ** 2) + (mx ** 4) / 12
    Ey = -m_ising.h * my + 0.5 * m_ising.theta * (my ** 2) + (my ** 4) / 12
    Eg = -0.5 * (n_ising.g + m_ising.g) * mx * my + 0.5 * n_ising.g * (mx ** 2) + 0.5 * m_ising.g * (my ** 2)
    return Ex + Ey + Eg


def calc_MI_Landau(ising_x, ising_y, beta):
    """
    Calculate the mutual information between X and Y in the 2-cell Landay model using the energy values
    """
    # Calculate the probability distribution of the energy values

    X = np.arange(ising_x.nc * 2) + 1
    Y = np.arange(ising_y.nc * 2) + 1
    mx = X / ising_x.nc - 1
    my = Y / ising_y.nc - 1

    mx_grid, my_grid = np.meshgrid(mx, my)
    E = calc_energy(mx_grid, my_grid, ising_x, ising_y)
    # kBT = 1 / 2 * np.sqrt((1 + ising_x.theta) * (1 + ising_y.theta) / (ising_x.nc * ising_y.nc))
    # beta = 1/kBT
    Pnm = np.exp(-beta * E)
    # Z is the double sum over Pnm
    Z = np.trapz(np.trapz(Pnm, mx, axis=0), my)
    Pnm = Pnm / Z
    # Plot heatmap of Pnm, with m on the x-axis and n on the y-axis
    # sns.heatmap(Pnm)
    # plt.xlim([0, 2*ising_x.nc])
    # plt.ylim([0, 2*ising_x.nc])
    # plt.show()
    # Calculate the entropy of the energy values
    Pn = np.trapz(Pnm, mx, axis=1)
    Pm = np.trapz(Pnm, my, axis=0)
    S_n = -np.trapz(xlogy(Pn, Pn), mx)
    S_m = -np.trapz(xlogy(Pm, Pm), my)
    S_nm = -np.trapz(np.trapz(xlogy(Pnm, Pnm), mx), my)
    # Calculate the mutual information between X and Y
    MI = S_n + S_m - S_nm
    return MI / np.log(2)


def run(ising_x, ising_y, beta):
    hs = np.linspace(-0.05, 0.05, 21)
    MI = np.zeros(len(hs))
    for ii, hh in enumerate(hs):
        print(f'Calculating MI for h={hh}')
        ising_x.h = hh
        ising_y.h = hh
        MI[ii] = calc_MI_Landau(ising_x, ising_y, beta)
    return hs, MI


if __name__ == "__main__":
    # Load the data
    setup_file = '../AEData/Raw/Landau_nc1000_steadystate_sweeph/setup_000.json'
    ising_x, ising_y, setup_data = utils.read_json(setup_file)
    betas = np.array([0.1, 1, 10]) * ising_x.nc * 2
    for beta in betas:
        print(f'Calculating MI for beta={beta}')
        hs, MI = run(ising_x, ising_y, beta)
        df = pd.DataFrame({'hx': hs, 'hy': hs,'MI': MI})
        df.to_csv(f'../AEData/Collected/collected_Landau_nc1000_analytic_beta_{beta}.tsv', sep="\t", index=False)
        plt.plot(hs, MI)
    plt.show()
    print('Finished')