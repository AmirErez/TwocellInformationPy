"""
Simulate the two-cell system as a Landau system at equilibrium according to Eq. 2-3 in
"Cell-to-Cell Information at a Feedback-Induced Bifurcation Point"
Erez et al., PHYSICAL REVIEW LETTERS 125, 048103 (2020)

Written by Amir Erez, July 2023
"""
import time
import numpy as np
from scipy import sparse
import utils
import pickle
# from utils import Ising


def calc_energy(m, n_ising, m_ising):
    mx = m[0]
    my = m[1]
    Ex = -n_ising.h * mx + 1 / 2 * n_ising.theta * (mx ** 2) + 1 / 12 * mx ** 4
    Ey = -m_ising.h * my + 1 / 2 * m_ising.theta * (my ** 2) + 1 / 12 * my ** 4  # onsite potential energy
    # generalized Ginzburg kinetic energy:
    Eg = -1 / 2 * (n_ising.g + m_ising.g) * mx * my + 1 / 2 * n_ising.g * (mx ** 2) + 1 / 2 * m_ising.g * (my ** 2)
    return Ex + Ey + Eg


def simulate_landau_2cell(m0, n_ising, m_ising, setup_data):
    print('Simulating mx and my in TDGL dynamics')
    st = time.time()
    prevM = np.array(m0)
    prev_energy = calc_energy(prevM, n_ising, m_ising)

    mx = []
    my = []  # x[1] / m_ising.nc - 1

    energy = []
    Pn = np.zeros((1, n_ising.nc * 3))[0]
    Pm = np.zeros((1, m_ising.nc * 3))[0]
    Pnm = np.zeros((n_ising.nc * 3, m_ising.nc * 3))

    is_burnin = True
    cnt = 0
    n_accept = 0
    kBT = 1 / 2 * np.sqrt((1 + n_ising.theta) * (1 + m_ising.theta) / (n_ising.nc*m_ising.nc))  # Analytic expression
    kBT = kBT / 10

    # Below justified from analytics
    # var_noise_0 = (1 + n_ising.theta) * (n_ising.nc)
    # var_noise_1 = (1 + m_ising.theta) * (m_ising.nc)
    var_noise = 2*kBT  # should be D=2Gamma kBT = 2kBT
    print('CHANGED HERE TOO')
    var_noise = var_noise * 50

    while cnt < setup_data['MAX_ITER']:
        cnt += 1
        if cnt % 100000 == 0:
            print('... Done {}'.format(cnt))

        if is_burnin and cnt >= setup_data['BURNIN']:
            is_burnin = False
            print('Finished burnin with {:.2f}% percent accepted'.format(n_accept / cnt * 100))
            cnt = 0
            n_accept = 0
        if not is_burnin:
            # n_ising.chang_h(0.1, 0.0, 2 * 1e3, changT)  ## Specifies protocol
            # m_ising.chang_h(0, -0.1, 2 * 1e3, changT)
            if setup_data['store_arrays']:
                mx.append(prevM[0])
                my.append(prevM[1])
                energy.append(prev_energy)

        change = np.random.normal(0, np.sqrt(var_noise))
        if np.random.rand() < 0.5:
            # change = np.random.normal(0, np.sqrt(var_noise_0))
            newM = np.asarray([prevM[0] + change, prevM[1]])
        else:
            # change = np.random.normal(0, np.sqrt(var_noise_1))
            newM = np.asarray([prevM[0], prevM[1] + change])

        new_energy = calc_energy(newM, n_ising, m_ising)
        if np.exp((prev_energy - new_energy) / kBT) > np.random.rand():
            prevM = newM.copy()
            prev_energy = new_energy
            n_accept += 1

        if not is_burnin:
            x_coord = int(np.round((prevM[0]+1)*n_ising.nc))
            y_coord = int(np.round((prevM[1]+1)*m_ising.nc))
            Pn[x_coord] += 1
            Pm[y_coord] += 1
            Pnm[x_coord][y_coord] += 1

    Pn = Pn / Pn.sum()
    Pm = Pm / Pm.sum()
    Pnm = Pnm / Pnm.sum()
    Pnm = sparse.csr_matrix(Pnm)

    print('Finished simulation with {:.2f}% percent accepted'.format(n_accept / cnt * 100))
    print('Elapsed time {}'.format(time.time() - st))

    if setup_data['store_arrays']:
        out = {'Pn': Pn, 'Pm': Pm, 'Pnm': Pnm,
               'mx': np.array(mx), 'my': np.array(my),
               # 'hx': np.array(hx), 'hy': np.array(hy),
               'E': np.array(energy), 'n_accept': n_accept,
               'cnt': cnt, 'params': setup_data}
    else:
        out = {'Pn': Pn, 'Pm': Pm, 'Pnm': Pnm,
               'n_accept': n_accept, 'cnt': cnt, 'params': setup_data}
    return out


if __name__ == "__main__":
    import argparse

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-o', type=str, required=True, help='Output file')
    PARSER.add_argument('-i', type=str, required=True, help='Setup file')
    s = PARSER.parse_args()
    out_file = s.o
    setup_file = s.i

    print('Input params: {} ; Output file: {}'.format(setup_file, out_file))
    ising_x, ising_y, setup_data = utils.read_json(setup_file)

    ret = simulate_landau_2cell(m0=[0, 0], n_ising=ising_x, m_ising=ising_y, setup_data=setup_data)

    pickle.dump(ret, open(out_file, 'wb'))
    print('Done')
