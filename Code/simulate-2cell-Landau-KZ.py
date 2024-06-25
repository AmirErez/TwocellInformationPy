"""
Simulate a KZ ramp of the two-cell system as a Landau system at according to Eq. 2-3 in
"Cell-to-Cell Information at a Feedback-Induced Bifurcation Point"
Erez et al., PHYSICAL REVIEW LETTERS 125, 048103 (2020)

Written by Amir Erez, July 2023
"""
import time
import numpy as np
import pickle
from utils import read_json
import os

def calc_energy(m, n_ising, m_ising):
    mx = m[0]
    my = m[1]
    Ex = -n_ising.h * mx + 1 / 2 * n_ising.theta * (mx ** 2) + 1 / 12 * mx ** 4
    Ey = -m_ising.h * my + 1 / 2 * m_ising.theta * (my ** 2) + 1 / 12 * my ** 4  # onsite potential energy
    # generalized Ginzburg kinetic energy:
    Eg = -1 / 2 * (n_ising.g + m_ising.g) * mx * my + 1 / 2 * n_ising.g * (mx ** 2) + 1 / 2 * m_ising.g * (my ** 2)
    return Ex + Ey + Eg


def simulate_landau_2cell_kz(m0, n_ising, m_ising, setup_data):
    print('Simulating mx and my in TDGL dynamics')
    st = time.time()
    prevM = np.array(m0)
    prev_energy = calc_energy(prevM, n_ising, m_ising)

    mx = [0]
    my = [0]
    hx = [n_ising.h]
    hy = [m_ising.h]
    energy = [prev_energy]
    theta_x = [n_ising.theta]
    theta_y = [m_ising.theta]

    is_burnin = True
    cnt = 0
    n_accept = 0
    kBT = 1 / 2 * np.sqrt((1 + n_ising.theta) * (1 + m_ising.theta) / (n_ising.nc*m_ising.nc))  # Analytic expression

    while is_burnin or cnt < setup_data['MAX_ITER']:
        cnt += 1
        if cnt % 10000 == 0:
            elapsed_realtime = time.time() - st
            print(f'... Done {cnt}, real time {elapsed_realtime:.2f}')

        if is_burnin and cnt >= setup_data['BURNIN']:
            is_burnin = False
            print('Finished burnin with {:.2f}% percent accepted'.format(n_accept / cnt * 100))
            cnt = 0
            n_accept = 0

        if not is_burnin:
            n_ising.update(cnt)
            m_ising.update(cnt)
            hx.append(n_ising.h)
            hy.append(m_ising.h)
            theta_x.append(n_ising.theta)
            theta_y.append(m_ising.theta)
            mx.append(prevM[0])
            my.append(prevM[1])
            energy.append(prev_energy)

        # Below justified from analytics
        # var_noise_0 = (1 + n_ising.theta) * (n_ising.nc)
        # var_noise_1 = (1 + m_ising.theta) * (m_ising.nc)
        var_noise = 2*kBT  # should be D=2Gamma kBT = 2kBT
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

    print(f'Finished in {time.time() - st} seconds')

    out = {'mx': np.array(mx), 'my': np.array(my),
           'hx': np.array(hx), 'hy': np.array(hy),
           'theta_x': np.array(theta_x), 'theta_y': np.array(theta_y),
           'E': np.array(energy), 'n_accept': n_accept, 'cnt': cnt}

    print('Finished simulation with {:.2f}% percent accepted'.format(n_accept / cnt * 100))
    print('Elapsed time {}'.format(time.time() - st))
    return out


if __name__ == "__main__":
    import argparse

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-o', type=str, required=True, help='Output file')
    PARSER.add_argument('-i', type=str, required=True, help='Setup file')
    PARSER.add_argument('-r', type=int, required=False, default=1, help='Number of replicates (default 1)')
    s = PARSER.parse_args()
    out_file = s.o
    setup_file = s.i
    n_replicates = s.r

    print(f'Input params: {setup_file} ; Output file: {out_file} ; Replicates: {n_replicates}')

    ising_x, ising_y, setup = read_json(setup_file)
    if setup['simulation_type'] != 'KZ':
        raise ValueError('Wrong simulation type. This is a "KZ" simulation only.')
    print(setup)

    all_ret = []
    for ii in range(n_replicates):
       print('----------------------------------------------')
       print(f'Staring replicate {ii+1}')
       unique_seed = (int(time.time() * 1000) + os.getpid()) % (2**32)
       np.random.seed(unique_seed)
       ret = simulate_landau_2cell_kz(m0=[0, 0], n_ising=ising_x, m_ising=ising_y, setup_data=setup)
       ret['params'] = setup
       all_ret.append(ret)

    with open(out_file, 'wb') as f:
        pickle.dump(all_ret, f)
    print('Done')

