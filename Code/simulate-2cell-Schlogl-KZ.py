from utils import Schlogl, read_json, rates_from_Schlogl
import time
import numpy as np
import pickle
import argparse
import random
import os

def simulate_schlogl_2cell_kz(x0, n_ising, m_ising, t_burnin, t_max, delta_t):
    st = time.time()
    prevX = x0.copy()
    cur_t = 0
    rxn_count = 0
    Xn = [x0[0]]
    Xm = [x0[1]]
    t = [0]
    hx = [n_ising.h]
    hy = [m_ising.h]
    theta_x = [n_ising.theta]
    theta_y = [m_ising.theta]

    rates = rates_from_Schlogl(Schlogl(n_ising), Schlogl(m_ising))
    stoich_matrix = np.array([[1, 0], [-1, 0], [-1, 1], [1, -1], [0, 1], [0, -1]])

    is_burnin = True
    while is_burnin or cur_t < t_max:
        rxn_count += 1
        if rxn_count % 100000 == 0:
            elapsed_realtime = time.time() - st
            print(f'... Done {rxn_count}, simulation time {cur_t:.2f}, real time {elapsed_realtime:.2f}')

        if is_burnin and cur_t >= t_burnin:
            is_burnin = False
            print('finished burnin')
            cur_t = 0
            rxn_count = 1

        if not is_burnin:
            n_ising.update(cur_t)
            m_ising.update(cur_t)
            rates = rates_from_Schlogl(Schlogl(n_ising), Schlogl(m_ising))
            if cur_t-t[-1] >= delta_t:
                Xn.append(prevX[0])
                Xm.append(prevX[1])
                hx.append(n_ising.h)
                hy.append(m_ising.h)
                theta_x.append(n_ising.theta)
                theta_y.append(m_ising.theta)
                t.append(cur_t)

        a = np.array([rates.kn1p + rates.kn2p * prevX[0] * (prevX[0] - 1),
                      rates.kn1m * prevX[0] + rates.kn2m * prevX[0] * (prevX[0] - 1) * (prevX[0] - 2),
                      rates.gamma_nm * prevX[0], rates.gamma_mn * prevX[1],
                      rates.km1p + rates.km2p * prevX[1] * (prevX[1] - 1),
                      rates.km1m * prevX[1] + rates.km2m * prevX[1] * (prevX[1] - 1) * (prevX[1] - 2)])

        a0 = np.sum(a)
        r = np.random.rand(1, 2)[0]
        tau = -np.log(r[0]) / a0

        mu = 0
        sa = a[0]
        r0 = r[1] * a0
        while sa < r0:
            mu += 1
            sa += a[mu]

        cur_t += tau
        prevX += stoich_matrix[mu]

    print(f'Finished in {time.time() - st} seconds')

    out = {'t': np.array(t), 'x': np.array(Xn), 'y': np.array(Xm),
           'hx': np.array(hx), 'hy': np.array(hy),
           'theta_x': np.array(theta_x),  'theta_y': np.array(theta_y),
           'rxn_count': rxn_count, 't_burnin': t_burnin, 't_max': t_max}
    return out


if __name__ == "__main__":
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
       unique_seed = (int(time.time() * 1000) + os.getpid()) % 2**32
       np.random.seed(unique_seed)
       ising_x.update(0)
       ising_y.update(0)
       x0 = [random.randint(0,ising_x.nc*2), random.randint(0,ising_y.nc*2)]
       ret = simulate_schlogl_2cell_kz(x0=x0, n_ising=ising_x, m_ising=ising_y,
                                       t_burnin=setup['BURNIN'], t_max=setup['MAX_ITER'], delta_t=setup['delta_t'])
       ret['params'] = setup
       all_ret.append(ret)
    with open(out_file, 'wb') as f:
        pickle.dump(all_ret, f)

    print('Done')
