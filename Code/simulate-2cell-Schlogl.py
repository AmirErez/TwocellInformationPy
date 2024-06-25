from utils import Schlogl, read_json, rates_from_Schlogl, Rates
import time
import numpy as np
import pickle
import argparse
from scipy import sparse


def simulate_schlogl_2cell(x0, n_ising, m_ising, setup_data):

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

    Schlogl_n = Schlogl(n_ising)
    Schlogl_m = Schlogl(m_ising)
    rates = rates_from_Schlogl(Schlogl_n, Schlogl_m)

    Pn = np.zeros((1, x0[0] * 10))[0]
    Pm = np.zeros((1, x0[0] * 10))[0]
    Pnm = np.zeros((x0[0] * 6, x0[0] * 6))

    stoich_matrix = np.array([[1, 0], [-1, 0], [-1, 1], [1, -1], [0, 1], [0, -1]])
    delta_t = 0.01

    is_burnin = True
    if setup_data['simulation_type'] != 'static':
        raise ValueError('Wrong simulation type. This is a "static" simulation only.')
    while cur_t < setup_data['MAX_ITER']:
        rxn_count += 1
        if rxn_count % 100000 == 0:
            elapsed_realtime = time.time() - st
            print(f'... Done {rxn_count}, simulation time {cur_t:.2f}, real time {elapsed_realtime:.2f}')

        if is_burnin and cur_t >= setup_data['BURNIN']:
            is_burnin = False
            print('finished burnin')
            cur_t = 0
            rxn_count = 1

        if not is_burnin and cur_t-t[-1] >= delta_t and setup_data['store_arrays']:
                Xn.append(prevX[0])
                Xm.append(prevX[1])
                hx.append(n_ising.get_h(cur_t))
                hy.append(m_ising.get_h(cur_t))
                theta_x.append(n_ising.get_theta(cur_t))
                theta_y.append(m_ising.get_theta(cur_t))
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
        Pn[prevX[0]] += tau
        Pm[prevX[1]] += tau
        Pnm[prevX[0]][prevX[1]] += tau
        prevX += stoich_matrix[mu]

    Pn = Pn / Pn.sum()
    Pm = Pm / Pm.sum()
    Pnm = sparse.csr_matrix(Pnm)
    Pnm = Pnm / Pnm.sum()
    print(f'Finished in {time.time() - st} seconds')

    # If store_arrays is true, return all the arrays, otherwise only return the first four
    if setup_data['store_arrays']:
        out = {'Pn': Pn, 'Pm': Pm, 'Pnm': Pnm,
               't': np.array(t), 'n': np.array(Xn), 'm': np.array(Xm),
               'hx': np.array(hx), 'hy': np.array(hy),
               'theta_x': np.array(theta_x),  'theta_y': np.array(theta_y),
               'rxn_count': rxn_count, 'params': setup_data}
    else:
        out = {'Pn': Pn, 'Pm': Pm, 'Pnm': Pnm,
               'rxn_count': rxn_count, 'params': setup_data}
    return out


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-o', type=str, required=True, help='Output file')
    PARSER.add_argument('-i', type=str, required=True, help='Setup file')
    s = PARSER.parse_args()
    out_file = s.o
    setup_file = s.i

    print(f'Input params: {setup_file} ; Output file: {out_file}')

    ising_x, ising_y, setup_data = read_json(setup_file)

    ret = simulate_schlogl_2cell(x0=[ising_x.nc, ising_y.nc],
                                 n_ising=ising_x, m_ising=ising_y,
                                 setup_data=setup_data)

    with open(out_file, 'wb') as f:
        pickle.dump(ret, f)

    print('Done')
