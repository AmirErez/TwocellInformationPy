"""
Mutually useful functions
"""

import numpy as np
from dataclasses import dataclass
from sklearn.metrics import mutual_info_score

class Ising:
    def __init__(self, nc, theta, h, g):

        self.nc = nc  # k2p/(3*k2m) ## Effective number of particles, sets the noise in the Landau theory
        self.theta = theta  # 3*k1m*k2m/(k2p**2)-1   ## -0.1 , 0 , 0.1
        self.h = h  # 9*k1p*k2m**2/(k2p)**3 - 3*k1m*k2m/k2p**2 + 2/3  ## h =0
        self.g = g  # 3*gamma*k2m/k2p**2
        self.get_h = None  # For dynamics
        self.get_theta = None  # For dynamics

    def update(self, cur_t):
        if self.get_h is not None and self.get_theta is not None:
            self.h = self.get_h(cur_t)
            self.theta = self.get_theta(cur_t)
        else:
            raise ValueError("get_h and get_theta not defined")

    # k1p,k1m , k2p , k2m,gamma
    # def change_h(self, Hi, Hf, tau_d, t):
    #     if self.g != 0:
    #         self.h = (Hi+(Hf-Hi)*t/tau_d)/self.g


class Schlogl:
    def __init__(self, ising):

        self.s = 3*(ising.nc-1)
        self.k2 = 3*ising.theta*ising.nc**2 + 3*(ising.nc-1)*(ising.nc-2)+1
        self.k = self.k2**0.5
        self.a = ((3*ising.theta+3*ising.h+1) * ising.nc ** 3 - 6*(ising.nc-1))/(self.k**2)
        self.N = ising.nc*10
        self.g = ising.g


@dataclass
class Rates:
    kn1p: float
    kn2p: float
    kn1m: float
    kn2m: float
    gamma_nm: float
    gamma_mn: float
    km1p: float
    km2p: float
    km1m: float
    km2m: float


def rates_from_Schlogl(Schlogl_n, Schlogl_m):
    kn1m = 1
    km1m = kn1m / (Schlogl_n.k ** 2) * (Schlogl_m.k ** 2)
    kn1p = Schlogl_n.a * kn1m
    km1p = Schlogl_m.a * km1m
    kn2m = kn1m / Schlogl_n.k ** 2
    km2m = km1m / Schlogl_m.k ** 2
    kn2p = kn2m * Schlogl_n.s
    km2p = km2m * Schlogl_m.s
    gamma_mn = Schlogl_m.g * km1m / 3 * ((Schlogl_m.s + 3) / Schlogl_m.k) ** 2
    gamma_nm = Schlogl_n.g * kn1m / 3 * ((Schlogl_n.s + 3) / Schlogl_n.k) ** 2
    return Rates(kn1m=kn1m, km1m=km1m, kn1p=kn1p, km1p=km1p,
                 kn2m=kn2m, km2m=km2m, kn2p=kn2p, km2p=km2p,
                 gamma_mn=gamma_mn, gamma_nm=gamma_nm)


def read_json(setup_file):
    import json
    with open(setup_file, "r") as in_file:
        setup_data = json.load(in_file)
    ising_x = Ising(setup_data['nc'], setup_data['theta_x'], setup_data['hx'],
                    setup_data['gx'])  # nc , theta_x , hx , gx
    ising_y = Ising(setup_data['nc'], setup_data['theta_y'], setup_data['hy'],
                    setup_data['gy'])  # nc , theta_y , hy , gy
    if setup_data['simulation_type'] == 'KZ':
        ising_x.get_h = eval(setup_data['get_hx'])
        ising_x.get_theta = eval(setup_data['get_theta_x'])
        ising_y.get_h = eval(setup_data['get_hy'])
        ising_y.get_theta = eval(setup_data['get_theta_y'])
    elif setup_data['simulation_type'] == 'static':
        ising_x.get_h = lambda t: ising_x.h
        ising_x.get_theta = lambda t: ising_x.theta
        ising_y.get_h = lambda t: ising_y.h
        ising_y.get_theta = lambda t: ising_y.theta

    print(f"Ising x; nc={setup_data['nc']}, theta={setup_data['theta_x']}, hx={setup_data['hx']}, gx={setup_data['gx']}")
    print(f"Ising y; nc={setup_data['nc']}, theta={setup_data['theta_y']}, hx={setup_data['hy']}, gx={setup_data['gy']}")
    return ising_x, ising_y, setup_data


def calculate_entropy(prob_dist):
    """
        Calculate the entropy of a probability distribution
        Check if prob_dist is a sparse CSR matrix
        """
    import numpy as np
    from scipy.sparse import isspmatrix_csr

    if isspmatrix_csr(prob_dist):
        non_zero_data = prob_dist.data
    else:
        non_zero_data = prob_dist.flatten()[prob_dist.flatten() > 0]

    # calculate entropy
    non_zero_data = non_zero_data / sum(non_zero_data)
    entropy = -np.sum(non_zero_data * np.log2(non_zero_data), where=non_zero_data > 0)
    return entropy


def calculate_MI(Px, Py, Pxy):
    Sx = calculate_entropy(Px)
    Sy = calculate_entropy(Py)
    Sxy = calculate_entropy(Pxy)
    MI = Sx + Sy - Sxy
    return MI


def mean_val(P):
    v = np.arange(len(P))
    return np.sum(v*P)

def var_val(P):
    v = np.arange(len(P))
    return np.sum((v**2)*P)-np.sum(v*P)**2


def calculate_MI_scipy_bits(x, y):
    return mutual_info_score(x, y) / np.log(2)

def calculate_MI_from_events(x, y):
    if len(x) != len(y):
        raise ValueError("Arrays x and y must have the same length.")

    if not (x.size and y.size):
        raise ValueError("Arrays x and y must not be empty.")

    x = np.rint(np.asarray(x))
    y = np.rint(np.asarray(y))
    mx = np.max([np.max(x), np.max(y)])
    edges = np.arange(0, mx + 1)  # +1 because the upper edge is exclusive

    hist, xedges, yedges = np.histogram2d(x, y, bins=[edges, edges], density=True)
    Pxy = hist
    # get Px and Py
    Px = Pxy.sum(axis=1)
    Py = Pxy.sum(axis=0)
    return calculate_MI(Px, Py, Pxy)

