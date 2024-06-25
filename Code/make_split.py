#!/usr/bin/env python
import json
import os
import numpy as np

# General params
config = {"nc": 1000,
          "hx": 0.0, "hy": 0.0,
          "theta_x": 0.0, "theta_y": 0.0,
          "gx": 1, "gy": 1,
          "BURNIN": 1000000,
          "MAX_ITER": 25000000,
          "simulation_type": "static",
          "simulation_method": "Landau",
          "store_arrays": False}

# ----------------------------------------------------------------------------------------
# Static, Landau simulation, scan theta
#config["simulation_method"] = "Landau"
#config['nc'] = 1000
#config['simulation_type'] = 'static'
#config["simulation_method"] = "Landau"
#config['BURNIN'] = 1000000
#config['MAX_ITER'] = 50000000
#theta_xs = np.linspace(-0.05, 0.05, 41)
#theta_ys = np.linspace(-0.05, 0.05, 41)
#hxs = np.zeros(theta_xs.shape)
#hys = np.zeros(theta_ys.shape)
#output_directory = "../AEData/Raw/Landau_nc1000_steadystate_sweeptheta"  # Specify the output directory here
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Static, Landau simulation, scan h
#config['nc'] = 1000
#config['simulation_type'] = 'static'
#config["simulation_method"] = "Landau"
#config['BURNIN'] = 1000000
#config['MAX_ITER'] = 200000000
#hxs = np.linspace(-0.05, 0.05, 81)
#hys = np.linspace(-0.05, 0.05, 81)
#theta_xs = np.zeros(hxs.shape)
#theta_ys = np.zeros(hys.shape)
#output_directory = "../AEData/Raw/Landau_nc1000_steadystate_sweeph"  # Specify the output directory here
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Static, Schlogl simulation, scan theta
# theta_xs = np.linspace(-0.2, 0.2, 41)
# theta_ys = np.linspace(-0.2, 0.2, 41)
# hxs = np.zeros(theta_xs.shape)
# hys = np.zeros(theta_ys.shape)
# output_directory = "../AEData/Raw/Schlogl_nc1000_steadystate_sweeptheta"  # Specify the output directory here
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Static, Schlogl simulation, scan h
config['nc'] = 1000
config['simulation_method'] = 'Schlogl'
#config['BURNIN'] = 3000
#config['MAX_ITER'] = 40000
config['BURNIN'] = 5000
config['MAX_ITER'] = 5000
hxs = np.linspace(-0.05, 0.05, 81)
hys = np.linspace(-0.05, 0.05, 81)
theta_xs = np.zeros(hxs.shape)
theta_ys = np.zeros(hys.shape)
output_directory = "../AEData/Raw/Schlogl_nc1000_steadystate_sweeph_shortt"  # Specify the output directory here
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# KZ, Schlogl, scan h
#config['simulation_type'] = 'KZ'
#config['simulation_method'] = 'Schlogl'
#output_directory = "../AEData/Raw/Schlogl_nc1000_KZ_prot1"  # Specify the output directory here
## Setup ramp
#hi = 0.05
#hf = -0.05

#config['hxi'] = hi
#config['hyi'] = hi
#config['hxf'] = hf
#config['hyf'] = hf

#config['delta_t'] = 0.01
#config['t_ramp'] = 1000
#config['t_start'] = 100
#config['BURNIN'] = 1000
#t_start = config['t_start']
#t_ramp = config['t_ramp']
#config['MAX_ITER'] = config['t_ramp']+2*config['t_start']

#config['get_hx'] = f'lambda t: {hi} if t<{t_start} else ({hi} +({hf}-{hi})*(t-{t_start})/{t_ramp}) if (t-{t_start})<{t_ramp} else {hf}'
#config['get_hy'] = f'lambda t: {hi} if t<{t_start} else ({hi} +({hf}-{hi})*(t-{t_start})/{t_ramp}) if (t-{t_start})<{t_ramp} else {hf}'
#config['get_theta_x'] = f'lambda t: 0'
#config['get_theta_y'] = f'lambda t: 0'
# ----------------------------------------------------------------------------------------


# ########################################### CODE ###########################################
def make_split(output_dir, theta_xs, theta_ys, hxs, hys):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for idx, (theta_x, theta_y, hx, hy) in enumerate(zip(theta_xs, theta_ys, hxs, hys)):
        config["theta_x"] = theta_x
        config["theta_y"] = theta_y
        config["hx"] = hx
        config["hy"] = hy

        file_name = f"setup_{idx:03d}.json"
        file_path = os.path.join(output_dir, file_name)

        with open(file_path, "w") as file:
            json.dump(config, file, indent=4)
        print(f"Created file: {file_path}")


def make_split_kz(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    idx = 0
    file_name = f"setup_{idx:03d}.json"
    file_path = os.path.join(output_dir, file_name)
    with open(file_path, "w") as file:
        json.dump(config, file, indent=4)
    print(f"Created file: {file_path}")


if __name__ == "__main__":
    if config['simulation_type'] == 'KZ':
        make_split_kz(output_directory)
    else:
        make_split(output_directory, theta_xs, theta_ys, hxs, hys)
