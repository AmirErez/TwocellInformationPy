import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import utils
import pickle
import glob
import os


def init_alldata():
    return {'all_E':       [],
            'all_mx':      [],
            'all_my':      [],
            'all_hx':      [],
            'all_hy':      [],
            'all_theta_x': [],
            'all_theta_y': [] }


def corrtime(M):
    n_realizations = M.shape[0]
    n_steps = M.shape[1]
    tau_b = 1000 # Batch time, assuming X.shape[1] >> tau_b >> tau_corr
    n_batches = int(np.floor(n_realizations/tau_b))
    n_realizations = n_batches*tau_b
    tau_corrs = np.zeros(n_steps)
    batch_means = np.zeros(n_batches)
    print(f'Calculating correlation time over {n_realizations} realizations using batch length {tau_b}')
    for i in range(n_steps):
        vartot = np.var(M[:n_realizations,i])
        for bb in range(n_batches):
            batch_means[bb] = np.mean(M[bb*tau_b:(bb*tau_b+i), i])
        var_b = np.var(batch_means)
        tau_corrs[i] = tau_b * var_b / 2 / vartot
        if i % 1000 == 0:
            print(f'{i+1}/{n_steps}')
    return tau_corrs


def corr(X, Y):
    cr = np.nanmean(X * Y, axis=1) - np.nanmean(X, axis=1) * np.nanmean(Y, axis=1)
    return cr / (np.nanstd(X, axis=1) * np.nanstd(Y, axis=1))


def MI(X, Y):
    MI = []
    n_steps = X.shape[0]
    print('Calculating MI')
    for i in range(n_steps):
        MI.append(utils.calculate_MI_scipy_bits(np.rint(X[i, :]), np.rint(Y[i, :])))
        # MI.append(utils.calculate_MI_from_events(X[i, :], Y[i, :]))
        if i % 1000 == 0:
            print(f'{i+1}/{n_steps}')
    return np.asarray(MI)


def run(indir, maxnum):
    i = 0
    glb = glob.glob(indir + "/out_*.pkl")
    n_runs = len(glb)
    if n_runs == 0:
        raise ValueError(f'Directory {indir} has no input files!')
    n_steps = None

    alldata = init_alldata()
    cnt = 0

    saved_params = None
    for gg, file_name in enumerate(glb):
        if cnt > maxnum:
            print(f'Reached {cnt} realizations, stopping')
            break
        print(f'Reading {file_name} {gg+1}/{n_runs}')
        with open(file_name, 'rb') as fh:
            manydata = pickle.load(fh)
            for rep, data in enumerate(manydata):
                cnt += 1
                print(f'   Replicate {rep+1}')
                if saved_params is None:
                    saved_params = data['params']
                alldata['all_mx'].append(data['mx'])
                alldata['all_my'].append(data['my'])
                alldata['all_hx'].append(data['hx'])
                alldata['all_hy'].append(data['hy'])
                alldata['all_theta_x'].append(data['theta_x'])
                alldata['all_theta_y'].append(data['theta_y'])
                alldata['all_E'].append(data['E'])

    alldata['params'] = saved_params
    all_X =  (np.asarray(alldata['all_mx']).T + 1) * alldata['params']['nc']
    all_Y =  (np.asarray(alldata['all_my']).T + 1) * alldata['params']['nc']

    out = {'mean_hx': np.nanmean(np.asarray(alldata['all_hx']), axis=0),
           'mean_hy': np.nanmean(np.asarray(alldata['all_hy']), axis=0),
           'mean_theta_x': np.nanmean(np.asarray(alldata['all_theta_x']), axis=0),
           'mean_theta_y': np.nanmean(np.asarray(alldata['all_theta_y']), axis=0),
           'corrtime_M': corrtime((np.asarray(alldata['all_mx'])+np.asarray(alldata['all_my']))/2),
           'corr': corr(all_X, all_Y),
           'MI': MI(all_X, all_Y),
           'mean_E': np.nanmean(np.asarray(alldata['all_E']), axis=0),
           'mean_mx': np.nanmean(np.asarray(alldata['all_mx']), axis=0),
           'mean_my': np.nanmean(np.asarray(alldata['all_my']), axis=0),
           'std_mx': np.nanstd(np.asarray(alldata['all_mx']), axis=0),
           'std_my': np.nanstd(np.asarray(alldata['all_my']), axis=0)}
    out['H'] = alldata['params']['gx'] * (out['mean_hx']+out['mean_hy']) + \
        0.5*(out['mean_hx']*out['mean_theta_y'] + out['mean_hy']*out['mean_theta_x'])
    out['T'] = out['mean_theta_x']*out['mean_theta_y'] + alldata['params']['gx']*(out['mean_theta_x']+out['mean_theta_y'])
    del all_X
    del all_Y
    return pd.DataFrame(out)


if __name__ == "__main__":
    import argparse
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-i', type=str, required=True, help='Input directory')
    PARSER.add_argument('--max', type=int, required=False, default=100000, help='Input directory')
    s = PARSER.parse_args()
    dirpath = s.i
    # Remove trailing '/' if present
    if dirpath[-1] == '/':
        dirpath = dirpath[:-1]
    # Check the directory exists
    if not os.path.isdir(dirpath):
        raise ValueError(f'{dirpath} is not a directory.')
    dirname = os.path.basename(dirpath)
    outfile = os.path.join(dirpath, f'collected_{dirname}.tsv')
    print(f'Input dir: {dirname} ; Output file: {outfile}')
    df = run(indir=dirpath, maxnum=s.max)
    df.to_csv(outfile, index=False, sep='\t', float_format='%.6f')
    print('Done')
