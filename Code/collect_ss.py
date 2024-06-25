#!python
import glob
import os
import pandas as pd
import argparse
import utils
import pickle


def process_file(file_path):
    """
    Process an individual file
    """
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    MI = utils.calculate_MI(Px=data['Pn'], Py=data['Pm'], Pxy=data['Pnm'])
    out = data['params'].copy()
    out['mx'] = (utils.mean_val(P=data['Pn'])-data['params']['nc'])/data['params']['nc']
    out['my'] = (utils.mean_val(P=data['Pm'])-data['params']['nc'])/data['params']['nc']
    out['varmx'] = utils.var_val(P=data['Pn'])/(data['params']['nc']**2)
    out['varmy'] = utils.var_val(P=data['Pm'])/(data['params']['nc']**2)
    out['MI'] = MI
    out['file'] = file_path
    return out


def process_files_in_directory(directory_path):
    """
    Process all .npz files in the specified directory
    """
    file_paths = glob.glob(os.path.join(directory_path, "*.pkl"))
    n = len(file_paths)

    print(f'Reading {n} files from {directory_path}')
    data_list = []
    for file_path in file_paths:
        print(f'Processing file {file_path}')
        data_list.append(process_file(file_path))

    return data_list


def save_data(data_list, output_file):
    """
    Save calculated data
    """
    df = pd.DataFrame(data_list)
    df.to_csv(output_file, sep='\t', index=False)


def main(directory_path, output_file):
    """
    Main function
    """
    data_list = process_files_in_directory(directory_path)
    save_data(data_list, output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process .npz files and save data to a .tsv file.')
    parser.add_argument('directory_path', type=str, help='Path to directory containing .npz files.')

    args = parser.parse_args()
    dirpath = args.directory_path
    # Remove trailing '/' if present
    if dirpath[-1] == '/':
        dirpath = dirpath[:-1]
    # Check the directory exists
    if not os.path.isdir(dirpath):
        raise ValueError(f'{dirpath} is not a directory.')
    dirname = os.path.basename(dirpath)
    outfile = os.path.join(dirpath, f'collected_{dirname}.tsv')
    main(args.directory_path, outfile)
