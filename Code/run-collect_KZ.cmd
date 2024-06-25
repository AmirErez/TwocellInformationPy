#!/bin/bash
#SBATCH -o logs/collect_%A.out
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH -t 23:59:00
#SBATCH --mem=512000

python -u collect_KZ_Schlogl.py -i ../Data/Raw/Schlogl_nc1000_KZ_prot1_rate4 --max 1000000
#python -u collect_KZ_Schlogl.py -i ../Data/Raw/Schlogl_nc1000_KZ_prot2_rate4 --max 1000000
#python -u collect_KZ_Schlogl.py -i ../Data/Raw/Schlogl_nc1000_KZ_prot3_rate4 --max 1000000

#python -u collect_KZ_Landau.py -i ../Data/Raw/Landau_nc1000_KZ_prot1_rate5/ --max 500000
#python -u collect_KZ_Landau.py -i ../Data/Raw/Landau_nc1000_KZ_prot2_rate5/ --max 500000
#python -u collect_KZ_Landau.py -i ../Data/Raw/Landau_nc1000_KZ_prot3_rate5/ --max 500000

