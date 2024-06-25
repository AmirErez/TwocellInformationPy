#!/bin/bash
#SBATCH --array=1-81
#SBATCH -o logs/schlogl_%A_%a.out
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH -t 23:59:00
#SBATCH --mem=8000

OFFSET=-1
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)

#outdir=../Data/Raw/Schlogl_nc1000_steadystate_sweeph
outdir=../Data/Raw/Schlogl_nc1000_steadystate_sweeph_shortt

files=( $( ls $outdir/setup*.json ) )
infile=${files[$LINE_NUM]}
# in bash, extract the suffix after setup_ and before .json
filenum=${infile#*setup_}
filenum=${filenum%.json}
outfile=$outdir/out_${filenum}.pkl
if [ -e $outfile ]; then
   echo "File $outfile exists. Quitting"
    exit 1
fi

echo "-----------------------------------------------------------"
echo $infile $outdir
python -u ./simulate-2cell-Schlogl.py -i $infile -o $outfile
if [ $? -ne 0 ]; then
    echo "Error. Exiting."
    exit 1
fi

echo "Done. Written $outfile."

