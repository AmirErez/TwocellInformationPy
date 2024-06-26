#!/bin/bash
#SBATCH --array=2-2000
#SBATCH -o logs/landau_%A_%a.out
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH -t 23:59:00
#SBATCH --mem=4000

if [ -z "$1" ]; then
    echo "Warning! No offset argument supplied. Setting to -1"
    OFFSET=-1
else
    OFFSET=$1
    echo "The offset is $OFFSET"
fi
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)

#outdir=../AEData/Raw/Landau_ncvar_KZ_prot1_rate5_return-betax10
outdir=../AEData/Raw/Landau_ncvar_KZ_prot1_rate5_return-betaover3

infile=$outdir/setup_000.json
outfile=$outdir/out_${LINE_NUM}.pkl
if [ -e $outfile ]; then
   echo "File $outfile exists. Quitting"
    exit 1
fi

echo "-----------------------------------------------------------"
echo $infile $outdir
python -u ./simulate-2cell-Landau-KZ-betaover.py -i $infile -o $outfile -r 100
if [ $? -ne 0 ]; then
    echo "Error. Exiting."
    exit 1
fi

echo "Done. Written $outfile."

