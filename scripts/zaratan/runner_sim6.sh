#!/bin/bash
#SBATCH -t 2000
#SBATCH -c 1
#SBATCH --mem-per-cpu=500000
#SBATCH --account=fagan-lab-cmns
. ~/.bashrc
module load r
Rscript zar_tda_sim6.R
~
