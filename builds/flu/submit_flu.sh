#!/bin/sh
#
# Reserve 1 CPUs for this job
#
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=05:59:00

source ~/.bashrc
ml MAFFT

python $@
