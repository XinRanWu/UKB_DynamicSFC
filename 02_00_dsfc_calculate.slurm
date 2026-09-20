#!/bin/bash

#SBATCH --job-name=dynamic_sfc
#SBATCH --partition=DCU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=250G
#SBATCH --output=logs/dynamic_sfc_%j.out
#SBATCH --error=logs/dynamic_sfc_%j.err

# -------------------------------------------------------------------------
# Example SLURM submission script for SFC/DSFC calculation
#
# This configuration reflects the original HPC environment used in the study.
# Partition names, MATLAB modules, memory requirements, and CPU allocations
# should be adapted to the user's computing environment.
# -------------------------------------------------------------------------

set -e

hostname

module purge
module load matlab/R2018b

export TZ=Asia/Shanghai
ulimit -s unlimited

mkdir -p logs

matlab -nodisplay -nosplash -nodesktop \
    -r "try, compute_sfc_dsfc_glasser; catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
