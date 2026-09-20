#!/bin/bash

#SBATCH --job-name=dynamic_sfc
#SBATCH --partition=DCU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=250G
#SBATCH --output=logs/dynamic_sfc_%j.out
#SBATCH --error=logs/dynamic_sfc_%j.err

# ============================================================================
# Example SLURM submission script for the Dynamic SFC analysis
#
# Original environment:
#   Linux HPC cluster
#   SLURM
#   MATLAB R2018b
#
# The partition name, MATLAB module, CPU allocation, memory allocation,
# and filesystem paths are cluster-specific and should be adapted to the
# user's computing environment.
# ============================================================================

set -e

echo "------------------------------------------------------------"
echo "Dynamic SFC analysis"
echo "------------------------------------------------------------"

echo "Host:"
hostname

echo "SLURM job ID:"
echo "${SLURM_JOB_ID}"

echo "Allocated CPUs:"
echo "${SLURM_CPUS_PER_TASK}"

echo "Start time:"
date

# ----------------------------------------------------------------------------
# MATLAB environment
# ----------------------------------------------------------------------------

module purge
module load matlab/R2018b

# Unlimited stack size was used in the original HPC environment.
ulimit -s unlimited

# ----------------------------------------------------------------------------
# Analysis parameters
# ----------------------------------------------------------------------------
#
# Window length can be passed from the command line:
#
#   sbatch submit_sfc_dsfc.slurm 30
#   sbatch submit_sfc_dsfc.slurm 40
#   sbatch submit_sfc_dsfc.slurm 50
#
# Default = 30 TR.
# ----------------------------------------------------------------------------

WINDOW_LENGTH=${1:-30}

echo "Window length: ${WINDOW_LENGTH} TR"

mkdir -p logs

# ----------------------------------------------------------------------------
# Run MATLAB
# ----------------------------------------------------------------------------

matlab \
    -nodisplay \
    -nosplash \
    -nodesktop \
    -r "try, \
        WINDOW_LENGTH=${WINDOW_LENGTH}; \
        WINDOW_STEP=5; \
        N_WORKERS=${SLURM_CPUS_PER_TASK}; \
        run('02_compute_sfc_dsfc_glasser.m'); \
        exit(0); \
        catch ME, \
        disp(getReport(ME,'extended')); \
        exit(1); \
        end"

echo "End time:"
date
