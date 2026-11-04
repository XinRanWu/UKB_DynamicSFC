#!/usr/bin/env bash
#SBATCH --job-name=ukb_sfc_dsfc
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/coupling_%j.out

set -euo pipefail
mkdir -p logs
matlab -batch "run('matlab/run_primary_pipeline.m')"
