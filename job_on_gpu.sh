#!/usr/bin/sh

## this must be run from directory where run.py exists.
## --workdir is not used in this file.

#SBATCH --job-name=recon_distmap_using_vae
#SBATCH --qos=csqos
##SBATCH --workdir=/scratch/akabir4/project_dir
#SBATCH --output=/scratch/akabir4/csbw_20_project/outputs/logs/recon_distmap_using_vae_log_1-%N-%j.output
#SBATCH --error=/scratch/akabir4/csbw_20_project/outputs/logs/recon_distmap_using_vae_log_1-%N-%j.error
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --gres=gpu:1
#SBATCH --partition=gpuq
#SBATCH --mem=32G

## python full_run_recon_distmap_using_ae.py
python full_run_recon_distmap_using_vae.py