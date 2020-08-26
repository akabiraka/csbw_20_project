#!/usr/bin/sh

#SBATCH --array=0-2050
#SBATCH --job-name=relaxProtocol_on_testset_translated_fragments
#SBATCH --qos=csqos
#SBATCH --output=/scratch/akabir4/csbw_20_project/outputs/logs/relaxProtocol_on_testset_fragments-%N-%A-%a.output
#SBATCH --error=/scratch/akabir4/csbw_20_project/outputs/logs/relaxProtocol_on_testset_fragments-%N-%A-%a.error

#SBATCH --partition=all-HiPri
#SBATCH --cpus-per-task=2
#SBATCH --mem=4000MB

# the array task is set in the environment variable $SLURM_ARRAY_TASK_ID in python you
# can scrape it with ID = int(os.environ["SLURM_ARRAY_TASK_ID"])

python evaluations/test_RelaxProtocol.py