#!/bin/bash -l
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cull0084@umn.edu
#SBATCH --job-name QH_D0XXXX.one_wags.slurm
#SBATCH -o slurm_logs/%j.zebra_zebref.prep_wags.out
#SBATCH -e slurm_logs/%j.zebra_zebref.prep_wags.err
#SBATCH -A fried255

set -e

conda activate snake73
cd $SLURM_SUBMIT_DIR

REF_DIR=/home/fried255/cull0084/.wags

snakemake -s prep_wags.smk \
    --use-singularity \
    --singularity-args "-B $PWD,$REF_DIR" \
    --profile slurm.go_wags/ \
    --configfile zebref_config.yaml

