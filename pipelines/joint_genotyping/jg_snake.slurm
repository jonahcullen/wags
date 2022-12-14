#!/bin/bash -l
#SBATCH -t 60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=cull0084@umn.edu 
#SBATCH --job-name jg_snake
#SBATCH -o %j.jg_snake.out
#SBATCH -e %j.jg_snake.err
#SBATCH -p small,amdsmall

set -e

conda activate snake718
cd $SLURM_SUBMIT_DIR


snakemake -s many_wags.smk \
    --use-singularity \
    --profile slurm.go_wags/ \
    --configfile canfam4_config.yaml \
    --keep-going

